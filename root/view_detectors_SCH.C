// view_detectors_SCH.C
// - BH2(15), BVH_U(22), SCH(64) 2D 히트 맵 확인용
// - SCH 히스토그램은 (1) 데이터 스캔 기반 자동 범위, 히트가 없으면
//   (2) DCGeomParam(E72)의 기하(센터=420, segW=11.5, segH=400, nSeg=64)에
//   여유 마진을 얹은 안전 범위를 사용

#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TBox.h"
#include "TLine.h"
#include "TParticle.h"
#include "TStyle.h"
#include "TInterpreter.h"
#include "TSystem.h"
#include <vector>
#include <algorithm>
#include <iostream>

// ---------- overlay helpers ----------
static void draw_bh2_overlay(double center_x) {
  const double seg_w = 14.0, seg_h = 100.0; const int n_segs = 15;
  const double total_w = n_segs*seg_w, start_x = center_x - total_w/2.0;
  auto *outline = new TBox(start_x, -seg_h/2., start_x+total_w, seg_h/2.);
  outline->SetFillStyle(0); outline->SetLineColor(kRed); outline->SetLineWidth(2); outline->Draw("SAME");
  for (int i=1;i<n_segs;++i){ const double x=start_x+i*seg_w; auto *l=new TLine(x,-seg_h/2.,x,seg_h/2.); l->SetLineColor(kRed); l->SetLineStyle(2); l->Draw("SAME");}
}

static void draw_bvh_u_overlay(double center_x) { // 22 seg
  const double seg_w = 10.0, seg_h = 140.0; const int n_segs = 22;
  const double total_w = n_segs*seg_w, start_x = center_x - total_w/2.0;
  auto *outline = new TBox(start_x, -seg_h/2., start_x+total_w, seg_h/2.);
  outline->SetFillStyle(0); outline->SetLineColor(kRed); outline->SetLineWidth(2); outline->Draw("SAME");
  for (int i=1;i<n_segs;++i){ const double x=start_x+i*seg_w; auto *l=new TLine(x,-seg_h/2.,x,seg_h/2.); l->SetLineColor(kRed); l->SetLineStyle(2); l->Draw("SAME");}
}

static void draw_sch_overlay(double center_x, double seg_w=11.5, double seg_h=400.0, int n_segs=64) {
  const double total_w = n_segs*seg_w, start_x = center_x - total_w/2.0;
  auto *outline = new TBox(start_x, -seg_h/2., start_x+total_w, seg_h/2.);
  outline->SetFillStyle(0); outline->SetLineColor(kRed); outline->SetLineWidth(2); outline->Draw("SAME");
  for (int i=1;i<n_segs;++i){ const double x=start_x+i*seg_w; auto *l=new TLine(x,-seg_h/2.,x,seg_h/2.); l->SetLineColor(kRed); l->SetLineStyle(2); l->Draw("SAME");}
}

// ---------- main ----------
void view_detectors_SCH(const char* fname="../E45_with_SCH.root",
                        double schCenterX=420.0,   // DCGeomParam 106행
                        double schSegW=11.5,       // DetSize SchSeg x
                        double schSegH=400.0,      // DetSize SchSeg y(전체)
                        int    schNSeg=64)         // NumOfSegSCH
{
  // TParticle 딕셔너리 준비(필수)
  gSystem->Load("libPhysics");
  gInterpreter->GenerateDictionary("vector<TParticle>","TParticle.h;vector");

  gStyle->SetOptStat(1111);
  gStyle->SetPalette(kBird);

  // 파일/트리
  TFile *f = TFile::Open(fname,"READ");
  if (!f || f->IsZombie()) { Error("view_detectors_SCH","cannot open file: %s", fname); return; }
  auto *tree = (TTree*)f->Get("g4hyptpc");
  if (!tree) { Error("view_detectors_SCH","tree g4hyptpc not found"); return; }

  // 브랜치 연결
  std::vector<TParticle> *bh2_hits=nullptr, *bvh_u_hits=nullptr, *sch_hits=nullptr;
  if (tree->GetBranch("BH2"))   tree->SetBranchAddress("BH2",&bh2_hits);
  else { Error("view_detectors_SCH","branch 'BH2' not found"); return; }

  if (tree->GetBranch("BVH_U")) tree->SetBranchAddress("BVH_U",&bvh_u_hits);
  else if (tree->GetBranch("BVH")) { tree->SetBranchAddress("BVH",&bvh_u_hits); Warning("view_detectors_SCH","using 'BVH' for BVH_U"); }

  if (tree->GetBranch("SCH"))   tree->SetBranchAddress("SCH",&sch_hits);
  else { Warning("view_detectors_SCH","branch 'SCH' not found; skipping SCH plot"); }

  // ---------- 1차 스캔: 실제 SCH 히트 범위 ----------
  double xmin=  1e9, xmax=-1e9, ymin=  1e9, ymax=-1e9;
  Long64_t N = tree->GetEntries();
  Long64_t Nscan = std::min<Long64_t>(N, 100000);
  Long64_t nSCHfilled = 0;

  for (Long64_t i=0; i<Nscan; ++i) {
    tree->GetEntry(i);
    if (!sch_hits) continue;
    for (const auto& h : *sch_hits) {
      const double x = h.Vx();
      const double y = h.Vy();
      if (!(x==x) || !(y==y)) continue; // NaN 방지
      xmin = std::min(xmin, x);
      xmax = std::max(xmax, x);
      ymin = std::min(ymin, y);
      ymax = std::max(ymax, y);
      ++nSCHfilled;
    }
  }

  // 데이터가 없으면 기하학적 예상 범위로 넉넉히
  if (nSCHfilled == 0) {
    const double x_half = (schSegW*schNSeg)/2.0;     // 64*11.5/2=368
    const double y_half = schSegH/2.0;               // 200
    const double xMargin = 100.0;                    // 여유
    const double yMargin = 100.0;
    xmin = schCenterX - x_half - xMargin;            // 420-368-100 = -48
    xmax = schCenterX + x_half + xMargin;            // 420+368+100 = 888
    ymin = -y_half - yMargin;                         // -300
    ymax =  y_half + yMargin;                         //  300
    Warning("view_detectors_SCH","No SCH hits in first %lld events. Using safe geometric ranges.", Nscan);
  }

  // 여유 마진 (자동 범위일 때도 살짝 확장)
  const double mx = 0.08*(xmax-xmin + 1.0);
  const double my = 0.08*(ymax-ymin + 1.0);

  // ---------- 히스토그램 ----------
  TH2D *h_hit_bh2   = new TH2D("h_hit_bh2",   "BH2 Hit Distribution;X [mm];Y [mm]",
                               250, -100, 150, 200, -100, 100);
  TH2D *h_hit_bvh_u = new TH2D("h_hit_bvh_u", "BVH_U Hit Distribution;X [mm];Y [mm]",
                               300, -100, 200, 200, -150, 150);

  // SCH는 넉넉 + 1~2mm/bin 수준
  const double xr = (xmax - xmin + 2*mx);
  const double yr = (ymax - ymin + 2*my);
  int nx = std::clamp((int)std::round(xr/1.5), 200, 2000);
  int ny = std::clamp((int)std::round(yr/1.5), 100, 1000);

  TH2D *h_hit_sch   = new TH2D("h_hit_sch",
                               "SCH Hit Distribution;X [mm];Y [mm]",
                               nx, xmin-mx, xmax+mx,
                               ny, ymin-my, ymax+my);

  // ---------- 본 채우기 ----------
  for (Long64_t i=0; i<N; ++i){
    tree->GetEntry(i);
    if (bh2_hits)   for (const auto& h:*bh2_hits)   h_hit_bh2->Fill(h.Vx(),h.Vy());
    if (bvh_u_hits) for (const auto& h:*bvh_u_hits) h_hit_bvh_u->Fill(h.Vx(),h.Vy());
    if (sch_hits)   for (const auto& h:*sch_hits)   h_hit_sch->Fill(h.Vx(),h.Vy());
  }

  // ---------- 그리기 ----------
  TCanvas *c = new TCanvas("c_hits","2D Hit Distributions with Overlays",1800,600);
  c->Divide(3,1);
  c->cd(1); h_hit_bh2->Draw("COLZ");   draw_bh2_overlay(35.0);   // DCGeomParam: BH2 X=35
  c->cd(2); h_hit_bvh_u->Draw("COLZ"); draw_bvh_u_overlay(40.0); // DCGeomParam: BVH_U X=40
  c->cd(3); h_hit_sch->Draw("COLZ");   draw_sch_overlay(schCenterX, schSegW, schSegH, schNSeg);

  c->SaveAs("hits_with_overlays_sch.png");

  // 디버그 출력
  printf("[SCH debug] scanned=%lld, hits=%lld, auto-range X:[%.1f, %.1f], Y:[%.1f, %.1f]\n",
         Nscan, nSCHfilled, xmin, xmax, ymin, ymax);
}

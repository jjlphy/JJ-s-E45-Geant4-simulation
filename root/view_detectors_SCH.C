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
#include <iostream>

// --- 오버레이 ---
void draw_bh2_overlay(double center_x) {
  const double seg_w = 14.0, seg_h = 100.0; const int n_segs = 15;
  const double total_w = n_segs*seg_w, start_x = center_x - total_w/2.0;
  auto *outline = new TBox(start_x, -seg_h/2., start_x+total_w, seg_h/2.);
  outline->SetFillStyle(0); outline->SetLineColor(kRed); outline->SetLineWidth(2); outline->Draw("SAME");
  for (int i=1;i<n_segs;++i){ double x=start_x+i*seg_w; auto *l=new TLine(x,-seg_h/2.,x,seg_h/2.); l->SetLineColor(kRed); l->SetLineStyle(2); l->Draw("SAME");}
}

void draw_bvh_u_overlay(double center_x) { // n=22
  const double seg_w = 10.0, seg_h = 140.0; const int n_segs = 22;
  const double total_w = n_segs*seg_w, start_x = center_x - total_w/2.0;
  auto *outline = new TBox(start_x, -seg_h/2., start_x+total_w, seg_h/2.);
  outline->SetFillStyle(0); outline->SetLineColor(kRed); outline->SetLineWidth(2); outline->Draw("SAME");
  for (int i=1;i<n_segs;++i){ double x=start_x+i*seg_w; auto *l=new TLine(x,-seg_h/2.,x,seg_h/2.); l->SetLineColor(kRed); l->SetLineStyle(2); l->Draw("SAME");}
}

// --- SCH 오버레이 (BVH_D 대체) ---
void draw_sch_overlay(double center_x) { // n=64, 11.5 x 400 mm
  const double seg_w = 11.5, seg_h = 400.0; const int n_segs = 64;
  const double total_w = n_segs*seg_w, start_x = center_x - total_w/2.0;
  auto *outline = new TBox(start_x, -seg_h/2., start_x+total_w, seg_h/2.);
  outline->SetFillStyle(0); outline->SetLineColor(kRed); outline->SetLineWidth(2); outline->Draw("SAME");
  for (int i=1;i<n_segs;++i){ double x=start_x+i*seg_w; auto *l=new TLine(x,-seg_h/2.,x,seg_h/2.); l->SetLineColor(kRed); l->SetLineStyle(2); l->Draw("SAME");}
}

void view_detectors_SCH() { // 함수 이름도 새로
  // 딕셔너리/라이브러리 로드
  gSystem->Load("libPhysics");
  gInterpreter->GenerateDictionary("vector<TParticle>","TParticle.h;vector");

  gStyle->SetOptStat(1111);
  gStyle->SetPalette(kBird);

  // ★ 입력 ROOT 파일명은 필요에 맞게 바꿔주세요
  TFile *f = TFile::Open("E45_SCH1.root","READ");
  if (!f || f->IsZombie()) { Error("view_detectors_SCH","cannot open file"); return; }
  auto *tree = (TTree*)f->Get("g4hyptpc");
  if (!tree) { Error("view_detectors_SCH","tree g4hyptpc not found"); return; }

  // 브랜치 목록(확인용)
  // tree->GetListOfBranches()->Print();

  std::vector<TParticle> *bh2_hits=nullptr, *bvh_u_hits=nullptr, *sch_hits=nullptr;

  // 브랜치 연결
  if (tree->GetBranch("BH2")) tree->SetBranchAddress("BH2",&bh2_hits);
  else { Error("view_detectors_SCH","branch 'BH2' not found"); return; }

  if (tree->GetBranch("BVH_U")) tree->SetBranchAddress("BVH_U",&bvh_u_hits);
  else if (tree->GetBranch("BVH")) { tree->SetBranchAddress("BVH",&bvh_u_hits); Warning("view_detectors_SCH","using 'BVH' for BVH_U"); }

  // ★ BVH_D → SCH
  if (tree->GetBranch("SCH")) tree->SetBranchAddress("SCH",&sch_hits);
  else { Warning("view_detectors_SCH","branch 'SCH' not found; skipping SCH plot"); }

  // 히스토그램 (SCH 범위는 중심 420 mm, 높이 ±200 mm 기준)
  TH2D *h_hit_bh2   = new TH2D("h_hit_bh2",   "BH2 Hit Distribution;X [mm];Y [mm]", 250, -100, 150, 200, -100, 100);
  TH2D *h_hit_bvh_u = new TH2D("h_hit_bvh_u", "BVH_U Hit Distribution;X [mm];Y [mm]", 300, -100, 200, 200, -150, 150);
  TH2D *h_hit_sch   = new TH2D("h_hit_sch",   "SCH Hit Distribution;X [mm];Y [mm]",  500, 100, 600, 200, -200, 200);

  const Long64_t n = tree->GetEntries();
  for (Long64_t i=0;i<n;++i){
    tree->GetEntry(i);
    if (bh2_hits)   for (const auto& h:*bh2_hits)   h_hit_bh2->Fill(h.Vx(),h.Vy());
    if (bvh_u_hits) for (const auto& h:*bvh_u_hits) h_hit_bvh_u->Fill(h.Vx(),h.Vy());
    if (sch_hits)   for (const auto& h:*sch_hits)   h_hit_sch->Fill(h.Vx(),h.Vy());
  }

  TCanvas *c = new TCanvas("c_hits","2D Hit Distributions with Overlays",1800,600);
  c->Divide(3,1);
  c->cd(1); h_hit_bh2->Draw("COLZ");   draw_bh2_overlay(35.0);    // BH2 중심 X = 35
  c->cd(2); h_hit_bvh_u->Draw("COLZ"); draw_bvh_u_overlay(40.0);  // BVH_U 중심 X = 40
  c->cd(3); h_hit_sch->Draw("COLZ");   draw_sch_overlay(420.0);   // SCH 중심 X = 420
  c->SaveAs("hits_with_overlays_sch.png");
}

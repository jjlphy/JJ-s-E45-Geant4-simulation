// -*- C++ -*-
// HTOF_2D_Hitpattern_10_28.C
// 2025-10-28 (for jaejin)
// 목적:
//   HTOF 34타일(0..33)을 사용해 2D 힛패턴 생성:
//     (A) π+  vs  π−
//     (B) p   vs  p̄
//     (C) π+  vs  p
//   그리고 파이(π± 합친) 1D 타일 점유도도 생성.
//   각 히스토그램은
//     1) BH2 seg 4–10 조건만 적용한 버전,
//     2) BH2 seg 4–10 && HTOF multiplicity≥2 버전
//   두 벌을 모두 그립니다.
//
// 가정/정책:
//   - 입력 TTree: <treename> (기본: "g4hyptpc")
//   - Branch:
//       vector<TParticle> BH2, HTOF   (필수)
//       vector<double>   BH2_edep, HTOF_edep (있으면 사용; 없으면 TParticle::Weight())
//   - BH2 seg ID 매핑: E72-like (x좌표→seg). (원 매크로와 동일)
//   - VALID HTOF hit: 타일별 ΣEdep ≥ threshold
//   - HTOF multiplicity: VALID 타일의 개수(입자종 구분 없이 union)
//   - 타겟 도달 여부는 별도 컷 없이, BH2 조건만 사용(요청사항 기준).
//
// 사용법(예):
//   root -l
//   .L HTOF_2D_Hitpattern_10_28.C+
//   HTOF_2D_Hitpattern_10_28("../rootfile/E45_Nov_piplusn_098.root","g4hyptpc", /*bh2_lo=*/4, /*bh2_hi=*/10, /*mipFrac=*/0.10, /*mipMeVperCm=*/2.0,/*BH2_thk_mm=*/5.0, /*HTOF_thk_mm=*/10.0,/*save=*/true, /*tag=*/"BH2_4_10");
//       "E45.root","g4hyptpc",
//       /*bh2_lo=*/4, /*bh2_hi=*/10,
//       /*mipFrac=*/0.10, /*mipMeVperCm=*/2.0,
//       /*BH2_thk_mm=*/5.0, /*HTOF_thk_mm=*/10.0,
//       /*save=*/true, /*tag=*/"BH2_4_10");
//
// 출력 파일(저장 옵션 true일 때):
//   HTOF_2D_pip_vs_pim_BH2_<tag>.png
//   HTOF_2D_pip_vs_pim_BH2_MPge2_<tag>.png
//   HTOF_2D_p_vs_pbar_BH2_<tag>.png
//   HTOF_2D_p_vs_pbar_BH2_MPge2_<tag>.png
//   HTOF_2D_pip_vs_p_BH2_<tag>.png
//   HTOF_2D_pip_vs_p_BH2_MPge2_<tag>.png
//   HTOF_1D_piOcc_BH2_<tag>.png
//   HTOF_1D_piOcc_BH2_MPge2_<tag>.png

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1I.h"
#include "TH2I.h"
#include "TSystem.h"
#include "TInterpreter.h"
#include "TParticle.h"
#include "TString.h"
#include "TLegend.h"
#include <vector>
#include <map>
#include <set>
#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>

namespace H2D_HP_10_28 {

// ----- dictionary for vector<TParticle> -----
struct DictGuard {
  DictGuard(){
    gSystem->Load("libPhysics");
    gInterpreter->GenerateDictionary("vector<TParticle>","TParticle.h;vector");
  }
};
static DictGuard _dg;

// ----- BH2 geometry mapping (E72-like) -----
static const int    kNBH2Seg   = 15;     // 0..14
static const double kBH2_x0    = -10.0;  // [mm]
static const double kBH2_pitch = 14.0;   // [mm]
static const double kBH2_total = kNBH2Seg * kBH2_pitch;

static inline int MapBH2_WorldToSeg(double x, double, double){
  const double xloc = x - kBH2_x0;
  const double xmin = -0.5*kBH2_total, xmax = 0.5*kBH2_total;
  if (xloc < xmin || xloc >= xmax) return -1;
  int idx = int(std::floor((xloc - xmin)/kBH2_pitch));
  if (idx<0) idx=0;
  if (idx>=kNBH2Seg) idx=kNBH2Seg-1;
  return idx;
}

static inline double MIPThrMeV(double frac,double dEdx,double thk_mm){
  return frac * dEdx * (thk_mm*0.1); // mm → cm
}

} // namespace H2D_HP_10_28


void HTOF_2D_Hitpattern_10_28(const char* filename="E45.root",
                               const char* treename="g4hyptpc",
                               int    bh2_lo=4, int bh2_hi=10,
                               double mipFrac=0.10,
                               double mipMeVperCm=2.0,
                               double BH2_thickness_mm=5.0,
                               double HTOF_thickness_mm=10.0,
                               bool   save=true,
                               const char* tag="BH2_4_10")
{
  using namespace H2D_HP_10_28;

  // ---------------- Open file & tree ----------------
  TFile* f=TFile::Open(filename,"READ");
  if(!f||f->IsZombie()){ std::cerr<<"[ERR] open "<<filename<<" failed\n"; return; }
  TTree* T=(TTree*)f->Get(treename);
  if(!T){ std::cerr<<"[ERR] tree "<<treename<<" not found\n"; return; }

  if(!T->GetBranch("BH2") || !T->GetBranch("HTOF")){
    std::cerr<<"[ERR] need BH2 & HTOF (vector<TParticle>) branches\n"; return;
  }

  // Optional edep branches
  const bool hasBH2edep  = (T->GetBranch("BH2_edep")  != nullptr);
  const bool hasHTOFedep = (T->GetBranch("HTOF_edep") != nullptr);

  // Set addresses
  std::vector<TParticle>* BH2=nullptr;   T->SetBranchAddress("BH2",&BH2);
  std::vector<TParticle>* HTOF=nullptr;  T->SetBranchAddress("HTOF",&HTOF);
  std::vector<double>* BH2_edep=nullptr;   if(hasBH2edep)  T->SetBranchAddress("BH2_edep",&BH2_edep);
  std::vector<double>* HTOF_edep=nullptr;  if(hasHTOFedep) T->SetBranchAddress("HTOF_edep",&HTOF_edep);

  // ---------------- Thresholds ----------------
  const double thrBH2  = MIPThrMeV(mipFrac,mipMeVperCm,BH2_thickness_mm);   // e.g. ~0.10 MeV
  const double thrHTOF = MIPThrMeV(mipFrac,mipMeVperCm,HTOF_thickness_mm);  // e.g. ~0.20 MeV

  std::cout<<std::fixed<<std::setprecision(3)
           <<"[INFO] BH2 thr="<<thrBH2<<" MeV, HTOF thr="<<thrHTOF<<" MeV | "
           <<"BH2 range = "<<bh2_lo<<"-"<<bh2_hi<<"\n";
  std::cout<<"[INFO] VALID HTOF tile = ΣEdep(tile) ≥ threshold; tile id = StatusCode(copy-no)\n";

  // ---------------- Histograms ----------------
  const int NT = 34; // HTOF tiles 0..33
  TH2I* h2_pip_pim_BH2        = new TH2I("h2_pip_pim_BH2",
    "HTOF 2D: #pi^{+} vs #pi^{-} | BH2 4-10;#pi^{+} tile;#pi^{-} tile", NT, -0.5, NT-0.5, NT, -0.5, NT-0.5);
  TH2I* h2_pip_pim_BH2_MPge2  = new TH2I("h2_pip_pim_BH2_MPge2",
    "HTOF 2D: #pi^{+} vs #pi^{-} | BH2 4-10 & MP#geq2;#pi^{+} tile;#pi^{-} tile", NT, -0.5, NT-0.5, NT, -0.5, NT-0.5);

  TH2I* h2_p_pbar_BH2         = new TH2I("h2_p_pbar_BH2",
    "HTOF 2D: p vs #bar{p} | BH2 4-10;p tile;#bar{p} tile", NT, -0.5, NT-0.5, NT, -0.5, NT-0.5);
  TH2I* h2_p_pbar_BH2_MPge2   = new TH2I("h2_p_pbar_BH2_MPge2",
    "HTOF 2D: p vs #bar{p} | BH2 4-10 & MP#geq2;p tile;#bar{p} tile", NT, -0.5, NT-0.5, NT, -0.5, NT-0.5);

  TH2I* h2_pip_p_BH2          = new TH2I("h2_pip_p_BH2",
    "HTOF 2D: #pi^{+} vs p | BH2 4-10;#pi^{+} tile;p tile", NT, -0.5, NT-0.5, NT, -0.5, NT-0.5);
  TH2I* h2_pip_p_BH2_MPge2    = new TH2I("h2_pip_p_BH2_MPge2",
    "HTOF 2D: #pi^{+} vs p | BH2 4-10 & MP#geq2;#pi^{+} tile;p tile", NT, -0.5, NT-0.5, NT, -0.5, NT-0.5);

  TH1I* h1_piOcc_BH2          = new TH1I("h1_piOcc_BH2",
    "HTOF 1D: #pi tile occupancy | BH2 4-10;tile;events", NT, -0.5, NT-0.5);
  TH1I* h1_piOcc_BH2_MPge2    = new TH1I("h1_piOcc_BH2_MPge2",
    "HTOF 1D: #pi tile occupancy | BH2 4-10 & MP#geq2;tile;events", NT, -0.5, NT-0.5);

  // ---------------- Event loop ----------------
  const Long64_t N = T->GetEntries();
  Long64_t N_total=0, N_bh2In=0, N_bh2In_mpge2=0;

  for(Long64_t ie=0; ie<N; ++ie){
    T->GetEntry(ie); N_total++;

    // ---- BH2 in-range (denominator) ----
    std::map<int,double> bh2E;
    if(BH2){
      const size_t n = BH2->size();
      for(size_t i=0;i<n;++i){
        const TParticle& p = BH2->at(i);
        const int sid = MapBH2_WorldToSeg(p.Vx(), p.Vy(), p.Vz());
        if(0<=sid && sid<kNBH2Seg){
          const double ed = (hasBH2edep && BH2_edep && i<BH2_edep->size())
                            ? BH2_edep->at(i) : p.GetWeight();
          bh2E[sid] += ed;
        }
      }
    }
    bool bh2_in_range=false;
    for(const auto& kv : bh2E){
      if(kv.second >= thrBH2){
        const int s = kv.first;
        if(bh2_lo<=s && s<=bh2_hi){ bh2_in_range=true; break; }
      }
    }
    if(!bh2_in_range) continue; // event rejected
    N_bh2In++;

    // ---- HTOF valid tiles per species (ΣE per tile >= thr) ----
    // PDG codes
    constexpr int PDG_PIP  = +211;
    constexpr int PDG_PIM  = -211;
    constexpr int PDG_P    = +2212;
    constexpr int PDG_PBAR = -2212;

    // 타일별 에너지 합(종별로 분리)
    std::map<int,double> E_tile_pip, E_tile_pim, E_tile_p, E_tile_pbar;
    std::map<int,double> E_tile_all; // multiplicity용 (종합)

    if(HTOF){
      const size_t n = HTOF->size();
      for(size_t i=0;i<n;++i){
        const TParticle& p = HTOF->at(i);
        const int tile = p.GetStatusCode(); // copy-no
        if(tile<0 || tile>=NT) continue;

        const double ed = (hasHTOFedep && HTOF_edep && i<HTOF_edep->size())
                          ? HTOF_edep->at(i) : p.GetWeight();

        // 전체 multiplicity 집계
        E_tile_all[tile] += ed;

        // 종별 집계
        const int pdg = p.GetPdgCode(); // HTSD에서 저장해 준 경우
        if      (pdg==PDG_PIP ) E_tile_pip [tile] += ed;
        else if (pdg==PDG_PIM ) E_tile_pim [tile] += ed;
        else if (pdg==PDG_P   ) E_tile_p   [tile] += ed;
        else if (pdg==PDG_PBAR) E_tile_pbar[tile] += ed;
      }
    }

    // VALID 타일 집합
    std::set<int> tiles_all, tiles_pip, tiles_pim, tiles_p, tiles_pbar;
    for(const auto& kv: E_tile_all) if(kv.second >= thrHTOF) tiles_all.insert(kv.first);
    for(const auto& kv: E_tile_pip) if(kv.second >= thrHTOF) tiles_pip.insert(kv.first);
    for(const auto& kv: E_tile_pim) if(kv.second >= thrHTOF) tiles_pim.insert(kv.first);
    for(const auto& kv: E_tile_p  ) if(kv.second >= thrHTOF) tiles_p  .insert(kv.first);
    for(const auto& kv: E_tile_pbar)if(kv.second >= thrHTOF) tiles_pbar.insert(kv.first);

    const int mp_all = (int)tiles_all.size();

    // ===== (1) BH2만 적용(무조건) =====
    // π+ vs π−
    for(int a: tiles_pip) for(int b: tiles_pim) h2_pip_pim_BH2->Fill(a,b);
    // p vs pbar
    for(int a: tiles_p) for(int b: tiles_pbar) h2_p_pbar_BH2->Fill(a,b);
    // π+ vs p
    for(int a: tiles_pip) for(int b: tiles_p) h2_pip_p_BH2->Fill(a,b);
    // π(± 합계) 1D
    {
      std::set<int> tiles_pi;
      tiles_pi.insert(tiles_pip.begin(), tiles_pip.end());
      tiles_pi.insert(tiles_pim.begin(), tiles_pim.end());
      for(int t: tiles_pi) h1_piOcc_BH2->Fill(t);
    }

    // ===== (2) BH2 && MP >= 2 =====
    if(mp_all >= 2){
      N_bh2In_mpge2++;
      for(int a: tiles_pip) for(int b: tiles_pim) h2_pip_pim_BH2_MPge2->Fill(a,b);
      for(int a: tiles_p)   for(int b: tiles_pbar) h2_p_pbar_BH2_MPge2->Fill(a,b);
      for(int a: tiles_pip) for(int b: tiles_p)    h2_pip_p_BH2_MPge2->Fill(a,b);
      {
        std::set<int> tiles_pi;
        tiles_pi.insert(tiles_pip.begin(), tiles_pip.end());
        tiles_pi.insert(tiles_pim.begin(), tiles_pim.end());
        for(int t: tiles_pi) h1_piOcc_BH2_MPge2->Fill(t);
      }
    }
  }

  // ---------------- Summary ----------------
  auto pct = [](long long a, long long b)->double{ return (b>0)? (100.0*(double)a/(double)b):0.0; };

  std::cout<<"\n==== Summary (BH2 "<<bh2_lo<<"-"<<bh2_hi<<") ====\n";
  std::cout<<"Total events                 : "<<N_total<<"\n";
  std::cout<<"BH2 in range (selected)      : "<<N_bh2In<<"\n";
  std::cout<<"BH2 in range & HTOF MP>=2    : "<<N_bh2In_mpge2
           <<"  ("<<std::setprecision(3)<<pct(N_bh2In_mpge2,N_bh2In)<<" % of BH2-selected)\n";

  // ---------------- Draw ----------------
  auto drawAndSave2D = [&](TH2I* h, const char* cname, const char* fname){
    TCanvas* c = new TCanvas(cname, cname, 800, 700);
    h->SetStats(0);
    h->GetXaxis()->SetNdivisions(34,false);
    h->GetYaxis()->SetNdivisions(34,false);
    h->Draw("COLZ");
    if(save) c->SaveAs(fname);
  };
  auto drawAndSave1D = [&](TH1I* h, const char* cname, const char* fname){
    TCanvas* c = new TCanvas(cname, cname, 900, 350);
    h->SetStats(0); h->SetLineWidth(2);
    h->Draw("hist");
    if(save) c->SaveAs(fname);
  };

  TString t(tag); if(t.IsNull()) t = Form("BH2_%d_%d", bh2_lo, bh2_hi);

  drawAndSave2D(h2_pip_pim_BH2,       "c_pip_pim_BH2",       "HTOF_2D_pip_vs_pim_BH2_"+t+".png");
  drawAndSave2D(h2_pip_pim_BH2_MPge2, "c_pip_pim_BH2_MPge2", "HTOF_2D_pip_vs_pim_BH2_MPge2_"+t+".png");

  drawAndSave2D(h2_p_pbar_BH2,        "c_p_pbar_BH2",        "HTOF_2D_p_vs_pbar_BH2_"+t+".png");
  drawAndSave2D(h2_p_pbar_BH2_MPge2,  "c_p_pbar_BH2_MPge2",  "HTOF_2D_p_vs_pbar_BH2_MPge2_"+t+".png");

  drawAndSave2D(h2_pip_p_BH2,         "c_pip_p_BH2",         "HTOF_2D_pip_vs_p_BH2_"+t+".png");
  drawAndSave2D(h2_pip_p_BH2_MPge2,   "c_pip_p_BH2_MPge2",   "HTOF_2D_pip_vs_p_BH2_MPge2_"+t+".png");

  drawAndSave1D(h1_piOcc_BH2,         "c_piOcc_BH2",         "HTOF_1D_piOcc_BH2_"+t+".png");
  drawAndSave1D(h1_piOcc_BH2_MPge2,   "c_piOcc_BH2_MPge2",   "HTOF_1D_piOcc_BH2_MPge2_"+t+".png");
}

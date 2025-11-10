// -*- C++ -*-
// HTOF_pair_and_tile_scan_14_27.C  (2025-11-09, histo added)
// Section 1: 인접 페어 (14,15) .. (26,27) %Beam (no-excl vs excl{2..5}) → Canvas A
// Section 2: 단일 타일 14..27  % of Σtiles(MP>=1) (no-excl vs excl{2..5}) → Canvas B
//
// 사용법:
//   root -l
//   .L HTOF_pair_and_tile_scan_14_27.C+
//   HTOF_pair_and_tile_scan_14_27("E45.root","g4hyptpc", 4,10, 0.10,2.0, 5.0,10.0);

#include "TFile.h"
#include "TTree.h"
#include "TParticle.h"
#include "TSystem.h"
#include "TInterpreter.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLine.h"
#include <vector>
#include <map>
#include <set>
#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>

namespace HPS {

static const int    kNBH2Seg   = 15;     // 0..14
static const double kBH2_x0    = -10.0;  // [mm]
static const double kBH2_pitch = 14.0;   // [mm]
static const double kBH2_total = kNBH2Seg * kBH2_pitch;
static const int    kNHTOF     = 34;     // 0..33

struct DictGuard {
  DictGuard(){
    gSystem->Load("libPhysics");
    gInterpreter->GenerateDictionary("vector<TParticle>","TParticle.h;vector");
  }
};
static DictGuard _dg;

static inline double MIPThrMeV(double frac,double dEdx,double thk_mm){
  return frac * dEdx * (thk_mm*0.1); // mm→cm
}
static inline int MapBH2_WorldToSeg(double x, double, double){
  const double xloc = x - kBH2_x0;
  const double xmin = -0.5*kBH2_total, xmax = 0.5*kBH2_total;
  if (xloc < xmin || xloc >= xmax) return -1;
  int idx = int(std::floor((xloc - xmin)/kBH2_pitch));
  if(idx<0) idx=0; if(idx>=kNBH2Seg) idx=kNBH2Seg-1;
  return idx;
}
static inline bool PairHit(const std::set<int>& tiles, int a, int b){
  return tiles.count(a) && tiles.count(b);
}

} // namespace HPS

void HTOF_pair_and_tile_scan_14_27(const char* filename,
                                   const char* treename="g4hyptpc",
                                   int bh2_lo=4, int bh2_hi=10,
                                   double mipFrac=0.10,
                                   double mipMeVperCm=2.0,
                                   double BH2_thickness_mm=5.0,
                                   double HTOF_thickness_mm=10.0)
{
  using namespace HPS;

  // ROOT style (가독성)
  gStyle->SetOptStat(0);
  gStyle->SetTitleFontSize(0.045);
  gStyle->SetLabelSize(0.035,"XY");
  gStyle->SetTitleSize(0.045,"XY");

  // open & tree
  TFile* f=TFile::Open(filename,"READ");
  if(!f||f->IsZombie()){ std::cerr<<"[ERR] open "<<filename<<" failed\n"; return; }
  TTree* T=(TTree*)f->Get(treename);
  if(!T){ std::cerr<<"[ERR] tree "<<treename<<" not found\n"; return; }
  if(!T->GetBranch("BH2") || !T->GetBranch("HTOF")){
    std::cerr<<"[ERR] need BH2 & HTOF (vector<TParticle>) branches\n"; return;
  }
  const bool hasBH2edep  = (T->GetBranch("BH2_edep")  != nullptr);
  const bool hasHTOFedep = (T->GetBranch("HTOF_edep") != nullptr);

  std::vector<TParticle>* BH2=nullptr;   T->SetBranchAddress("BH2",&BH2);
  std::vector<TParticle>* HTOF=nullptr;  T->SetBranchAddress("HTOF",&HTOF);
  std::vector<double>* BH2_edep=nullptr; if(hasBH2edep)  T->SetBranchAddress("BH2_edep",&BH2_edep);
  std::vector<double>* HTOF_edep=nullptr;if(hasHTOFedep) T->SetBranchAddress("HTOF_edep",&HTOF_edep);

  // thresholds
  const double thrBH2  = MIPThrMeV(mipFrac,mipMeVperCm,BH2_thickness_mm);
  const double thrHTOF = MIPThrMeV(mipFrac,mipMeVperCm,HTOF_thickness_mm);

  // 인접 페어 (14,15) .. (26,27) 13개
  std::vector<std::pair<int,int>> pairs;
  for(int a=14; a<=26; ++a) pairs.emplace_back(a, a+1);

  // 단일 타일 14..27 (14개)
  std::vector<int> tiles_list;
  for(int t=14; t<=27; ++t) tiles_list.push_back(t);

  // 카운터
  long long N_total=0, N_beam=0;

  // Section1: pairs
  std::vector<long long> pair_noexcl(pairs.size(), 0);
  std::vector<long long> pair_excl  (pairs.size(), 0);

  // Section2: tiles + denominators
  std::vector<long long> tile_noexcl(tiles_list.size(), 0);
  std::vector<long long> tile_excl  (tiles_list.size(), 0);
  long long denom_noexcl = 0;
  long long denom_excl   = 0;

  const std::set<int> EXCL_2_5 = {2,3,4,5};

  const Long64_t N=T->GetEntries();
  for(Long64_t ie=0; ie<N; ++ie){
    T->GetEntry(ie); N_total++;

    // ---- BH2: Beam in-range ----
    std::map<int,double> bh2E;
    if(BH2){
      for(size_t i=0;i<BH2->size();++i){
        const TParticle& p=BH2->at(i);
        int sid = MapBH2_WorldToSeg(p.Vx(),p.Vy(),p.Vz());
        if(0<=sid && sid<kNBH2Seg){
          const double ed=(hasBH2edep && BH2_edep && i<BH2_edep->size())? BH2_edep->at(i) : p.GetWeight();
          bh2E[sid]+=ed;
        }
      }
    }
    std::set<int> bh2Valid;
    for(const auto& kv:bh2E) if(kv.second>=thrBH2) bh2Valid.insert(kv.first);
    bool beam=false; for(int s:bh2Valid){ if(bh2_lo<=s && s<=bh2_hi){ beam=true; break; } }
    if(!beam) continue;
    N_beam++;

    // ---- HTOF 누적 (no-excl / excl{2..5}) ----
    std::map<int,double> acc_noexcl, acc_excl;
    if(HTOF){
      for(size_t i=0;i<HTOF->size();++i){
        const TParticle& p=HTOF->at(i);
        const int tid = p.GetStatusCode();
        if(!(0<=tid && tid<kNHTOF)) continue;
        const double ed=(hasHTOFedep && HTOF_edep && i<HTOF_edep->size())? HTOF_edep->at(i) : p.GetWeight();
        acc_noexcl[tid] += ed;
        if(!EXCL_2_5.count(tid)) acc_excl[tid] += ed;
      }
    }
    std::set<int> tiles_noexcl, tiles_excl;
    for(const auto& kv:acc_noexcl) if(kv.second>=thrHTOF) tiles_noexcl.insert(kv.first);
    for(const auto& kv:acc_excl)   if(kv.second>=thrHTOF) tiles_excl.insert(kv.first);

    // ---- Section 1: 페어 카운트 (Beam 기준) ----
    for(size_t ip=0; ip<pairs.size(); ++ip){
      const int a = pairs[ip].first;
      const int b = pairs[ip].second;
      if(PairHit(tiles_noexcl,a,b)) pair_noexcl[ip]++;
      if(PairHit(tiles_excl,  a,b)) pair_excl[ip]++;
    }

    // ---- Section 2: 단일 타일 점유도 (MP>=1에서 집계) ----
    if(!tiles_noexcl.empty()) denom_noexcl += (long long)tiles_noexcl.size();
    if(!tiles_excl.empty())   denom_excl   += (long long)tiles_excl.size();

    for(size_t it=0; it<tiles_list.size(); ++it){
      int t = tiles_list[it];
      if(tiles_noexcl.count(t)) tile_noexcl[it]++;
      if(tiles_excl.count(t))   tile_excl[it]++;
    }
  }

  // ---- 터미널 출력 (요약) ----
  auto pct_beam = [&](long long x)->double{
    return (N_beam>0)? (100.0*double(x)/double(N_beam)) : 0.0;
  };
  auto pct_denom = [&](long long x, long long denom)->double{
    return (denom>0)? (100.0*double(x)/double(denom)) : 0.0;
  };

  std::cout<<std::fixed<<std::setprecision(3);
  std::cout<<"[INFO] BH2 thr="<<thrBH2<<" MeV, HTOF thr="<<thrHTOF<<" MeV\n";
  std::cout<<"[INFO] BH2 window: "<<bh2_lo<<"–"<<bh2_hi<<"\n";
  std::cout<<"Total events : "<<N_total<<"\n";
  std::cout<<"Beam (BH2 in-range) : "<<N_beam<<"\n\n";

  std::cout<<"== Section 1: Adjacent-pair hits among Beam ==\n";
  std::cout<<" Pair    no-excl(count, %Beam)      excl{2–5}(count, %Beam)\n";
  for(size_t ip=0; ip<pairs.size(); ++ip){
    const auto& pr = pairs[ip];
    std::cout<<" ("<<std::setw(2)<<pr.first<<","<<std::setw(2)<<pr.second<<") : "
             <<std::setw(8)<<pair_noexcl[ip]<<" ("<<std::setw(7)<<pct_beam(pair_noexcl[ip])<<"%)"
             <<"    "
             <<std::setw(8)<<pair_excl[ip]  <<" ("<<std::setw(7)<<pct_beam(pair_excl[ip])  <<"%)\n";
  }
  std::cout<<"\n";

  std::cout<<"== Section 2: Single-tile occupancy among Beam (MP>=1 only)\n";
  std::cout<<" Denominator (no-excl) = Σ valid tiles over events with MP>=1 = "<<denom_noexcl<<"\n";
  std::cout<<" Denominator (excl{2–5}) = Σ valid tiles over events with MP>=1 = "<<denom_excl<<"\n";
  std::cout<<" Tile   no-excl(count, %of Σtiles_MP>=1)    excl{2–5}(count, %of Σtiles_MP>=1)\n";
  for(size_t it=0; it<tiles_list.size(); ++it){
    int t = tiles_list[it];
    std::cout<<"  "<<std::setw(2)<<t<<"   "
             <<std::setw(8)<<tile_noexcl[it]<<" ("<<std::setw(7)<<pct_denom(tile_noexcl[it], denom_noexcl)<<"%)"
             <<"            "
             <<std::setw(8)<<tile_excl[it]  <<" ("<<std::setw(7)<<pct_denom(tile_excl[it],   denom_excl)  <<"%)\n";
  }

  // ====== 히스토그램 생성 ======
  // Section 1: Pair %Beam
  const int nPairs = (int)pairs.size();
  TH1D* hPair_noex = new TH1D("hPair_noex","Section 1: Adjacent pairs %Beam;HTOF pair;Fraction of Beam [%]", nPairs, 0.5, nPairs+0.5);
  TH1D* hPair_excl = new TH1D("hPair_excl","Section 1: Adjacent pairs %Beam;HTOF pair;Fraction of Beam [%]", nPairs, 0.5, nPairs+0.5);
  for(int i=1;i<=nPairs;++i){
    const int a = pairs[i-1].first;
    const int b = pairs[i-1].second;
    hPair_noex->GetXaxis()->SetBinLabel(i, Form("(%d,%d)",a,b));
    hPair_excl->GetXaxis()->SetBinLabel(i, Form("(%d,%d)",a,b));
    hPair_noex->SetBinContent(i, pct_beam(pair_noexcl[i-1]));
    hPair_excl->SetBinContent(i, pct_beam(pair_excl[i-1]));
  }
  hPair_noex->SetLineColor(kBlue+1);
  hPair_noex->SetMarkerColor(kBlue+1);
  hPair_noex->SetMarkerStyle(20);
  hPair_excl->SetLineColor(kRed+1);
  hPair_excl->SetMarkerColor(kRed+1);
  hPair_excl->SetMarkerStyle(24);

  TCanvas* cA = new TCanvas("cA","Pairs %Beam", 1100, 420);
  cA->SetLeftMargin(0.08); cA->SetBottomMargin(0.25); // 라벨 겹침 방지
  hPair_noex->SetMinimum(0);
  hPair_noex->Draw("HIST");
  hPair_excl->Draw("HIST SAME");

  auto legA = new TLegend(0.78,0.78,0.98,0.98);
  legA->AddEntry(hPair_noex,"no-excl","l");
  legA->AddEntry(hPair_excl,"excl{2-5}","l");
  legA->Draw();

  cA->SaveAs("HTOF_pairs_14_27_percentBeam.png");

  // Section 2: Tiles % of Σtiles(MP>=1)
  const int nTiles = (int)tiles_list.size();
  TH1D* hTile_noex = new TH1D("hTile_noex","Section 2: Tile share % of #Sigma tiles (MP#geq1);HTOF tile;Share [%]", nTiles, 0.5, nTiles+0.5);
  TH1D* hTile_excl = new TH1D("hTile_excl","Section 2: Tile share % of #Sigma tiles (MP#geq1);HTOF tile;Share [%]", nTiles, 0.5, nTiles+0.5);
  for(int i=1;i<=nTiles;++i){
    const int t = tiles_list[i-1];
    hTile_noex->GetXaxis()->SetBinLabel(i, Form("%d",t));
    hTile_excl->GetXaxis()->SetBinLabel(i, Form("%d",t));
    hTile_noex->SetBinContent(i, pct_denom(tile_noexcl[i-1], denom_noexcl));
    hTile_excl->SetBinContent(i,   pct_denom(tile_excl[i-1],   denom_excl));
  }
  hTile_noex->SetLineColor(kBlue+1);
  hTile_noex->SetMarkerColor(kBlue+1);
  hTile_noex->SetMarkerStyle(20);
  hTile_excl->SetLineColor(kRed+1);
  hTile_excl->SetMarkerColor(kRed+1);
  hTile_excl->SetMarkerStyle(24);

  TCanvas* cB = new TCanvas("cB","Tiles % of sum", 1100, 420);
  cB->SetLeftMargin(0.08); cB->SetBottomMargin(0.18);
  hTile_noex->SetMinimum(0);
  hTile_noex->Draw("HIST");
  hTile_excl->Draw("HIST SAME");

  auto legB = new TLegend(0.78,0.78,0.98,0.98);
  legB->AddEntry(hTile_noex,"no-excl","l");
  legB->AddEntry(hTile_excl,"excl{2-5}","l");
  legB->Draw();

  cB->SaveAs("HTOF_tiles_14_27_shareOfSum.png");

  std::cout<<"\n[Saved] HTOF_pairs_14_27_percentBeam.png\n";
  std::cout<<"[Saved] HTOF_tiles_14_27_shareOfSum.png\n";
}

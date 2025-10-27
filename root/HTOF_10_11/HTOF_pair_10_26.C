// HTOF_pair_10_26.C
// Count HTOF adjacent-pair hits within events that pass:
//   (1) BH2 in [bh2_lo, bh2_hi]  AND  HTOF multiplicity >= 2   --> MP>=2 set
//   (2) BH2 in [bh2_lo, bh2_hi]  AND  HTOF multiplicity >  2   --> MP>2 set
//
// Pairs counted (inclusive per event; unique per event for each pair):
//   (13,14), (14,15), (15,16), (16,17), (17,18), (18,19), (19,20),
//   (20,21), (21,22), (22,23), (23,24), (24,25), (25,26)
//
// Notes / Policies:
//   - HTOF tile ID = TParticle::StatusCode() (copy-no written in SD)
//   - VALID hit: sum(edep) >= threshold (use *_edep if exists else Weight())
//   - BH2 seg id is mapped from world x using E72-like geometry numbers here.

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1I.h"
#include "TLegend.h"
#include "TSystem.h"
#include "TInterpreter.h"
#include "TParticle.h"
#include "TString.h"
#include <vector>
#include <map>
#include <set>
#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>

namespace P10_26 {

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

} // namespace P10_26


void HTOF_pair_10_26(const char* filename="E45.root",
                      const char* treename="g4hyptpc",
                      int    bh2_lo=4, int bh2_hi=10,
                      double mipFrac=0.10,
                      double mipMeVperCm=2.0,
                      double BH2_thickness_mm=5.0,
                      double HTOF_thickness_mm=10.0,
                      bool   save=true,
                      const char* tag="BH2_4_10")
{
  using namespace P10_26;

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

  // ---------------- Thresholds & pairs ----------------
  const double thrBH2  = MIPThrMeV(mipFrac,mipMeVperCm,BH2_thickness_mm);   // ~0.10 MeV
  const double thrHTOF = MIPThrMeV(mipFrac,mipMeVperCm,HTOF_thickness_mm);  // ~0.20 MeV

  std::cout<<std::fixed<<std::setprecision(3)
           <<"[INFO] BH2 thr="<<thrBH2<<" MeV, HTOF thr="<<thrHTOF<<" MeV | "
           <<"BH2 range = "<<bh2_lo<<"-"<<bh2_hi<<"\n";
  std::cout<<"[INFO] HTOF tile = StatusCode(copy-no). VALID hit = sum(edep) >= thr.\n";

  // Pairs to check (13..26 adjacent)
  const std::pair<int,int> pairs[] = {
    {13,14},{14,15},{15,16},{16,17},{17,18},{18,19},{19,20},
    {20,21},{21,22},{22,23},{23,24},{24,25},{25,26}
  };
  const int NP = sizeof(pairs)/sizeof(pairs[0]);

  // ---------------- Accumulators ----------------
  Long64_t N_total=0, N_bh2In=0;
  Long64_t N_mp_ge2=0, N_mp_gt2=0;

  // Pair counters (unique-per-event per set)
  std::vector<Long64_t> pairCnt_ge2(NP,0), pairCnt_gt2(NP,0);

  // For optional histograms
  TH1I* hPairs_ge2 = new TH1I("hPairs_ge2",
      Form("HTOF pairs in events with MP#geq2 (BH2 %d-%d);Pair index;Events",bh2_lo,bh2_hi),
      NP, 0.5, NP+0.5);
  TH1I* hPairs_gt2 = new TH1I("hPairs_gt2",
      Form("HTOF pairs in events with MP>2 (BH2 %d-%d);Pair index;Events",bh2_lo,bh2_hi),
      NP, 0.5, NP+0.5);
  hPairs_ge2->SetDirectory(nullptr);
  hPairs_gt2->SetDirectory(nullptr);

  // X-axis labels like "(13,14)" ...
  for(int i=1;i<=NP;++i){
    hPairs_ge2->GetXaxis()->SetBinLabel(i, Form("(%d,%d)", pairs[i-1].first, pairs[i-1].second));
    hPairs_gt2->GetXaxis()->SetBinLabel(i, Form("(%d,%d)", pairs[i-1].first, pairs[i-1].second));
  }

  // ---------------- Event loop ----------------
  const Long64_t N = T->GetEntries();
  for(Long64_t ie=0; ie<N; ++ie){
    T->GetEntry(ie); N_total++;

    // ---- BH2 valid segments (for event pre-selection) ----
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

    // ---- HTOF valid tiles by copy-no (StatusCode) ----
    std::map<int,double> tileE;
    if(HTOF){
      const size_t n = HTOF->size();
      for(size_t i=0;i<n;++i){
        const TParticle& p = HTOF->at(i);
        const int tile = p.GetStatusCode();      // copy-no based tile ID
        if(tile < 0 || tile > 33) continue;
        const double ed = (hasHTOFedep && HTOF_edep && i<HTOF_edep->size())
                          ? HTOF_edep->at(i) : p.GetWeight();
        tileE[tile] += ed;
      }
    }

    // Build valid-tile set & multiplicity
    std::set<int> tilesValid;
    for(const auto& kv : tileE){
      if(kv.second >= thrHTOF) tilesValid.insert(kv.first);
    }
    const int mult = (int)tilesValid.size();

    // ---- MP ≥ 2 set ----
    if(mult >= 2){
      N_mp_ge2++;
      // Count each pair once per event if both tiles present
      for(int ip=0; ip<NP; ++ip){
        const auto [a,b] = pairs[ip];
        if(tilesValid.count(a) && tilesValid.count(b)){
          pairCnt_ge2[ip]++; hPairs_ge2->Fill(ip+1);
        }
      }
    }

    // ---- MP > 2 set ----
    if(mult > 2){
      N_mp_gt2++;
      for(int ip=0; ip<NP; ++ip){
        const auto [a,b] = pairs[ip];
        if(tilesValid.count(a) && tilesValid.count(b)){
          pairCnt_gt2[ip]++; hPairs_gt2->Fill(ip+1);
        }
      }
    }
  }

  // ---------------- Print summary ----------------
  auto pct = [](long long a, long long b)->double{ return (b>0)? (100.0*(double)a/(double)b):0.0; };

  std::cout<<"\n==== Summary (BH2 "<<bh2_lo<<"-"<<bh2_hi<<") ====\n";
  std::cout<<"Total events            : "<<N_total<<"\n";
  std::cout<<"BH2 in range (selected) : "<<N_bh2In<<"\n";

  std::cout<<std::setprecision(3);
  std::cout<<"HTOF multiplicity >= 2  : "<<N_mp_ge2<<"  ("<<pct(N_mp_ge2,N_bh2In)<<" % of BH2-selected)\n";
  for(int ip=0; ip<NP; ++ip){
    const auto& pr = pairs[ip];
    const long long n = pairCnt_ge2[ip];
    const double p_mp2 = pct(n, N_mp_ge2);
    const double p_bh2 = pct(n, N_bh2In);
    std::cout<<"  Pair "<<std::setw(2)<<pr.first<<","<<std::setw(2)<<pr.second
             <<" : "<<n<<"  ("<<p_mp2<<" % of MP>=2, "<<p_bh2<<" % of BH2-selected)\n";
  }

  std::cout<<"HTOF multiplicity >  2  : "<<N_mp_gt2<<"  ("<<pct(N_mp_gt2,N_bh2In)<<" % of BH2-selected)\n";
  for(int ip=0; ip<NP; ++ip){
    const auto& pr = pairs[ip];
    const long long n = pairCnt_gt2[ip];
    const double p_mp3 = pct(n, N_mp_gt2); // MP>2 모임 대비
    const double p_bh2 = pct(n, N_bh2In);  // BH2-selected 대비
    std::cout<<"  Pair "<<std::setw(2)<<pr.first<<","<<std::setw(2)<<pr.second
             <<" : "<<n<<"  ("<<p_mp3<<" % of MP>2, "<<p_bh2<<" % of BH2-selected)\n";
  }

  // ---------------- Draw (optional) ----------------
  TCanvas* c1 = new TCanvas("cPairs_ge2","Pairs in MP#geq2", 1000, 450);
  hPairs_ge2->SetLineWidth(2);
  hPairs_ge2->LabelsOption("v","X");
  hPairs_ge2->Draw("hist");

  TCanvas* c2 = new TCanvas("cPairs_gt2","Pairs in MP>2", 1000, 450);
  hPairs_gt2->SetLineWidth(2);
  hPairs_gt2->LabelsOption("v","X");
  hPairs_gt2->Draw("hist");

  if(save){
    TString t(tag);
    if(t.IsNull()) t = Form("BH2_%d_%d", bh2_lo, bh2_hi);
    c1->SaveAs("HTOF_pairs_MPge2_"+t+".png");
    c2->SaveAs("HTOF_pairs_MPgt2_"+t+".png");
  }
}

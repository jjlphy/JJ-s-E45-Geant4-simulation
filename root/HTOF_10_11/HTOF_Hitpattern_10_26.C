// HTOF_HitPat34_copyNo.C — Draw 34-tile HTOF hit pattern using copy-no in StatusCode()
// Only events with BH2 valid-hit in [bh2_lo, bh2_hi] are counted.
//
// Usage:
//   root -l
//   .L HTOF_HitPat34_copyNo.C+
//   HTOF_HitPat34_copyNo("E45.root","g4hyptpc", /*bh2_lo=*/4, /*bh2_hi=*/10,
//                        /*mipFrac=*/0.10, /*mipMeVperCm=*/2.0,
//                        /*BH2_thk_mm=*/5.0, /*HTOF_thk_mm=*/10.0,
//                        /*UNIQUE_PER_EVENT=*/true, /*save=*/true, "BH2_4_10");

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1D.h"
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

namespace HP34 {

// ----- dictionary for vector<TParticle> -----
struct DictGuard {
  DictGuard(){
    gSystem->Load("libPhysics");
    gInterpreter->GenerateDictionary("vector<TParticle>","TParticle.h;vector");
  }
};
static DictGuard _dg;

// ----- geometry (BH2 mapping; E72-like defaults as in your code) -----
static const int    kNBH2Seg   = 15;     // 0..14
static const double kBH2_x0    = -10.0;  // [mm]
static const double kBH2_pitch = 14.0;   // [mm]
static const double kBH2_total = kNBH2Seg * kBH2_pitch;

// world (x,*,*) → BH2 seg id
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

} // namespace HP34


void HTOF_HitPat34_copyNo(const char* filename="E45.root",
                          const char* treename="g4hyptpc",
                          int    bh2_lo=4, int bh2_hi=10,
                          double mipFrac=0.10,
                          double mipMeVperCm=2.0,
                          double BH2_thickness_mm=5.0,
                          double HTOF_thickness_mm=10.0,
                          bool   UNIQUE_PER_EVENT=true,
                          bool   save=true,
                          const char* tag="BH2_4_10")
{
  using namespace HP34;

  // Open file & tree
  TFile* f=TFile::Open(filename,"READ");
  if(!f||f->IsZombie()){ std::cerr<<"[ERR] open "<<filename<<" failed\n"; return; }
  TTree* T=(TTree*)f->Get(treename);
  if(!T){ std::cerr<<"[ERR] tree "<<treename<<" not found\n"; return; }

  // Required branches
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

  // Thresholds (policy)
  const double thrBH2  = MIPThrMeV(mipFrac,mipMeVperCm,BH2_thickness_mm);   // ~0.10 MeV
  const double thrHTOF = MIPThrMeV(mipFrac,mipMeVperCm,HTOF_thickness_mm);  // ~0.20 MeV

  std::cout<<std::fixed<<std::setprecision(3)
           <<"[INFO] BH2 thr="<<thrBH2<<" MeV, HTOF thr="<<thrHTOF<<" MeV | "
           <<"BH2 range = "<<bh2_lo<<"-"<<bh2_hi<<"\n";
  std::cout<<"[INFO] ID policy: HTOF tile = TParticle::StatusCode() (copy-no). "
              "VALID hit = sum(edep) >= thr (use *_edep if exists else Weight).\n";

  // Histogram (34 tiles: 0..33)
  TH1D* h = new TH1D("hHTOF_tile34_copyNo",
                     Form("HTOF 34-tile Hit Pattern (copy-no) | BH2 %d-%d;Tile ID (0..33);Counts",
                          bh2_lo,bh2_hi),
                     34, -0.5, 33.5);
  h->SetDirectory(nullptr);

  Long64_t N_tot=0, N_bh2In=0, N_filled=0;

  const Long64_t N = T->GetEntries();
  for(Long64_t ie=0; ie<N; ++ie){
    T->GetEntry(ie);
    N_tot++;

    // ---------- BH2 valid segs (for event filter) ----------
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
    if(!bh2_in_range) continue; // event vetoed by BH2 condition
    N_bh2In++;

    // ---------- HTOF: collect tile energy by copy-no (StatusCode) ----------
    std::map<int,double> tileE;
    if(HTOF){
      const size_t n = HTOF->size();
      for(size_t i=0;i<n;++i){
        const TParticle& p = HTOF->at(i);
        const int tile = p.GetStatusCode();         // ★ copy-no based tile ID
        if(tile<0 || tile>33) continue;             // guard
        const double ed = (hasHTOFedep && HTOF_edep && i<HTOF_edep->size())
                          ? HTOF_edep->at(i) : p.GetWeight();
        tileE[tile] += ed;
      }
    }

    // ---------- fill histogram ----------
    std::set<int> filled_this_event;
    for(const auto& kv : tileE){
      const int    tile = kv.first;
      const double sumE = kv.second;
      if(sumE < thrHTOF) continue;                  // VALID hit only
      if(UNIQUE_PER_EVENT){
        if(filled_this_event.insert(tile).second){
          h->Fill(tile); N_filled++;
        }
      }else{
        h->Fill(tile); N_filled++;
      }
    }
  }

  // Draw
  TCanvas* c = new TCanvas("cHTOF34_copyNo","HTOF 34-tile Hit Pattern (copy-no)", 900, 500);
  h->SetLineWidth(2);
  h->Draw("hist");
  c->Update();

  if(save){
    TString t(tag);
    if(t.IsNull()) t = Form("BH2_%d_%d", bh2_lo, bh2_hi);
    c->SaveAs("HTOF_tile34_copyNo_"+t+".png");
  }

  std::cout<<"[SUM] Total events     : "<<N_tot<<"\n";
  std::cout<<"[SUM] BH2-in-range evt : "<<N_bh2In<<"\n";
  std::cout<<"[SUM] Filled entries   : "<<N_filled
           <<"  ("<<(h->GetEntries())<<" histogram entries)\n";
}

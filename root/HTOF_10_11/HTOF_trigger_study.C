// HTOF_trigger_study.C
// Logic:
//   Trig1 = (BH2 in [lo,hi]) AND (HTOF multiplicity >= 2)
//   Trig2 = Trig1 AND NOT BeamVeto
//   BeamVeto = (tile20 AND tile21) OR (tile22 AND tile23)   // 포함관계: ≥2에서도 적용
//
// Denominator for % : #events with BH2 in [lo,hi]
// Policies: HTOF ID = TParticle::StatusCode(), BH2 ID = coordinate→segment.
//           edep: *_edep branch if exists else TParticle::Weight().
//           VALID hit requires sum(edep) >= threshold (no fallback).
// Cuts: default 0.1 MIP with dE/dx=2 MeV/cm, BH2=5 mm, HTOF=10 mm.

#include "TFile.h"
#include "TTree.h"
#include "TParticle.h"
#include "TSystem.h"
#include "TInterpreter.h"
#include "TH1I.h"
#include "TCanvas.h"
#include "TString.h"
#include <vector>
#include <map>
#include <set>
#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>

namespace HTS {

// geometry (E72-like)
static const int    kNBH2Seg   = 15;     // 0..14
static const double kBH2_x0    = -10.0;  // [mm]
static const double kBH2_pitch = 14.0;   // [mm]
static const double kBH2_total = kNBH2Seg * kBH2_pitch;
static const int    kNHTOF     = 34;     // 0..33

struct DictGuard {
  DictGuard(){ gSystem->Load("libPhysics");
               gInterpreter->GenerateDictionary("vector<TParticle>","TParticle.h;vector");}
};
static DictGuard _dg;

static inline double MIPThrMeV(double frac,double dEdx,double thk_mm){
  return frac * dEdx * (thk_mm*0.1); // mm→cm
}

// x(world) → BH2 seg id
static int MapBH2_WorldToSeg(double x, double, double){
  const double xloc = x - kBH2_x0;
  const double xmin = -0.5*kBH2_total, xmax = 0.5*kBH2_total;
  if (xloc < xmin || xloc >= xmax) return -1;
  int idx = int(std::floor((xloc - xmin)/kBH2_pitch));
  if(idx<0) idx=0; if(idx>=kNBH2Seg) idx=kNBH2Seg-1;
  return idx;
}

static inline std::string JoinTuple(const std::set<int>& S){
  std::string out="(";
  bool first=true; for(int v:S){ if(!first) out+=','; first=false; out+=std::to_string(v); }
  out+=')'; return out;
}

} // namespace HTS


void HTOF_trigger_study(const char* filename,
                        const char* treename="g4hyptpc",
                        int bh2_lo=4, int bh2_hi=10,
                        double mipFrac=0.1,
                        double mipMeVperCm=2.0,
                        double BH2_thickness_mm=5.0,
                        double HTOF_thickness_mm=10.0,
                        bool save=true,
                        const char* tag="BH2_4_10")
{
  using namespace HTS;

  // open
  TFile* f=TFile::Open(filename,"READ");
  if(!f||f->IsZombie()){ std::cerr<<"[ERR] open "<<filename<<" failed\n"; return; }
  TTree* T=(TTree*)f->Get(treename);
  if(!T){ std::cerr<<"[ERR] tree "<<treename<<" not found\n"; return; }

  if(!T->GetBranch("BH2") || !T->GetBranch("HTOF")){
    std::cerr<<"[ERR] need BH2 & HTOF (vector<TParticle>)\n"; return;
  }
  const bool hasBH2edep  = (T->GetBranch("BH2_edep")  != nullptr);
  const bool hasHTOFedep = (T->GetBranch("HTOF_edep") != nullptr);

  std::vector<TParticle>* BH2=nullptr;   T->SetBranchAddress("BH2",&BH2);
  std::vector<TParticle>* HTOF=nullptr;  T->SetBranchAddress("HTOF",&HTOF);
  std::vector<double>* BH2_edep=nullptr;   if(hasBH2edep)  T->SetBranchAddress("BH2_edep",&BH2_edep);
  std::vector<double>* HTOF_edep=nullptr;  if(hasHTOFedep) T->SetBranchAddress("HTOF_edep",&HTOF_edep);

  // thresholds
  const double thrBH2  = MIPThrMeV(mipFrac,mipMeVperCm,BH2_thickness_mm);   // e.g. 0.10 MeV
  const double thrHTOF = MIPThrMeV(mipFrac,mipMeVperCm,HTOF_thickness_mm);  // e.g. 0.20 MeV
  std::cout<<std::fixed<<std::setprecision(3)
           <<"[INFO] BH2 thr="<<thrBH2<<" MeV, HTOF thr="<<thrHTOF<<" MeV | "
           <<"BH2 range = "<<bh2_lo<<"-"<<bh2_hi<<"\n";
  std::cout<<"[INFO] ID policy: HTOF=StatusCode, BH2=coord-map; edep: *_edep→Weight; valid if sum>=thr.\n";

  // counters
  Long64_t N_total=0, N_bh2In=0;
  Long64_t N_trig0=0; // (any BH2) & HTOF mult>=2  (참고용)
  Long64_t N_trig1=0; // (BH2 in range) & HTOF mult>=2
  Long64_t N_trig2=0; // trig1 & !beamveto

  // collect Trig2 survivors' HTOF combinations
  std::map<std::string, Long64_t> comboCounts;

  const Long64_t N=T->GetEntries();
  for(Long64_t ie=0; ie<N; ++ie){
    T->GetEntry(ie); N_total++;

    // --- BH2 valid (sum edep ≥ thr)
    std::map<int,double> bh2E;
    if(BH2){
      for(size_t i=0;i<BH2->size();++i){
        const TParticle& p=BH2->at(i);
        int sid = MapBH2_WorldToSeg(p.Vx(),p.Vy(),p.Vz());
        if(0<=sid && sid<kNBH2Seg){
          const double ed = (hasBH2edep && BH2_edep && i<BH2_edep->size()) ? BH2_edep->at(i)
                                                                            : p.GetWeight();
          bh2E[sid]+=ed;
        }
      }
    }
    std::set<int> bh2Valid;
    for(auto& kv:bh2E) if(kv.second>=thrBH2) bh2Valid.insert(kv.first);

    // --- HTOF valid (ID=StatusCode, sum edep ≥ thr)
    std::map<int,double> htofE;
    if(HTOF){
      for(size_t i=0;i<HTOF->size();++i){
        const TParticle& p=HTOF->at(i);
        int tid = p.GetStatusCode();
        if(0<=tid && tid<kNHTOF){
          const double ed = (hasHTOFedep && HTOF_edep && i<HTOF_edep->size()) ? HTOF_edep->at(i)
                                                                               : p.GetWeight();
          htofE[tid]+=ed;
        }
      }
    }
    std::set<int> htofValid;
    for(auto& kv:htofE) if(kv.second>=thrHTOF) htofValid.insert(kv.first);
    const int mult = (int)htofValid.size();

    // Trig0 (참고): any BH2? → 여기선 BH2 유효가 비어도 HTOF mult>=2면 세긴다 (원 정의 그대로)
    if(mult>=2) N_trig0++;

    // Range denom
    bool bh2_in_range=false;
    for(int s:bh2Valid){ if(bh2_lo<=s && s<=bh2_hi){ bh2_in_range=true; break; } }
    if(bh2_in_range) N_bh2In++;

    // Trig1
    if(bh2_in_range && mult>=2){
      N_trig1++;

      // BeamVeto
      const bool has20 = htofValid.count(20);
      const bool has21 = htofValid.count(21);
      const bool has22 = htofValid.count(22);
      const bool has23 = htofValid.count(23);
      const bool beamVeto = (has20 && has21) || (has22 && has23);

      if(!beamVeto){
        N_trig2++;
        // record combo
        const std::string key = JoinTuple(htofValid);
        comboCounts[key]++;
      }
    }
  }

  auto pct = [](LongLong_t a, LongLong_t b)->double{ return (b>0)? (100.0*(double)a/(double)b):0.0; };

  // ---- print summary ----
  std::cout<<"\n==== Trigger summary (BH2 "<<bh2_lo<<"-"<<bh2_hi<<") ====\n";
  std::cout<<"Total events             : "<<N_total<<"\n";
  std::cout<<"BH2 in range (denom)     : "<<N_bh2In<<"\n";
  std::cout<<"Trig0 (any BH2 & m>=2)   : "<<N_trig0<<"  ("<<std::setprecision(3)<<pct(N_trig0,N_total)<<" % of Total)\n";
  std::cout<<"Trig1 (range & m>=2)     : "<<N_trig1<<"  ("<<pct(N_trig1,N_bh2In)<<" % of BH2-range)\n";
  std::cout<<"Trig2 (Trig1 & !Veto)    : "<<N_trig2<<"  ("<<pct(N_trig2,N_bh2In)<<" % of BH2-range)\n";

  // ---- build histogram of Trig2 HTOF combos (dynamic bins) ----
  const int nComb = (int)comboCounts.size();
  TH1I* hComb = new TH1I("hHTOFComb_trig2",
                         Form("HTOF combinations (Trig2 survivors) | BH2 %d-%d;Combo;Events",
                              bh2_lo,bh2_hi),
                         nComb, 0.5, nComb+0.5);
  hComb->SetDirectory(nullptr);

  int bin=1;
  for(const auto& kv:comboCounts){
    hComb->GetXaxis()->SetBinLabel(bin, kv.first.c_str());
    hComb->SetBinContent(bin, (double)kv.second);
    ++bin;
  }

  // draw/save
  TCanvas* c=new TCanvas("cTrig2Comb","Trig2 HTOF combinations", std::max(900, 18*std::max(1,nComb)), 600);
  hComb->LabelsOption("v","X");
  hComb->Draw("hist");
  if(save){
    TString t(tag);
    c->SaveAs("Trig2_HTOFcombo_"+t+".png");
  }
}

// -*- C++ -*-
// Trigger_study_plusbeam_10_28.C  (2025-10-28 for jaejin)
// Sections:
//   (Sec1) BH2 Seg 4–10
//   (Sec2) Narrow Beam(1): BH2 Seg 4–9
//   (Sec3) Narrow Beam(2): BH2 Seg 5–10
//
// Trigger logic per section:
//   Beam      = (BH2 in [lo,hi])                         // 분모
//   Trig1     = Beam && (HTOF multiplicity >= 2)
//   Trig2(*)  = Trig1 && !(BeamVeto(*))
//
// BeamVeto variants for **plus-beam** (HTOF tile copy-no.):
//   tight      : (17,18)&&(18,19)              → {17,18,19}   모두 존재하면 veto
//   fit        : tight + (19,20)               → {17,18,19,20}
//   wide       : fit   + (16,17)               → {16,17,18,19,20}
//   ultra-wide : wide  + (20,21)               → {16,17,18,19,20,21}
//
// Denominator for % : #events with BH2 in [lo,hi]  (Beam)
//
// Cuts: default 0.1 MIP with dE/dx=2 MeV/cm, BH2=5 mm, HTOF=10 mm.

#include "TFile.h"
#include "TTree.h"
#include "TParticle.h"
#include "TSystem.h"
#include "TInterpreter.h"
#include "TString.h"
#include <vector>
#include <map>
#include <set>
#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>

namespace TSP {

// ---------- geometry (E72-like) ----------
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

// mm → MeV threshold
static inline double MIPThrMeV(double frac,double dEdx,double thk_mm){
  return frac * dEdx * (thk_mm*0.1); // mm→cm
}

// x(world) → BH2 seg id
static inline int MapBH2_WorldToSeg(double x, double, double){
  const double xloc = x - kBH2_x0;
  const double xmin = -0.5*kBH2_total, xmax = 0.5*kBH2_total;
  if (xloc < xmin || xloc >= xmax) return -1;
  int idx = int(std::floor((xloc - xmin)/kBH2_pitch));
  if(idx<0) idx=0; if(idx>=kNBH2Seg) idx=kNBH2Seg-1;
  return idx;
}

// pretty %
static inline double pct(long long a, long long b){
  return (b>0)? (100.0*double(a)/double(b)):0.0;
}

// BeamVeto variants (based on presence of tile IDs)
// plus-beam spec: tight={17,18,19}, fit={17,18,19,20}, wide={16..20}, ultra={16..21}
struct VetoSet { bool tight=false, fit=false, wide=false, ultra=false; };

static inline VetoSet EvalVeto(const std::set<int>& tiles){
  auto has=[&](int t){ return tiles.count(t)>0; };
  const bool h16=has(16), h17=has(17), h18=has(18), h19=has(19), h20=has(20), h21=has(21);

  VetoSet v;
  v.tight = (h17 && h18 && h19);
  v.fit   = (h17 && h18 && h19 && h20);
  v.wide  = (h16 && h17 && h18 && h19 && h20);
  v.ultra = (h16 && h17 && h18 && h19 && h20 && h21);
  return v;
}

// One pass over the tree for a given BH2 range
struct Counts {
  long long N_total=0;
  long long N_beam =0;     // BH2 in [lo,hi]
  long long N_trig1=0;     // Beam && MP>=2
  long long N_trig2_tight=0, N_trig2_fit=0, N_trig2_wide=0, N_trig2_ultra=0; // survivors (!veto)
};

static Counts ProcessRange(TTree* T,
                           int bh2_lo, int bh2_hi,
                           double thrBH2, double thrHTOF,
                           bool hasBH2edep, bool hasHTOFedep)
{
  std::vector<TParticle>* BH2=nullptr;   T->SetBranchAddress("BH2",&BH2);
  std::vector<TParticle>* HTOF=nullptr;  T->SetBranchAddress("HTOF",&HTOF);
  std::vector<double>* BH2_edep=nullptr; if(hasBH2edep)  T->SetBranchAddress("BH2_edep",&BH2_edep);
  std::vector<double>* HTOF_edep=nullptr;if(hasHTOFedep) T->SetBranchAddress("HTOF_edep",&HTOF_edep);

  Counts C;
  const Long64_t N=T->GetEntries();
  for(Long64_t ie=0; ie<N; ++ie){
    T->GetEntry(ie); C.N_total++;

    // ---- BH2 valid segments (Σedep≥thr)
    std::map<int,double> bh2E;
    if(BH2){
      for(size_t i=0;i<BH2->size();++i){
        const TParticle& p = BH2->at(i);
        int sid = MapBH2_WorldToSeg(p.Vx(),p.Vy(),p.Vz());
        if(0<=sid && sid<kNBH2Seg){
          const double ed = (hasBH2edep && BH2_edep && i<BH2_edep->size()) ? BH2_edep->at(i) : p.GetWeight();
          bh2E[sid]+=ed;
        }
      }
    }
    std::set<int> bh2Valid;
    for(const auto& kv:bh2E) if(kv.second>=thrBH2) bh2Valid.insert(kv.first);

    bool beam=false;
    for(int s:bh2Valid){ if(bh2_lo<=s && s<=bh2_hi){ beam=true; break; } }
    if(beam) C.N_beam++;

    // ---- HTOF valid tiles (ID=StatusCode, Σedep≥thr)
    std::map<int,double> htofE;
    if(HTOF){
      for(size_t i=0;i<HTOF->size();++i){
        const TParticle& p = HTOF->at(i);
        int tid = p.GetStatusCode();
        if(0<=tid && tid<kNHTOF){
          const double ed = (hasHTOFedep && HTOF_edep && i<HTOF_edep->size()) ? HTOF_edep->at(i) : p.GetWeight();
          htofE[tid]+=ed;
        }
      }
    }
    std::set<int> tiles;
    for(const auto& kv:htofE) if(kv.second>=thrHTOF) tiles.insert(kv.first);
    const int mult = (int)tiles.size();

    // ---- Trig1 & Trig2(*)
    if(beam && mult>=2){
      C.N_trig1++;
      const auto v = EvalVeto(tiles);
      if(!v.tight) C.N_trig2_tight++;
      if(!v.fit)   C.N_trig2_fit++;
      if(!v.wide)  C.N_trig2_wide++;
      if(!v.ultra) C.N_trig2_ultra++;
    }
  }
  return C;
}

// Extra KPIs per section:
//  (1) 모든 Trig의 Beam 기준 비율 [%]
//  (2) Beam Veto efficiency = 100 - 100*(Trig2/Beam)  [소수점 3자리]
//  (3) Beam induced background (counts) = round(1e6 * Trig2/Beam)
static void PrintSection(const char* title, const Counts& C, int lo, int hi){
  const double BeamNorm = 1e6; // Beam=1,000,000 가정
  auto bg = [&](long long trig2)->long long{
    if(C.N_beam<=0) return 0;
    return (long long) llround( BeamNorm * (double)trig2 / (double)C.N_beam );
  };
  auto p = [&](long long x)->double{
    return pct(x, C.N_beam);
  };
  auto eff = [&](long long trig2)->double{
    return 100.0 - p(trig2);  // Beam veto efficiency (%)
  };

  std::cout<<"\n== "<<title<<" | BH2 "<<lo<<"–"<<hi<<" ==\n";
  std::cout<<"Total events                 : "<<C.N_total<<"\n";
  std::cout<<"Beam  (BH2 in-range)         : "<<C.N_beam<<"\n";

  // (1) Beam 기준 비율
  std::cout<<std::fixed;
  std::cout<<"[Beam-based %]\n";
  std::cout<<"  Trig1 (Beam && MP>=2)      : "<<std::setprecision(3)<<p(C.N_trig1)<<" %\n";
  std::cout<<"  Trig2 tight  (!veto)       : "<<std::setprecision(3)<<p(C.N_trig2_tight)<<" %\n";
  std::cout<<"  Trig2 fit    (!veto)       : "<<std::setprecision(3)<<p(C.N_trig2_fit)  <<" %\n";
  std::cout<<"  Trig2 wide   (!veto)       : "<<std::setprecision(3)<<p(C.N_trig2_wide) <<" %\n";
  std::cout<<"  Trig2 ultra  (!veto)       : "<<std::setprecision(3)<<p(C.N_trig2_ultra)<<" %\n";

  // (2) Beam Veto efficiency (소수점 3자리)
  std::cout<<"[Beam Veto efficiency (%), 100 - Trig2/Beam*100]\n";
  std::cout<<"  tight                      : "<<std::setprecision(3)<<eff(C.N_trig2_tight)<<"\n";
  std::cout<<"  fit                        : "<<std::setprecision(3)<<eff(C.N_trig2_fit)  <<"\n";
  std::cout<<"  wide                       : "<<std::setprecision(3)<<eff(C.N_trig2_wide) <<"\n";
  std::cout<<"  ultra                      : "<<std::setprecision(3)<<eff(C.N_trig2_ultra)<<"\n";

  // (3) Beam-induced background (Beam=1,000,000 가정)
  std::cout<<"[Beam-induced background (counts) @ Beam=1,000,000]\n";
  std::cout<<"  tight                      : "<<bg(C.N_trig2_tight)<<"\n";
  std::cout<<"  fit                        : "<<bg(C.N_trig2_fit)  <<"\n";
  std::cout<<"  wide                       : "<<bg(C.N_trig2_wide) <<"\n";
  std::cout<<"  ultra                      : "<<bg(C.N_trig2_ultra)<<"\n";
}

} // namespace TSP


void Trigger_study_plusbeam_10_28(const char* filename,
                                  const char* treename="g4hyptpc",
                                  // Section 1 default BH2 4–10
                                  int bh2_lo=4, int bh2_hi=10,
                                  double mipFrac=0.10,
                                  double mipMeVperCm=2.0,
                                  double BH2_thickness_mm=5.0,
                                  double HTOF_thickness_mm=10.0)
{
  using namespace TSP;

  // open & branches
  TFile* f=TFile::Open(filename,"READ");
  if(!f||f->IsZombie()){ std::cerr<<"[ERR] open "<<filename<<" failed\n"; return; }
  TTree* T=(TTree*)f->Get(treename);
  if(!T){ std::cerr<<"[ERR] tree "<<treename<<" not found\n"; return; }
  if(!T->GetBranch("BH2") || !T->GetBranch("HTOF")){
    std::cerr<<"[ERR] need BH2 & HTOF (vector<TParticle>) branches\n"; return;
  }
  const bool hasBH2edep  = (T->GetBranch("BH2_edep")  != nullptr);
  const bool hasHTOFedep = (T->GetBranch("HTOF_edep") != nullptr);

  // thresholds
  const double thrBH2  = MIPThrMeV(mipFrac,mipMeVperCm,BH2_thickness_mm);   // e.g. 0.10 MeV
  const double thrHTOF = MIPThrMeV(mipFrac,mipMeVperCm,HTOF_thickness_mm);  // e.g. 0.20 MeV
  std::cout<<std::fixed<<std::setprecision(3)
           <<"[INFO] BH2 thr="<<thrBH2<<" MeV, HTOF thr="<<thrHTOF<<" MeV\n";
  std::cout<<"[INFO] HTOF ID = StatusCode, BH2 ID = coord→seg; edep: *_edep else Weight; VALID if sum>=thr.\n";

  // -------- Section 1: BH2 4–10 --------
  const int s1_lo = bh2_lo, s1_hi = bh2_hi;         // default 4–10
  // -------- Section 2: BH2 4–9  --------
  const int s2_lo = 4,      s2_hi = 9;
  // -------- Section 3: BH2 5–10 --------
  const int s3_lo = 5,      s3_hi = 10;

  auto C1 = ProcessRange(T, s1_lo, s1_hi, thrBH2, thrHTOF, hasBH2edep, hasHTOFedep);
  auto C2 = ProcessRange(T, s2_lo, s2_hi, thrBH2, thrHTOF, hasBH2edep, hasHTOFedep);
  auto C3 = ProcessRange(T, s3_lo, s3_hi, thrBH2, thrHTOF, hasBH2edep, hasHTOFedep);

  PrintSection("Section 1 (Full Beam)",   C1, s1_lo, s1_hi);
  PrintSection("Section 2 (Narrow 1)",    C2, s2_lo, s2_hi);
  PrintSection("Section 3 (Narrow 2)",    C3, s3_lo, s3_hi);
}

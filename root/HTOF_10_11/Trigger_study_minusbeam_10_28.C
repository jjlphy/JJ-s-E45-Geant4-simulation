// -*- C++ -*-
// Trigger_study_minusbeam_10_28.C (2025-10-28 for jaejin)
// Minus-beam version using veto windows centered on tiles 20–24.
//
// Sections:
//   Sec1: BH2 4–10   | Sec2: BH2 4–9   | Sec3: BH2 5–10
//
// Trigger logic per section:
//   Beam   = (BH2 in [lo,hi])
//   Trig1  = Beam && (HTOF multiplicity >= 2)              // after excluding tiles
//   Trig2* = Trig1 && !BeamVeto(*)
//
// BeamVeto variants (copy-no):
//   tight      : (20,21) && (21,22)         → {20,21,22} all present
//   fit        : tight   && (22,23)         → {20,21,22,23}
//   wide       : fit     && (19,20)         → {19,20,21,22,23}
//   ultra-wide : wide    && (23,24)         → {19,20,21,22,23,24}
//
// Denominator for % : #events with BH2 in [lo,hi] (Beam)
//
// Extras:
//  - Prints ABS counts + % (Beam-based), Beam Veto efficiency, Beam-induced background (@Beam=1e6)
//  - Consistency check: Trig1 == Trig2_survivors + Veto_fired (per variant)
//  - excludeTilesCSV: comma list (e.g., "20,21") to exclude tiles from multiplicity & veto
//
// Usage:
//   root -l
//   .L Trigger_study_minusbeam_10_28.C+
//   Trigger_study_minusbeam_10_28("../rootfile/E45_Nov_beamminus_098.root","g4hyptpc", 4,10, 0.10,2.0, 5.0,10.0, "");
//                                 "g4hyptpc", 4,10, 0.10,2.0, 5.0,10.0, "");
//
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
#include <sstream>
#include <cmath>
#include <algorithm>

namespace TSM {

static const int    kNBH2Seg   = 15;     // 0..14
static const double kBH2_x0    = -10.0;  // [mm]
static const double kBH2_pitch = 14.0;   // [mm]
static const double kBH2_total = kNBH2Seg * kBH2_pitch;
static const int    kNHTOF     = 34;     // 0..33

struct DictGuard {
  DictGuard(){ gSystem->Load("libPhysics");
               gInterpreter->GenerateDictionary("vector<TParticle>","TParticle.h;vector"); }
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
static inline double pct(long long a, long long b){ return (b>0)? (100.0*double(a)/double(b)) : 0.0; }

// parse "a,b,c" → {a,b,c}
static std::set<int> ParseCSVInt(const char* csv){
  std::set<int> s; if(!csv) return s;
  std::stringstream ss(csv); std::string tok;
  while(std::getline(ss, tok, ',')){ if(tok.empty()) continue; s.insert(std::stoi(tok)); }
  return s;
}

// ---- minus-beam veto centered on 20–24 ----
struct VetoSet { bool tight=false, fit=false, wide=false, ultra=false; };
static inline VetoSet EvalVeto(const std::set<int>& tiles){
  auto has=[&](int t){ return tiles.count(t)>0; };
  const bool h19=has(19), h20=has(20), h21=has(21), h22=has(22), h23=has(23), h24=has(24);

  VetoSet v;
  // tight: (20,21)&&(21,22) → {20,21,22}
  v.tight = (h20 && h21 && h22);
  // fit  : tight&&(22,23)   → {20,21,22,23}
  v.fit   = v.tight && h23;
  // wide : fit&&(19,20)     → {19,20,21,22,23}
  v.wide  = v.fit && h19;
  // ultra: wide&&(23,24)    → {19,20,21,22,23,24}
  v.ultra = v.wide && h24;
  return v;
}

struct Counts {
  long long N_total=0;
  long long N_beam =0;     // BH2 in [lo,hi]
  long long N_trig1=0;     // Beam && MP>=2 (after exclusion)
  // survivors (!veto)
  long long N_trig2_tight=0, N_trig2_fit=0, N_trig2_wide=0, N_trig2_ultra=0;
  // veto fired counts
  long long N_veto_tight=0, N_veto_fit=0, N_veto_wide=0, N_veto_ultra=0;
  // tile presence (optional)
  long long have19=0, have20=0, have21=0, have22=0, have23=0, have24=0;
};

static Counts ProcessRange(TTree* T,
                           int bh2_lo, int bh2_hi,
                           double thrBH2, double thrHTOF,
                           bool hasBH2edep, bool hasHTOFedep,
                           const std::set<int>& exclude)
{
  std::vector<TParticle>* BH2=nullptr;   T->SetBranchAddress("BH2",&BH2);
  std::vector<TParticle>* HTOF=nullptr;  T->SetBranchAddress("HTOF",&HTOF);
  std::vector<double>* BH2_edep=nullptr; if(hasBH2edep)  T->SetBranchAddress("BH2_edep",&BH2_edep);
  std::vector<double>* HTOF_edep=nullptr;if(hasHTOFedep) T->SetBranchAddress("HTOF_edep",&HTOF_edep);

  Counts C;
  const Long64_t N=T->GetEntries();
  for(Long64_t ie=0; ie<N; ++ie){
    T->GetEntry(ie); C.N_total++;

    // ---- BH2 valid (Σedep≥thr) ----
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
    if(beam) C.N_beam++;

    // ---- HTOF valid (Σedep≥thr) with exclusion ----
    std::map<int,double> htofE;
    if(HTOF){
      for(size_t i=0;i<HTOF->size();++i){
        const TParticle& p=HTOF->at(i);
        int tid = p.GetStatusCode();
        if(0<=tid && tid<kNHTOF){
          if(exclude.count(tid)) continue;
          const double ed=(hasHTOFedep && HTOF_edep && i<HTOF_edep->size())? HTOF_edep->at(i) : p.GetWeight();
          htofE[tid]+=ed;
        }
      }
    }
    std::set<int> tiles;
    for(const auto& kv:htofE) if(kv.second>=thrHTOF) tiles.insert(kv.first);
    const int mult = (int)tiles.size();

    // presence (diagnostic)
    if(tiles.count(19)) C.have19++;
    if(tiles.count(20)) C.have20++;
    if(tiles.count(21)) C.have21++;
    if(tiles.count(22)) C.have22++;
    if(tiles.count(23)) C.have23++;
    if(tiles.count(24)) C.have24++;

    // ---- Trig1 & Trig2(*)
    if(beam && mult>=2){
      C.N_trig1++;
      const auto v = EvalVeto(tiles);
      if(!v.tight) C.N_trig2_tight++; else C.N_veto_tight++;
      if(!v.fit)   C.N_trig2_fit++;   else C.N_veto_fit++;
      if(!v.wide)  C.N_trig2_wide++;  else C.N_veto_wide++;
      if(!v.ultra) C.N_trig2_ultra++; else C.N_veto_ultra++;
    }
  }
  return C;
}

static void PrintSection(const char* title, const Counts& C, int lo, int hi){
  const double BeamNorm = 1e6; // Beam=1,000,000 가정
  auto bg = [&](long long trig2)->long long{
    if(C.N_beam<=0) return 0;
    return (long long) llround( BeamNorm * (double)trig2 / (double)C.N_beam );
  };
  auto p   = [&](long long x)->double{ return pct(x, C.N_beam); };
  auto eff = [&](long long trig2)->double{ return 100.0 - p(trig2); };

  std::cout<<"\n== "<<title<<" | BH2 "<<lo<<"–"<<hi<<" ==\n";
  std::cout<<"Total events                 : "<<C.N_total<<"\n";
  std::cout<<"Beam  (BH2 in-range)         : "<<C.N_beam<<"\n";

  std::cout<<std::fixed<<std::setprecision(3);
  std::cout<<"[Beam-based counts & %]\n";
  std::cout<<"  Trig1 (Beam && MP>=2)      : "<<C.N_trig1<<"  ("<<p(C.N_trig1)<<" %)\n";
  std::cout<<"  Trig2 tight  (!veto)       : "<<C.N_trig2_tight<<"  ("<<p(C.N_trig2_tight)<<" %)\n";
  std::cout<<"  Trig2 fit    (!veto)       : "<<C.N_trig2_fit  <<"  ("<<p(C.N_trig2_fit)  <<" %)\n";
  std::cout<<"  Trig2 wide   (!veto)       : "<<C.N_trig2_wide <<"  ("<<p(C.N_trig2_wide) <<" %)\n";
  std::cout<<"  Trig2 ultra  (!veto)       : "<<C.N_trig2_ultra<<"  ("<<p(C.N_trig2_ultra)<<" %)\n";

  std::cout<<"[Beam Veto efficiency (%), = 100 - Trig2/Beam*100]\n";
  std::cout<<"  tight                      : "<<eff(C.N_trig2_tight)<<"\n";
  std::cout<<"  fit                        : "<<eff(C.N_trig2_fit)  <<"\n";
  std::cout<<"  wide                       : "<<eff(C.N_trig2_wide) <<"\n";
  std::cout<<"  ultra                      : "<<eff(C.N_trig2_ultra)<<"\n";

  std::cout<<"[Beam-induced background (counts) @ Beam=1,000,000]\n";
  std::cout<<"  tight                      : "<<bg(C.N_trig2_tight)<<"\n";
  std::cout<<"  fit                        : "<<bg(C.N_trig2_fit)  <<"\n";
  std::cout<<"  wide                       : "<<bg(C.N_trig2_wide) <<"\n";
  std::cout<<"  ultra                      : "<<bg(C.N_trig2_ultra)<<"\n";

  auto chk = [&](const char* name, long long surv, long long veto){
    const long long sum = surv + veto;
    std::cout<<"[CHECK] "<<name<<": Trig1 == survivors + veto ?  "
             <<C.N_trig1<<" == "<<surv<<" + "<<veto<<"  → "
             <<((C.N_trig1==sum) ? "OK" : "MISMATCH!")<<"\n";
  };
  std::cout<<"[Veto fired counts]\n";
  std::cout<<"  veto_tight                 : "<<C.N_veto_tight<<"\n";
  std::cout<<"  veto_fit                   : "<<C.N_veto_fit  <<"\n";
  std::cout<<"  veto_wide                  : "<<C.N_veto_wide <<"\n";
  std::cout<<"  veto_ultra                 : "<<C.N_veto_ultra<<"\n";
  chk("tight", C.N_trig2_tight, C.N_veto_tight);
  chk("fit",   C.N_trig2_fit,   C.N_veto_fit);
  chk("wide",  C.N_trig2_wide,  C.N_veto_wide);
  chk("ultra", C.N_trig2_ultra, C.N_veto_ultra);

  std::cout<<"[Tile presence among Beam events (approx, for diagnosis)]\n";
  std::cout<<"  have19="<<C.have19<<", have20="<<C.have20<<", have21="<<C.have21
           <<", have22="<<C.have22<<", have23="<<C.have23<<", have24="<<C.have24<<"\n";
}

} // namespace TSM


void Trigger_study_minusbeam_10_28(const char* filename,
                                   const char* treename="g4hyptpc",
                                   int bh2_lo=4, int bh2_hi=10,
                                   double mipFrac=0.10,
                                   double mipMeVperCm=2.0,
                                   double BH2_thickness_mm=5.0,
                                   double HTOF_thickness_mm=10.0,
                                   const char* excludeTilesCSV="" )
{
  using namespace TSM;

  // open
  TFile* f=TFile::Open(filename,"READ");
  if(!f||f->IsZombie()){ std::cerr<<"[ERR] open "<<filename<<" failed\n"; return; }
  TTree* T=(TTree*)f->Get(treename);
  if(!T){ std::cerr<<"[ERR] tree "<<treename<<" not found\n"; return; }
  if(!T->GetBranch("BH2") || !T->GetBranch("HTOF")){
    std::cerr<<"[ERR] need BH2 & HTOF (vector<TParticle>) branches\n"; return;
  }
  const bool hasBH2edep  = (T->GetBranch("BH2_edep")  != nullptr);
  const bool hasHTOFedep = (T->GetBranch("HTOF_edep") != nullptr);

  // thresholds & exclude
  const double thrBH2  = MIPThrMeV(mipFrac,mipMeVperCm,BH2_thickness_mm);
  const double thrHTOF = MIPThrMeV(mipFrac,mipMeVperCm,HTOF_thickness_mm);
  auto exclude = ParseCSVInt(excludeTilesCSV);

  std::cout<<std::fixed<<std::setprecision(3)
           <<"[INFO] BH2 thr="<<thrBH2<<" MeV, HTOF thr="<<thrHTOF<<" MeV\n";
  if(!exclude.empty()){
    std::cout<<"[INFO] Excluding HTOF tiles: ";
    bool first=true; for(int t:exclude){ if(!first) std::cout<<","; first=false; std::cout<<t; }
    std::cout<<"\n";
  }

  // sections
  const int s1_lo = bh2_lo, s1_hi = bh2_hi; // Full Beam
  const int s2_lo = 4,      s2_hi = 9;      // Narrow 1
  const int s3_lo = 5,      s3_hi = 10;     // Narrow 2

  auto C1 = TSM::ProcessRange(T, s1_lo, s1_hi, thrBH2, thrHTOF, hasBH2edep, hasHTOFedep, exclude);
  auto C2 = TSM::ProcessRange(T, s2_lo, s2_hi, thrBH2, thrHTOF, hasBH2edep, hasHTOFedep, exclude);
  auto C3 = TSM::ProcessRange(T, s3_lo, s3_hi, thrBH2, thrHTOF, hasBH2edep, hasHTOFedep, exclude);

  TSM::PrintSection("Section 1 (Full Beam)", C1, s1_lo, s1_hi);
  TSM::PrintSection("Section 2 (Narrow 1)",  C2, s2_lo, s2_hi);
  TSM::PrintSection("Section 3 (Narrow 2)",  C3, s3_lo, s3_hi);
}

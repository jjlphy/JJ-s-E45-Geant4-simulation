// -*- C++ -*-
// Trigger_study_minusreaction_10_28.C  (2025-10-28 for jaejin; minus-beam OR-pair veto)
//
// Veto (OR over adjacent HTOF pairs; copy-no):
//   tight      : (20,21) OR (21,22)
//   fit        : tight    OR (22,23)
//   wide       : fit      OR (19,20)
//   ultra-wide : wide     OR (23,24)
//
// Sections:
//   Sec1: BH2 4–10      (Full Beam)
//   Sec2: BH2 4–9       (Narrow 1)
//   Sec3: BH2 5–10      (Narrow 2)
//
// Denominator for % : Beam = #events with BH2 in [lo,hi]
// Overkill(%) = 100 × (BeamVeto / Beam)
// Notes:
//   - HTOF tile ID = TParticle::StatusCode()  (0..33)
//   - BH2 seg from world (x) via E72-like mapping
//   - VALID hit: sum(edep) ≥ threshold; energy = *_edep if exists else Weight()
//   - Optional analysis-side exclusion: excludeTilesCSV (e.g. "20,21")
//
// Usage:
//   root -l
//   .L Trigger_study_minusreaction_10_28.C+
//   Trigger_study_minusreaction_10_28("../rootfile/E45_Nov_beamminus_105.root","g4hyptpc",
//                                     4,10, 0.10,2.0, 5.0,10.0, "");
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

namespace TSMR {

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

static inline double pct(long long a, long long b){
  return (b>0)? (100.0 * double(a) / double(b)) : 0.0;
}

// parse "a,b,c" → {a,b,c}
static std::set<int> ParseCSVInt(const char* csv){
  std::set<int> s; if(!csv) return s;
  std::stringstream ss(csv); std::string tok;
  while(std::getline(ss, tok, ',')){
    if(tok.empty()) continue;
    s.insert(std::stoi(tok));
  }
  return s;
}

// ---- Veto variants (minus-beam OR-pair logic) ----
struct VetoFlags { bool tight=false, fit=false, wide=false, ultra=false; };

static inline VetoFlags EvalVetoMinus_ORpairs(const std::set<int>& tiles){
  const bool p19_20 = PairHit(tiles,19,20);
  const bool p20_21 = PairHit(tiles,20,21);
  const bool p21_22 = PairHit(tiles,21,22);
  const bool p22_23 = PairHit(tiles,22,23);
  const bool p23_24 = PairHit(tiles,23,24);

  VetoFlags v;
  v.tight = (p20_21 || p21_22);
  v.fit   = (v.tight || p22_23);
  v.wide  = (v.fit   || p19_20);
  v.ultra = (v.wide  || p23_24);
  return v;
}

struct Counts {
  long long N_total=0;
  long long N_beam =0;      // BH2 in [lo,hi]
  long long N_trig1=0;      // Beam && (HTOF m≥2)  [참고용]

  // Veto fired (within Beam; 실무적으로 m≥2가 아니면 인접쌍이 성립X이므로
  // 아래 카운트는 Beam && m≥2 조건 위에서 평가됨)
  long long veto_tight=0, veto_fit=0, veto_wide=0, veto_ultra=0;

  // 진단용: 해당 페어 존재 여부 누적(Beam&&m≥2에서 집계)
  long long have19_20=0, have20_21=0, have21_22=0, have22_23=0, have23_24=0;
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
  const Long64_t N = T->GetEntries();
  for(Long64_t ie=0; ie<N; ++ie){
    T->GetEntry(ie); C.N_total++;

    // ---- BH2: 유효 세그먼트
    std::map<int,double> bh2E;
    if(BH2){
      for(size_t i=0;i<BH2->size();++i){
        const TParticle& p = BH2->at(i);
        int sid = MapBH2_WorldToSeg(p.Vx(),p.Vy(),p.Vz());
        if(0<=sid && sid<kNBH2Seg){
          const double ed = (hasBH2edep && BH2_edep && i<BH2_edep->size())
                            ? BH2_edep->at(i) : p.GetWeight();
          bh2E[sid]+=ed;
        }
      }
    }
    std::set<int> bh2Valid;
    for(const auto& kv:bh2E) if(kv.second>=thrBH2) bh2Valid.insert(kv.first);

    bool beam=false;
    for(int s:bh2Valid){ if(bh2_lo<=s && s<=bh2_hi){ beam=true; break; } }
    if(beam) C.N_beam++;

    // ---- HTOF: 유효 타일 (배제 적용)
    std::map<int,double> htofE;
    if(HTOF){
      for(size_t i=0;i<HTOF->size();++i){
        const TParticle& p = HTOF->at(i);
        int tid = p.GetStatusCode();
        if(0<=tid && tid<kNHTOF){
          if(exclude.count(tid)) continue;
          const double ed = (hasHTOFedep && HTOF_edep && i<HTOF_edep->size())
                            ? HTOF_edep->at(i) : p.GetWeight();
          htofE[tid]+=ed;
        }
      }
    }
    std::set<int> tiles;
    for(const auto& kv:htofE) if(kv.second>=thrHTOF) tiles.insert(kv.first);
    const int mult = (int)tiles.size();

    // ---- Trig1 (참고): Beam && m≥2
    if(beam && mult>=2) C.N_trig1++;

    // ---- BeamVeto OR-pair 평가 (Beam && m≥2 위에서만 의미 있음)
    if(beam && mult>=2){
      // 페어 존재 통계(진단)
      const bool p19_20 = PairHit(tiles,19,20);
      const bool p20_21 = PairHit(tiles,20,21);
      const bool p21_22 = PairHit(tiles,21,22);
      const bool p22_23 = PairHit(tiles,22,23);
      const bool p23_24 = PairHit(tiles,23,24);
      if(p19_20) C.have19_20++;
      if(p20_21) C.have20_21++;
      if(p21_22) C.have21_22++;
      if(p22_23) C.have22_23++;
      if(p23_24) C.have23_24++;

      const VetoFlags v = EvalVetoMinus_ORpairs(tiles);
      if(v.tight) C.veto_tight++;
      if(v.fit)   C.veto_fit++;
      if(v.wide)  C.veto_wide++;
      if(v.ultra) C.veto_ultra++;
    }
  }
  return C;
}

static void PrintSection(const char* title, const Counts& C, int lo, int hi){
  std::cout<<"\n== "<<title<<" | BH2 "<<lo<<"–"<<hi<<" ==\n";
  std::cout<<"Total events                : "<<C.N_total<<"\n";
  std::cout<<"Beam  (BH2 in-range)        : "<<C.N_beam<<"\n";
  std::cout<<"Trig1 (Beam && m>=2)        : "<<C.N_trig1<<"\n";

  // Beam 기준 절대개수 + %
  auto pbeam = [&](long long x){ return pct(x, C.N_beam); };

  std::cout<<std::fixed<<std::setprecision(3);
  std::cout<<"[BeamVeto fired | Beam-based counts & %]\n";
  std::cout<<"  tight      : "<<C.veto_tight<<"  ("<<pbeam(C.veto_tight)<<" % of Beam)\n";
  std::cout<<"  fit        : "<<C.veto_fit  <<"  ("<<pbeam(C.veto_fit)  <<" % of Beam)\n";
  std::cout<<"  wide       : "<<C.veto_wide <<"  ("<<pbeam(C.veto_wide) <<" % of Beam)\n";
  std::cout<<"  ultra      : "<<C.veto_ultra<<"  ("<<pbeam(C.veto_ultra)<<" % of Beam)\n";

  // Overkill(%) = 100 × BeamVeto / Beam  (요청 정의)
  std::cout<<"[Overkill (%) = 100 × BeamVeto / Beam]\n";
  std::cout<<"  tight      : "<<pbeam(C.veto_tight)<<"\n";
  std::cout<<"  fit        : "<<pbeam(C.veto_fit)  <<"\n";
  std::cout<<"  wide       : "<<pbeam(C.veto_wide) <<"\n";
  std::cout<<"  ultra      : "<<pbeam(C.veto_ultra)<<"\n";

  // 진단(선택적): 각 인접 페어가 Beam&&m≥2에서 얼마나 자주 나타났는지
  std::cout<<"[Diag] Pair presence within (Beam && m>=2):\n";
  std::cout<<"  (19,20)="<<C.have19_20<<", (20,21)="<<C.have20_21
           <<", (21,22)="<<C.have21_22<<", (22,23)="<<C.have22_23
           <<", (23,24)="<<C.have23_24<<"\n";
}

} // namespace TSMR


void Trigger_study_minusreaction_10_28(const char* filename,
                                       const char* treename="g4hyptpc",
                                       int bh2_lo=4, int bh2_hi=10,
                                       double mipFrac=0.10,
                                       double mipMeVperCm=2.0,
                                       double BH2_thickness_mm=5.0,
                                       double HTOF_thickness_mm=10.0,
                                       const char* excludeTilesCSV="" ) // e.g. "20,21"
{
  using namespace TSMR;

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
  const int s1_lo = bh2_lo, s1_hi = bh2_hi; // Full
  const int s2_lo = 4,      s2_hi = 9;      // Narrow 1
  const int s3_lo = 5,      s3_hi = 10;     // Narrow 2

  auto C1 = ProcessRange(T, s1_lo, s1_hi, thrBH2, thrHTOF, hasBH2edep, hasHTOFedep, exclude);
  auto C2 = ProcessRange(T, s2_lo, s2_hi, thrBH2, thrHTOF, hasBH2edep, hasHTOFedep, exclude);
  auto C3 = ProcessRange(T, s3_lo, s3_hi, thrBH2, thrHTOF, hasBH2edep, hasHTOFedep, exclude);

  PrintSection("Section 1 (Full Beam)",   C1, s1_lo, s1_hi);
  PrintSection("Section 2 (Narrow 1)",    C2, s2_lo, s2_hi);
  PrintSection("Section 3 (Narrow 2)",    C3, s3_lo, s3_hi);
}

// -*- C++ -*-
// Trigger_study_minusreaction_10_28_HTOFselect.C  (2025-10-28 for jaejin; minus-beam OR-pair veto)
//   - MP>=2 판정에서만 HTOF {1,2,3,4}를 제외 (veto/페어 판정은 full-valid로 유지)
//   - "배제 때문에 잃는" 이벤트를 2종 집계:
//       (A) Lost_MP:  Beam && (MP_base>=2) && (MP_excl<2)
//       (B) Lost_Trig2_*: Beam && (Trig2_base_*==true) && (MP_excl<2)  // *_=tight/fit/wide/ultra
//   - 스케일: Beam(BH2 4–10)=1,000,000 으로 환산한 개수(각 섹션에 표기)
//
// Veto (OR over adjacent HTOF pairs; copy-no):
//   tight      : (20,21) OR (21,22)
//   fit        : tight    OR (22,23)
//   wide       : fit      OR (19,20)
//   ultra-wide : wide     OR (23,24)
//
// Sections:
//   Sec1: BH2 4–10      (Full Beam)   [MP-only exclusion = {1,2,3,4}]
//   Sec2: BH2 4–9       (Narrow 1)    [MP-only exclusion = {1,2,3,4}]
//
//   Sec3: BH2 5–10      (Narrow 2)    [MP-only exclusion = {1,2,3,4}]
//
// Denominator for % : Beam = #events with BH2 in [lo,hi]
// Overkill(%) = 100 × (BeamVeto / Beam)   // 기존 출력 유지
// Notes:
//   - HTOF tile ID = TParticle::StatusCode()  (0..33)
//   - BH2 seg from world (x) via E72-like mapping
//   - VALID hit: sum(edep) ≥ threshold; energy = *_edep if exists else Weight()
//   - Optional analysis-side exclusion: excludeTilesCSV (e.g. "20,21") → full-exclude(멀티/비토 공통)
//
// Usage:
//   root -l
//   .L Trigger_study_minusreaction_10_28_HTOFselect.C+
//   Trigger_study_minusreaction_10_28_HTOFselect("../rootfile/E45_Nov_beamminus_105.root","g4hyptpc",
//                                                4,10, 0.10,2.0, 5.0,10.0, "");
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

namespace TSMRsel {

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

// parse "a,b,c" → {a,b,c}  (full-exclude)
static std::set<int> ParseCSVInt(const char* csv){
  std::set<int> s; if(!csv) return s;
  std::stringstream ss(csv); std::string tok;
  while(std::getline(ss, tok, ',')){ if(tok.empty()) continue; s.insert(std::stoi(tok)); }
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
  // base tallies
  long long N_total=0;
  long long N_beam =0;      // BH2 in [lo,hi]
  long long N_trig1_base=0; // Beam && (MP_base>=2)

  // veto fired (Beam && MP_base>=2에서 평가)
  long long veto_tight=0, veto_fit=0, veto_wide=0, veto_ultra=0;

  // "배제 때문에 잃는" 이벤트
  long long Lost_MP=0;                 // Beam && MP_base>=2 && MP_excl<2
  long long Lost_Trig2_tight=0;        // Beam && Trig2_base_tight && MP_excl<2
  long long Lost_Trig2_fit=0;          // Beam && Trig2_base_fit   && MP_excl<2
  long long Lost_Trig2_wide=0;         // Beam && Trig2_base_wide  && MP_excl<2
  long long Lost_Trig2_ultra=0;        // Beam && Trig2_base_ultra && MP_excl<2

  // diag: pair presence among Beam && MP_base>=2
  long long have19_20=0, have20_21=0, have21_22=0, have22_23=0, have23_24=0;
};

static Counts ProcessRange_HTOFselect(TTree* T,
                           int bh2_lo, int bh2_hi,
                           double thrBH2, double thrHTOF,
                           bool hasBH2edep, bool hasHTOFedep,
                           const std::set<int>& fullExclude,   // multiplicity & veto 공통 제외
                           const std::set<int>& mpOnlyExclude) // MP 계산에서만 제외 ({1,2,3,4})
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

    bool beam=false; for(int s:bh2Valid){ if(bh2_lo<=s && s<=bh2_hi){ beam=true; break; } }
    if(beam) C.N_beam++;

    // ---- HTOF: 두 경로로 집계
    //   tiles_full : fullExclude만 적용 → base veto/페어/MP_base용
    //   tiles_mp   : fullExclude + mpOnlyExclude 적용 → MP_excl용
    std::map<int,double> htofE_full, htofE_mp;
    if(HTOF){
      for(size_t i=0;i<HTOF->size();++i){
        const TParticle& p = HTOF->at(i);
        int tid = p.GetStatusCode();
        if(!(0<=tid && tid<kNHTOF)) continue;

        if(fullExclude.count(tid)) continue; // 완전 제외

        const double ed = (hasHTOFedep && HTOF_edep && i<HTOF_edep->size())
                          ? HTOF_edep->at(i) : p.GetWeight();

        htofE_full[tid]+=ed;                         // 항상 누적
        if(!mpOnlyExclude.count(tid)) htofE_mp[tid]+=ed; // MP-only 제외 반영
      }
    }
    std::set<int> tiles_full, tiles_mp;
    for(const auto& kv:htofE_full) if(kv.second>=thrHTOF) tiles_full.insert(kv.first);
    for(const auto& kv:htofE_mp)   if(kv.second>=thrHTOF) tiles_mp.insert(kv.first);

    const int MP_base = (int)tiles_full.size();
    const int MP_excl = (int)tiles_mp.size();

    // ---- base Trig1 / Veto / Trig2
    if(beam && MP_base>=2){
      C.N_trig1_base++;

      // 진단: 페어 존재
      const bool p19_20 = PairHit(tiles_full,19,20);
      const bool p20_21 = PairHit(tiles_full,20,21);
      const bool p21_22 = PairHit(tiles_full,21,22);
      const bool p22_23 = PairHit(tiles_full,22,23);
      const bool p23_24 = PairHit(tiles_full,23,24);
      if(p19_20) C.have19_20++;
      if(p20_21) C.have20_21++;
      if(p21_22) C.have21_22++;
      if(p22_23) C.have22_23++;
      if(p23_24) C.have23_24++;

      const VetoFlags v = EvalVetoMinus_ORpairs(tiles_full);
      if(v.tight) C.veto_tight++;
      if(v.fit)   C.veto_fit++;
      if(v.wide)  C.veto_wide++;
      if(v.ultra) C.veto_ultra++;

      // (A) Lost_MP: base에서는 MP>=2였는데, mp-only 제외 후 MP<2로 떨어짐
      if(MP_excl < 2) C.Lost_MP++;

      // (B) Lost_Trig2_*: base에서 각 variant Trig2 생존했으나, MP_excl<2라서 생존 불가
      const bool Trig2_base_tight = (!v.tight);
      const bool Trig2_base_fit   = (!v.fit);
      const bool Trig2_base_wide  = (!v.wide);
      const bool Trig2_base_ultra = (!v.ultra);

      if(MP_excl < 2){
        if(Trig2_base_tight) ++C.Lost_Trig2_tight;
        if(Trig2_base_fit)   ++C.Lost_Trig2_fit;
        if(Trig2_base_wide)  ++C.Lost_Trig2_wide;
        if(Trig2_base_ultra) ++C.Lost_Trig2_ultra;
      }
    }
  }
  return C;
}

struct PrintCfg {
  long long baseBeamSec1=0; // 섹션1(BH2 4–10)의 Beam(스케일 기준)
};

static void PrintSection(const char* title, const Counts& C, int lo, int hi, const PrintCfg& cfg){
  auto pbeam = [&](long long x){ return pct(x, C.N_beam); };
  auto scale_to_1M_Sec1 = [&](long long x)->long long{
    if(cfg.baseBeamSec1<=0) return 0;
    return (long long) llround( 1e6 * (double)x / (double)cfg.baseBeamSec1 );
  };
  auto scale_to_1M_THIS = [&](long long x)->long long{
    if(C.N_beam<=0) return 0;
    return (long long) llround( 1e6 * (double)x / (double)C.N_beam );
  };

  std::cout<<"\n== "<<title<<" | BH2 "<<lo<<"–"<<hi<<" ==\n";
  std::cout<<"Total events                 : "<<C.N_total<<"\n";
  std::cout<<"Beam  (BH2 in-range)         : "<<C.N_beam<<"\n";
  std::cout<<"Trig1_base (Beam && MP_base>=2) : "<<C.N_trig1_base<<"\n";

  std::cout<<std::fixed<<std::setprecision(3);
  std::cout<<"[BeamVeto fired | Beam-based %]\n";
  std::cout<<"  tight      : "<<C.veto_tight<<"  ("<<pbeam(C.veto_tight)<<" % of Beam)\n";
  std::cout<<"  fit        : "<<C.veto_fit  <<"  ("<<pbeam(C.veto_fit)  <<" % of Beam)\n";
  std::cout<<"  wide       : "<<C.veto_wide <<"  ("<<pbeam(C.veto_wide) <<" % of Beam)\n";
  std::cout<<"  ultra      : "<<C.veto_ultra<<"  ("<<pbeam(C.veto_ultra)<<" % of Beam)\n";

  std::cout<<"[Overkill (%) = 100 × BeamVeto / Beam]  // base 정의 유지\n";
  std::cout<<"  tight      : "<<pbeam(C.veto_tight)<<"\n";
  std::cout<<"  fit        : "<<pbeam(C.veto_fit)  <<"\n";
  std::cout<<"  wide       : "<<pbeam(C.veto_wide) <<"\n";
  std::cout<<"  ultra      : "<<pbeam(C.veto_ultra)<<"\n";

  // --- NEW: 배제 때문에 잃는 이벤트 ---
  std::cout<<"[LOSS due to MP-only exclusion of HTOF {1,2,3,4}]  (Beam-based %, + scaled counts)\n";
  std::cout<<"  Lost_MP (Beam && MP_base>=2 → MP_excl<2)\n";
  std::cout<<"    count="<<C.Lost_MP<<", frac="<<pbeam(C.Lost_MP)<<" %"
           <<", scaled@Beam(THIS=1e6)="<<scale_to_1M_THIS(C.Lost_MP)
           <<", scaled@Beam(4-10=1e6)="<<scale_to_1M_Sec1(C.Lost_MP)<<"\n";

  std::cout<<"  Lost_Trig2_tight (base Trig2_tight survived → MP_excl<2 → drop)\n";
  std::cout<<"    count="<<C.Lost_Trig2_tight<<", frac="<<pbeam(C.Lost_Trig2_tight)<<" %"
           <<", scaled@THIS="<<scale_to_1M_THIS(C.Lost_Trig2_tight)
           <<", scaled@Sec1="<<scale_to_1M_Sec1(C.Lost_Trig2_tight)<<"\n";

  std::cout<<"  Lost_Trig2_fit   (base Trig2_fit   survived → MP_excl<2 → drop)\n";
  std::cout<<"    count="<<C.Lost_Trig2_fit<<", frac="<<pbeam(C.Lost_Trig2_fit)<<" %"
           <<", scaled@THIS="<<scale_to_1M_THIS(C.Lost_Trig2_fit)
           <<", scaled@Sec1="<<scale_to_1M_Sec1(C.Lost_Trig2_fit)<<"\n";

  std::cout<<"  Lost_Trig2_wide  (base Trig2_wide  survived → MP_excl<2 → drop)\n";
  std::cout<<"    count="<<C.Lost_Trig2_wide<<", frac="<<pbeam(C.Lost_Trig2_wide)<<" %"
           <<", scaled@THIS="<<scale_to_1M_THIS(C.Lost_Trig2_wide)
           <<", scaled@Sec1="<<scale_to_1M_Sec1(C.Lost_Trig2_wide)<<"\n";

  std::cout<<"  Lost_Trig2_ultra (base Trig2_ultra survived → MP_excl<2 → drop)\n";
  std::cout<<"    count="<<C.Lost_Trig2_ultra<<", frac="<<pbeam(C.Lost_Trig2_ultra)<<" %"
           <<", scaled@THIS="<<scale_to_1M_THIS(C.Lost_Trig2_ultra)
           <<", scaled@Sec1="<<scale_to_1M_Sec1(C.Lost_Trig2_ultra)<<"\n";

  // 진단: 페어 존재 빈도
  std::cout<<"[Diag] Pair presence within (Beam && MP_base>=2):\n";
  std::cout<<"  (19,20)="<<C.have19_20<<", (20,21)="<<C.have20_21
           <<", (21,22)="<<C.have21_22<<", (22,23)="<<C.have22_23
           <<", (23,24)="<<C.have23_24<<"\n";
}

} // namespace TSMRsel


void Trigger_study_minusreaction_10_28_HTOFselect(const char* filename,
                                       const char* treename="g4hyptpc",
                                       int bh2_lo=4, int bh2_hi=10,
                                       double mipFrac=0.10,
                                       double mipMeVperCm=2.0,
                                       double BH2_thickness_mm=5.0,
                                       double HTOF_thickness_mm=10.0,
                                       const char* excludeTilesCSV="" ) // e.g. "20,21"
{
  using namespace TSMRsel;

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

  // thresholds & excludes
  const double thrBH2  = MIPThrMeV(mipFrac,mipMeVperCm,BH2_thickness_mm);
  const double thrHTOF = MIPThrMeV(mipFrac,mipMeVperCm,HTOF_thickness_mm);
  auto fullExclude = ParseCSVInt(excludeTilesCSV);
  const std::set<int> mpOnlyExclude = {1,2,3,4}; // 요청: MP계산에서만 1–4 제외

  std::cout<<std::fixed<<std::setprecision(3)
           <<"[INFO] BH2 thr="<<thrBH2<<" MeV, HTOF thr="<<thrHTOF<<" MeV\n";
  if(!fullExclude.empty()){
    std::cout<<"[INFO] FULL-EXCLUDE (mult & veto): ";
    bool first=true; for(int t:fullExclude){ if(!first) std::cout<<","; first=false; std::cout<<t; }
    std::cout<<"\n";
  }
  std::cout<<"[INFO] MP-only EXCLUDE for multiplicity: 1,2,3,4 (veto는 full-valid로 평가)\n";

  // sections
  const int s1_lo = bh2_lo, s1_hi = bh2_hi; // Full (4–10)
  const int s2_lo = 4,      s2_hi = 9;      // Narrow 1
  const int s3_lo = 5,      s3_hi = 10;     // Narrow 2

  auto C1 = ProcessRange_HTOFselect(T, s1_lo, s1_hi, thrBH2, thrHTOF, hasBH2edep, hasHTOFedep,
                                    fullExclude, mpOnlyExclude);
  auto C2 = ProcessRange_HTOFselect(T, s2_lo, s2_hi, thrBH2, thrHTOF, hasBH2edep, hasHTOFedep,
                                    fullExclude, mpOnlyExclude);
  auto C3 = ProcessRange_HTOFselect(T, s3_lo, s3_hi, thrBH2, thrHTOF, hasBH2edep, hasHTOFedep,
                                    fullExclude, mpOnlyExclude);

  // 스케일 기준(Beam 4–10 = 1,000,000)
  TSMRsel::PrintCfg cfg;
  cfg.baseBeamSec1 = C1.N_beam;

  TSMRsel::PrintSection("Section 1 (Full Beam, MP-only excl {1,2,3,4})", C1, s1_lo, s1_hi, cfg);
  TSMRsel::PrintSection("Section 2 (Narrow 1, MP-only excl {1,2,3,4})",  C2, s2_lo, s2_hi, cfg);
  TSMRsel::PrintSection("Section 3 (Narrow 2, MP-only excl {1,2,3,4})",  C3, s3_lo, s3_hi, cfg);
}

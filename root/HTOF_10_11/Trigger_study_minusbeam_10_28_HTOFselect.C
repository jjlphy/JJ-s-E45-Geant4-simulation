// -*- C++ -*-
// Trigger_study_minusbeam_10_28_HTOFselect.C  (2025-10-28 for jaejin)
// Minus-beam version with OR-based BeamVeto over adjacent HTOF pairs.
//
// NEW in this revision:
//   - Added Section0 (baseline, no MP-only exclusion) and reports:
//       * N( Beam && MP>=2 && hits any of {1..4} )
//       * N( Beam && MP>=2 && hits any of {0..5} )
//     with absolute counts and Beam-based %.
//   - Sections 1–3 : MP-only exclusion {1,2,3,4} (veto uses full set)
//   - Sections 4–6 : MP-only exclusion {0,1,2,3,4,5} (veto uses full set)
//
// Sections (now 7 in total):
//   Sec0: BH2 4–10 (Baseline, no MP-only exclusion) → report MP≥2∧hit(1–4) and MP≥2∧hit(0–5)
//   Sec1: BH2 4–10 (Full Beam)   , MP-excl {1,2,3,4}
//   Sec2: BH2 4–9  (Narrow 1)    , MP-excl {1,2,3,4}
//   Sec3: BH2 5–10 (Narrow 2)    , MP-excl {1,2,3,4}
//   Sec4: BH2 4–10 (Full Beam)   , MP-excl {0,1,2,3,4,5}
//   Sec5: BH2 4–9  (Narrow 1)    , MP-excl {0,1,2,3,4,5}
//   Sec6: BH2 5–10 (Narrow 2)    , MP-excl {0,1,2,3,4,5}
//
// Trigger logic (각 섹션 공통):
//   Beam   = (BH2 in [lo,hi])
//   Trig1  = Beam && (HTOF multiplicity >= 2)   // using MP-only exclusion set of the section
//   Trig2* = Trig1 && !BeamVeto(*)
//
// BeamVeto (인접 페어 = 두 타일 동시 히트; 콤마는 AND):
//   pair(a,b) := (tile a) && (tile b)
//   tight      : pair(20,21) OR pair(21,22)
//   fit        : tight       OR pair(22,23)
//   wide       : fit         OR pair(19,20)
//   ultra-wide : wide        OR pair(23,24)
//   ※ 비토는 "전체 유효 타일"로 판정 (MP-only 제외는 적용하지 않음)
//
// Denominator for % : Beam(= BH2 in-range) 이벤트 개수
// 출력(각 섹션):
//  - 절대 개수 + Beam 기준 %
//  - Beam Veto efficiency = 100% - (Trig2/Beam)*100   (소수점 3자리)
//  - Beam-induced background (counts) @ Beam=1,000,000
//  - Section0 추가 리포트: MP≥2∧hit(1–4), MP≥2∧hit(0–5) (절대/퍼센트)
//
// 옵션 excludeTilesCSV: 완전 제외(멀티/비토 둘 다에서) 디버그용 추가 마스크.
//   예: "20,21" → 해당 타일은 멀티·비토 모두에서 무시.
//   본 요청의 MP-only 제외는 여기에 더해 따로 적용됨.
//
// 사용:
//   root -l
//   .L Trigger_study_minusbeam_10_28_HTOFselect.C+
//   Trigger_study_minusbeam_10_28_HTOFselect("../rootfile/E45_Nov_beamminus_098.root","g4hyptpc",
//                                            4,10, 0.10,2.0, 5.0,10.0, "");
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

namespace TSMsel {

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

// parse "a,b,c" → {a,b,c} (완전제외용)
static std::set<int> ParseCSVInt(const char* csv){
  std::set<int> s; if(!csv) return s;
  std::stringstream ss(csv); std::string tok;
  while(std::getline(ss, tok, ',')){ if(tok.empty()) continue; s.insert(std::stoi(tok)); }
  return s;
}

// ===== OR-logic BeamVeto over adjacent pairs (minus-beam, centered 19–24) =====
static inline bool PairHit(const std::set<int>& tiles, int a, int b){
  return tiles.count(a) && tiles.count(b);
}
struct VetoSet { bool tight=false, fit=false, wide=false, ultra=false; };
static inline VetoSet EvalVeto_ORpairs(const std::set<int>& tiles){
  VetoSet v;
  const bool p19_20 = PairHit(tiles,19,20);
  const bool p20_21 = PairHit(tiles,20,21);
  const bool p21_22 = PairHit(tiles,21,22);
  const bool p22_23 = PairHit(tiles,22,23);
  const bool p23_24 = PairHit(tiles,23,24);

  v.tight = (p20_21 || p21_22);
  v.fit   = (v.tight || p22_23);
  v.wide  = (v.fit   || p19_20);
  v.ultra = (v.wide  || p23_24);
  return v;
}

struct Counts {
  long long N_total=0;
  long long N_beam =0;     // BH2 in [lo,hi]
  long long N_trig1=0;     // Beam && MP>=2 (using MP-only exclusion)
  // survivors (!veto)
  long long N_trig2_tight=0, N_trig2_fit=0, N_trig2_wide=0, N_trig2_ultra=0;
  // veto fired counts
  long long N_veto_tight=0, N_veto_fit=0, N_veto_wide=0, N_veto_ultra=0;
  // tile presence (diagnostic on full valid set)
  long long have19=0, have20=0, have21=0, have22=0, have23=0, have24=0;

  // Section0 전용 리포트: MP≥2 ∧ hit(1–4) / hit(0–5)
  long long N_trig1_with_1_4 = 0;
  long long N_trig1_with_0_5 = 0;
};

static Counts ProcessRange_HTOFselect(TTree* T,
                           int bh2_lo, int bh2_hi,
                           double thrBH2, double thrHTOF,
                           bool hasBH2edep, bool hasHTOFedep,
                           const std::set<int>& fullExclude,        // 완전제외(멀티/비토 공통)
                           const std::set<int>& mpOnlyExclude)      // MP-only 제외(비토에는 미적용)
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

    // ---- HTOF valid accumulation (two paths):
    // path A: "fullValid" = 완전제외(fullExclude)만 적용 → Veto/집합표시용
    // path B: "mpValid"   = 완전제외 + MP-only 제외(mpOnlyExclude) → multiplicity 판정용
    std::map<int,double> htofE_full, htofE_mp;
    if(HTOF){
      for(size_t i=0;i<HTOF->size();++i){
        const TParticle& p=HTOF->at(i);
        int tid = p.GetStatusCode();
        if(!(0<=tid && tid<kNHTOF)) continue;

        // 완전 제외 마스크: 둘 다에서 제외
        if(fullExclude.count(tid)) continue;

        const double ed=(hasHTOFedep && HTOF_edep && i<HTOF_edep->size())? HTOF_edep->at(i) : p.GetWeight();
        // fullValid에 먼저 누적
        htofE_full[tid]+=ed;
        // mpValid는 mpOnlyExclude에 걸리면 배제
        if(!mpOnlyExclude.count(tid)) htofE_mp[tid]+=ed;
      }
    }

    std::set<int> tiles_full, tiles_mp;
    for(const auto& kv:htofE_full) if(kv.second>=thrHTOF) tiles_full.insert(kv.first);
    for(const auto& kv:htofE_mp)   if(kv.second>=thrHTOF) tiles_mp.insert(kv.first);

    const int mult = (int)tiles_mp.size();  // MP>=2는 mp-only exclusion이 반영된 집합으로 평가

    // presence (diagnostic; fullValid 기준)
    if(tiles_full.count(19)) C.have19++;
    if(tiles_full.count(20)) C.have20++;
    if(tiles_full.count(21)) C.have21++;
    if(tiles_full.count(22)) C.have22++;
    if(tiles_full.count(23)) C.have23++;
    if(tiles_full.count(24)) C.have24++;

    // ---- Trig1 & Trig2(*)  (Veto는 tiles_full로 평가!)
    if(beam && mult>=2){
      C.N_trig1++;
      const auto v = EvalVeto_ORpairs(tiles_full);
      if(!v.tight) C.N_trig2_tight++; else C.N_veto_tight++;
      if(!v.fit)   C.N_trig2_fit++;   else C.N_veto_fit++;
      if(!v.wide)  C.N_trig2_wide++;  else C.N_veto_wide++;
      if(!v.ultra) C.N_trig2_ultra++; else C.N_veto_ultra++;

      // ----- Section0 전용 집계: MP≥2 ∧ hit(1–4), MP≥2 ∧ hit(0–5) -----
      bool hit_1_4 = false, hit_0_5 = false;
      for(int t=1; t<=4; ++t) if(tiles_full.count(t)) { hit_1_4=true; break; }
      for(int t=0; t<=5; ++t) if(tiles_full.count(t)) { hit_0_5=true; break; }
      if(hit_1_4) C.N_trig1_with_1_4++;
      if(hit_0_5) C.N_trig1_with_0_5++;
    }
  }
  return C;
}

static void PrintSection(const char* title, const Counts& C, int lo, int hi,
                         const std::set<int>& mpOnlyExclude)
{
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

  // MP-only exclusion 표시
  if(!mpOnlyExclude.empty()){
    std::cout<<"[MP-only exclusion for multiplicity] ";
    bool first=true; for(int t:mpOnlyExclude){ if(!first) std::cout<<","; first=false; std::cout<<t; }
    std::cout<<"\n";
  }

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

  std::cout<<"[Tile presence among Beam events (full-valid, for diagnosis)]\n";
  std::cout<<"  have19="<<C.have19<<", have20="<<C.have20<<", have21="<<C.have21
           <<", have22="<<C.have22<<", have23="<<C.have23<<", have24="<<C.have24<<"\n";
}

// ---- Section0 전용 출력 (Baseline; MP-only exclusion 없음) ----
static void PrintSection0(const Counts& C0, int lo, int hi){
  auto p = [&](long long x)->double{ return pct(x, C0.N_beam); };

  std::cout<<"\n== Section 0 (Baseline; no MP-only exclusion) | BH2 "<<lo<<"–"<<hi<<" ==\n";
  std::cout<<"Total events                 : "<<C0.N_total<<"\n";
  std::cout<<"Beam  (BH2 in-range)         : "<<C0.N_beam<<"\n";
  std::cout<<std::fixed<<std::setprecision(3);
  std::cout<<"[MP>=2 while passing specific HTOF groups]\n";
  std::cout<<"  MP>=2 & hit(1–4)           : "<<C0.N_trig1_with_1_4<<"  ("<<p(C0.N_trig1_with_1_4)<<" % of Beam)\n";
  std::cout<<"  MP>=2 & hit(0–5)           : "<<C0.N_trig1_with_0_5<<"  ("<<p(C0.N_trig1_with_0_5)<<" % of Beam)\n";
}

} // namespace TSMsel


void Trigger_study_minusbeam_10_28_HTOFselect(const char* filename,
                                   const char* treename="g4hyptpc",
                                   int bh2_lo=4, int bh2_hi=10,
                                   double mipFrac=0.10,
                                   double mipMeVperCm=2.0,
                                   double BH2_thickness_mm=5.0,
                                   double HTOF_thickness_mm=10.0,
                                   const char* excludeTilesCSV="" ) // 완전제외(멀티/비토 공통)
{
  using namespace TSMsel;

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

  // thresholds & masks
  const double thrBH2  = MIPThrMeV(mipFrac,mipMeVperCm,BH2_thickness_mm);
  const double thrHTOF = MIPThrMeV(mipFrac,mipMeVperCm,HTOF_thickness_mm);
  auto fullExclude = ParseCSVInt(excludeTilesCSV);

  std::cout<<std::fixed<<std::setprecision(3)
           <<"[INFO] BH2 thr="<<thrBH2<<" MeV, HTOF thr="<<thrHTOF<<" MeV\n";
  if(!fullExclude.empty()){
    std::cout<<"[INFO] FULL EXCLUDE (both mult & veto): ";
    bool first=true; for(int t:fullExclude){ if(!first) std::cout<<","; first=false; std::cout<<t; }
    std::cout<<"\n";
  }

  // Sections BH2 ranges
  const int S1_lo = bh2_lo, S1_hi = bh2_hi; // Full Beam
  const int S2_lo = 4,      S2_hi = 9;      // Narrow 1
  const int S3_lo = 5,      S3_hi = 10;     // Narrow 2

  // MP-only exclusion sets
  const std::set<int> MP_EXCL_1_3 = {1,2,3,4};
  const std::set<int> MP_EXCL_4_6 = {0,1,2,3,4,5};

  // ---- Section 0 (Baseline; no MP-only exclusion) ----
  const std::set<int> MP_EXCL_0 = {}; // none
  auto C0 = ProcessRange_HTOFselect(T, S1_lo, S1_hi, thrBH2, thrHTOF,
                                    hasBH2edep, hasHTOFedep, fullExclude, MP_EXCL_0);
  PrintSection0(C0, S1_lo, S1_hi);

  // ---- Sections 1–3 (MP-only exclude {1,2,3,4}) ----
  auto C1 = ProcessRange_HTOFselect(T, S1_lo, S1_hi, thrBH2, thrHTOF, hasBH2edep, hasHTOFedep,
                                    fullExclude, MP_EXCL_1_3);
  auto C2 = ProcessRange_HTOFselect(T, S2_lo, S2_hi, thrBH2, thrHTOF, hasBH2edep, hasHTOFedep,
                                    fullExclude, MP_EXCL_1_3);
  auto C3 = ProcessRange_HTOFselect(T, S3_lo, S3_hi, thrBH2, thrHTOF, hasBH2edep, hasHTOFedep,
                                    fullExclude, MP_EXCL_1_3);

  // ---- Sections 4–6 (MP-only exclude {0,1,2,3,4,5}) ----
  auto C4 = ProcessRange_HTOFselect(T, S1_lo, S1_hi, thrBH2, thrHTOF, hasBH2edep, hasHTOFedep,
                                    fullExclude, MP_EXCL_4_6);
  auto C5 = ProcessRange_HTOFselect(T, S2_lo, S2_hi, thrBH2, thrHTOF, hasBH2edep, hasHTOFedep,
                                    fullExclude, MP_EXCL_4_6);
  auto C6 = ProcessRange_HTOFselect(T, S3_lo, S3_hi, thrBH2, thrHTOF, hasBH2edep, hasHTOFedep,
                                    fullExclude, MP_EXCL_4_6);

  // Print 1–6
  PrintSection("Section 1 (Full Beam) MP-excl{1,2,3,4}", C1, S1_lo, S1_hi, MP_EXCL_1_3);
  PrintSection("Section 2 (Narrow 1) MP-excl{1,2,3,4}",  C2, S2_lo, S2_hi, MP_EXCL_1_3);
  PrintSection("Section 3 (Narrow 2) MP-excl{1,2,3,4}",  C3, S3_lo, S3_hi, MP_EXCL_1_3);

  PrintSection("Section 4 (Full Beam) MP-excl{0,1,2,3,4,5}", C4, S1_lo, S1_hi, MP_EXCL_4_6);
  PrintSection("Section 5 (Narrow 1) MP-excl{0,1,2,3,4,5}",  C5, S2_lo, S2_hi, MP_EXCL_4_6);
  PrintSection("Section 6 (Narrow 2) MP-excl{0,1,2,3,4,5}",  C6, S3_lo, S3_hi, MP_EXCL_4_6);
}

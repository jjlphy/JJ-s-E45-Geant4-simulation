// -*- C++ -*-
// Trigger_study_plusbeam_10_28_HTOFselect.C  (2025-10-29)
// - Plus-beam (π+), OR-based BeamVeto over adjacent pairs around 16–21
// - Baseline(Section0), Sections1–6
// - Survived combo listing + ROOT-CAUSE diagnostics (EM/Secondary/Primary-π+ only)
//
// Usage:
//   root -l
//   .L Trigger_study_plusbeam_10_28_HTOFselect.C+
//   Trigger_study_plusbeam_10_28_HTOFselect("E45.root","g4hyptpc",
//       4,10, 0.10,2.0, 5.0,10.0, /*excludeTilesCSV=*/"",
//       /*primaryOriginCode=*/0, /*excludePDGCSV=*/"", /*TOPN=*/20);

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
#include <cctype>

namespace TSPsel {

static const int    kNBH2Seg   = 15;     // 0..14
static const double kBH2_x0    = -10.0;  // [mm]
static const double kBH2_pitch = 14.0;   // [mm]
static const double kBH2_total = kNBH2Seg * kBH2_pitch;
static const int    kNHTOF     = 34;     // 0..33

// PDG codes
static const int PDG_E_MINUS = 11;
static const int PDG_E_PLUS  = -11;
static const int PDG_GAMMA   = 22;
static const int PDG_PI_PLUS = +211;

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

// parse CSV "a,b,c" → set<int>
static std::set<int> ParseCSVInt(const char* csv){
  std::set<int> s; if(!csv) return s;
  std::stringstream ss(csv); std::string tok;
  while(std::getline(ss, tok, ',')){ if(tok.empty()) continue; s.insert(std::stoi(tok)); }
  return s;
}
// parse CSV "11,-11,22" → set<int> (PDG)
static std::set<int> ParseCSVPDG(const char* csv){
  std::set<int> s; if(!csv) return s;
  std::stringstream ss(csv); std::string tok;
  while(std::getline(ss, tok, ',')){ if(tok.empty()) continue; s.insert(std::stoi(tok)); }
  return s;
}

// set<int> → "5,20,23"
static std::string KeyFromSet(const std::set<int>& tiles){
  if(tiles.empty()) return std::string("-");
  std::ostringstream os; bool first=true;
  for(int t: tiles){ if(!first) os<<","; first=false; os<<t; }
  return os.str();
}

// ===== OR-logic BeamVeto over adjacent pairs (16–21 중심, π+) =====
static inline bool PairHit(const std::set<int>& tiles, int a, int b){
  return tiles.count(a) && tiles.count(b);
}
struct VetoSet { bool tight=false, fit=false, wide=false, ultra=false; };
static inline VetoSet EvalVeto_ORpairs(const std::set<int>& tiles){
  VetoSet v;
  const bool p16_17 = PairHit(tiles,16,17);
  const bool p17_18 = PairHit(tiles,17,18);
  const bool p18_19 = PairHit(tiles,18,19);
  const bool p19_20 = PairHit(tiles,19,20);
  const bool p20_21 = PairHit(tiles,20,21);

  v.tight = (p17_18 || p18_19);
  v.fit   = (v.tight || p19_20);
  v.wide  = (v.fit   || p16_17);
  v.ultra = (v.wide  || p20_21);
  return v;
}

// ===== per-tile composition flags (유효 타일 판단 이후에 쓰기) =====
struct TileInfo {
  double sumE = 0.0;           // Σedep (MeV) for the tile
  bool   hasEM = false;        // any e± or γ contributor present?
  bool   hasSecondary = false; // any contributor with origin != primary?
  bool   allPiPrimary = true;  // all contributors are (origin=primary && PDG = π+)?
};

struct Counts {
  long long N_total=0;
  long long N_beam =0;     // BH2 in [lo,hi]
  long long N_trig1=0;     // Beam && MP>=2 (using MP-only exclusion)
  long long N_trig2_tight=0, N_trig2_fit=0, N_trig2_wide=0, N_trig2_ultra=0;
  long long N_veto_tight=0,  N_veto_fit=0,  N_veto_wide=0,  N_veto_ultra=0;

  // diagnostics (full-valid)
  long long have16=0, have17=0, have18=0, have19=0, have20=0, have21=0;

  // Section0-only
  long long N_trig1_with_1_4=0, N_trig1_with_0_5=0;

  // survived combo tallies
  std::map<std::string,long long> Comb_tight, Comb_fit, Comb_wide, Comb_ultra;

  // root-cause tallies (by veto region)
  struct RC {
    long long ALL_PI_PRIMARY=0;   // 모든 유효 타일이 Primary π+ 만
    long long HAS_EM=0;           // e±/γ 기여 포함
    long long HAS_SECONDARY=0;    // 비-Primary 기여 포함
    double minSumE = 1e99;        // 생존 이벤트들의 타일 ΣEdep 중 최소/최대
    double maxSumE = 0.0;
    std::map<std::string,long long> top_allPi, top_hasEM, top_hasSec;
  };
  RC rc_tight, rc_fit, rc_wide, rc_ultra;
};

static void UpdateEnergyExtrema(Counts::RC& rc, const std::map<int,TileInfo>& tilesInfo){
  for(const auto& kv:tilesInfo){
    rc.minSumE = std::min(rc.minSumE, kv.second.sumE);
    rc.maxSumE = std::max(rc.maxSumE, kv.second.sumE);
  }
}

static void BumpRC(Counts::RC& rc,
                   const std::string& key,
                   bool allPiPrimary, bool hasEM, bool hasSec)
{
  if(allPiPrimary) rc.ALL_PI_PRIMARY++;
  if(hasEM)        rc.HAS_EM++;
  if(hasSec)       rc.HAS_SECONDARY++;

  if(allPiPrimary) rc.top_allPi[key]++;
  if(hasEM)        rc.top_hasEM[key]++;
  if(hasSec)       rc.top_hasSec[key]++;
}

// 정렬 출력 (Top-N)
static void PrintPairs(const char* head, const std::map<std::string,long long>& M, int TOPN){
  std::vector<std::pair<std::string,long long>> v(M.begin(), M.end());
  std::sort(v.begin(), v.end(),
            [](auto& a, auto& b){
              if(a.second!=b.second) return a.second>b.second;
              return a.first<b.first;
            });
  std::cout<<head<<"\n";
  if(v.empty()){ std::cout<<"  (none)\n"; return; }
  int n=0; for(const auto& kv: v){ if(n++>=TOPN) break;
    std::cout<<"  comb("<<kv.first<<"): "<<kv.second<<"\n";
  }
}

// 메인 처리
static Counts ProcessRange_HTOFselect(TTree* T,
                           int bh2_lo, int bh2_hi,
                           double thrBH2, double thrHTOF,
                           bool hasBH2edep, bool hasHTOFedep,
                           const std::set<int>& fullExclude,
                           const std::set<int>& mpOnlyExclude,
                           int primaryOriginCode,
                           const std::set<int>& excludePDG)
{
  std::vector<TParticle>* BH2=nullptr;   T->SetBranchAddress("BH2",&BH2);
  std::vector<TParticle>* HTOF=nullptr;  T->SetBranchAddress("HTOF",&HTOF);
  std::vector<double>* BH2_edep=nullptr; if(hasBH2edep)  T->SetBranchAddress("BH2_edep",&BH2_edep);
  std::vector<double>* HTOF_edep=nullptr;if(hasHTOFedep) T->SetBranchAddress("HTOF_edep",&HTOF_edep);

  Counts C;
  const Long64_t N=T->GetEntries();
  for(Long64_t ie=0; ie<N; ++ie){
    T->GetEntry(ie); C.N_total++;

    // ---- BH2 valid ----
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

    // ---- HTOF accumulation (full vs mp) + composition flags ----
    std::map<int,TileInfo> info_full; // per tile → ΣE, flags
    std::map<int,double>   accum_mp;  // multiplicity set만

    if(HTOF){
      for(size_t i=0;i<HTOF->size();++i){
        const TParticle& p=HTOF->at(i);
        const int tid = p.GetStatusCode();
        if(!(0<=tid && tid<kNHTOF)) continue;
        if(fullExclude.count(tid))  continue;

        const int pdg = p.GetPdgCode();
        if(!excludePDG.empty() && excludePDG.count(pdg)) continue;

        const double ed=(hasHTOFedep && HTOF_edep && i<HTOF_edep->size())? HTOF_edep->at(i) : p.GetWeight();

        // --- full-valid 타일 구성 정보 누적 ---
        auto& tf = info_full[tid];
        tf.sumE += ed;

        const bool isEM = (pdg==PDG_GAMMA || pdg==PDG_E_MINUS || pdg==PDG_E_PLUS);
        if(isEM) tf.hasEM = true;

        const unsigned int uid = static_cast<unsigned int>(p.GetUniqueID());
        const unsigned int pri = static_cast<unsigned int>(primaryOriginCode);
        if(uid!=pri) tf.hasSecondary = true;

        const bool isPrimaryPi = (uid==pri && pdg==PDG_PI_PLUS);
        if(!isPrimaryPi) tf.allPiPrimary = false;

        // --- mp-valid (MP-only exclusion 적용) ---
        if(!mpOnlyExclude.count(tid)) accum_mp[tid] += ed;
      }
    }

    // 타일 유효 판정
    std::set<int> tiles_full, tiles_mp;
    for(const auto& kv:info_full) if(kv.second.sumE>=thrHTOF) tiles_full.insert(kv.first);
    for(const auto& kv:accum_mp)  if(kv.second      >=thrHTOF) tiles_mp.insert(kv.first);

    // presence (diagnostic) — π+ 중심영역(16–21)
    if(tiles_full.count(16)) C.have16++;
    if(tiles_full.count(17)) C.have17++;
    if(tiles_full.count(18)) C.have18++;
    if(tiles_full.count(19)) C.have19++;
    if(tiles_full.count(20)) C.have20++;
    if(tiles_full.count(21)) C.have21++;

    // ---- Triggers ----
    if(beam && (int)tiles_mp.size()>=2){
      C.N_trig1++;
      const auto v = EvalVeto_ORpairs(tiles_full);
      const std::string key = KeyFromSet(tiles_full);

      bool anyEM=false, anySecondary=false, allPiPrimary=true;
      for(int t: tiles_full){
        const auto& tf = info_full[t];
        anyEM        = anyEM        || tf.hasEM;
        anySecondary = anySecondary || tf.hasSecondary;
        allPiPrimary = allPiPrimary && tf.allPiPrimary;
      }

      auto bump_all = [&](Counts::RC& rc){
        UpdateEnergyExtrema(rc, info_full);
        BumpRC(rc, key, allPiPrimary, anyEM, anySecondary);
      };

      if(!v.tight){ C.N_trig2_tight++; C.Comb_tight[key]++;  bump_all(C.rc_tight); }  else { C.N_veto_tight++; }
      if(!v.fit)  { C.N_trig2_fit++;   C.Comb_fit[key]++;    bump_all(C.rc_fit);   }  else { C.N_veto_fit++;   }
      if(!v.wide) { C.N_trig2_wide++;  C.Comb_wide[key]++;   bump_all(C.rc_wide);  }  else { C.N_veto_wide++;  }
      if(!v.ultra){ C.N_trig2_ultra++; C.Comb_ultra[key]++;  bump_all(C.rc_ultra); }  else { C.N_veto_ultra++; }

      // Section0용 (baseline 통계 확인용)
      bool hit_1_4=false, hit_0_5=false;
      for(int t=1; t<=4; ++t) if(tiles_full.count(t)) { hit_1_4=true; break; }
      for(int t=0; t<=5; ++t) if(tiles_full.count(t)) { hit_0_5=true; break; }
      if(hit_1_4) C.N_trig1_with_1_4++;
      if(hit_0_5) C.N_trig1_with_0_5++;
    }
  }
  return C;
}

// 공통 섹션 프린트
static void PrintSection(const char* title, const Counts& C, int lo, int hi,
                         const std::set<int>& mpOnlyExclude, int TOPN)
{
  const double BeamNorm = 1e6;
  auto bg = [&](long long trig2)->long long{
    if(C.N_beam<=0) return 0;
    return (long long) llround( BeamNorm * (double)trig2 / (double)C.N_beam );
  };
  auto p   = [&](long long x)->double{ return pct(x, C.N_beam); };
  auto eff = [&](long long trig2)->double{ return 100.0 - p(trig2); };

  std::cout<<"\n== "<<title<<" | BH2 "<<lo<<"–"<<hi<<" ==\n";
  std::cout<<"Total events                 : "<<C.N_total<<"\n";
  std::cout<<"Beam  (BH2 in-range)         : "<<C.N_beam<<"\n";
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

  std::cout<<"[Tile presence among Beam events (full-valid, 16–21 diag)]\n";
  std::cout<<"  have16="<<C.have16<<", have17="<<C.have17<<", have18="<<C.have18
           <<", have19="<<C.have19<<", have20="<<C.have20<<", have21="<<C.have21<<"\n";

  // ---- Survived combo lists (Top-N) ----
  PrintPairs("== Survived combinations (tight) ==", C.Comb_tight, TOPN);
  PrintPairs("== Survived combinations (fit) ==",   C.Comb_fit,   TOPN);
  PrintPairs("== Survived combinations (wide) ==",  C.Comb_wide,  TOPN);
  PrintPairs("== Survived combinations (ultra) ==", C.Comb_ultra, TOPN);

  // ---- Root-cause summary ----
  auto pr_rc = [&](const char* head, const Counts::RC& rc){
    std::cout<<head<<"\n";
    std::cout<<"  ALL_PI_PRIMARY : "<<rc.ALL_PI_PRIMARY<<"\n";
    std::cout<<"  HAS_EM         : "<<rc.HAS_EM<<"\n";
    std::cout<<"  HAS_SECONDARY  : "<<rc.HAS_SECONDARY<<"\n";
    if(rc.minSumE<=rc.maxSumE){
      std::cout<<std::setprecision(3)
               <<"  ΣEdep over tiles  : min="<<rc.minSumE<<" MeV, max="<<rc.maxSumE<<" MeV\n";
    }
    PrintPairs("    Top combos (ALL_PI_PRIMARY)", rc.top_allPi, 10);
    PrintPairs("    Top combos (HAS_EM)",         rc.top_hasEM, 10);
    PrintPairs("    Top combos (HAS_SECONDARY)",  rc.top_hasSec,10);
  };
  pr_rc("[Root-cause tallies | tight]", C.rc_tight);
  pr_rc("[Root-cause tallies | fit  ]", C.rc_fit);
  pr_rc("[Root-cause tallies | wide ]", C.rc_wide);
  pr_rc("[Root-cause tallies | ultra]", C.rc_ultra);
}

// Section0 전용
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

} // namespace TSPsel


void Trigger_study_plusbeam_10_28_HTOFselect(const char* filename,
                                   const char* treename="g4hyptpc",
                                   int bh2_lo=4, int bh2_hi=10,
                                   double mipFrac=0.10,
                                   double mipMeVperCm=2.0,
                                   double BH2_thickness_mm=5.0,
                                   double HTOF_thickness_mm=10.0,
                                   const char* excludeTilesCSV="" ,
                                   int primaryOriginCode=0,
                                   const char* excludePDGCSV="" ,
                                   int TOPN=20)
{
  using namespace TSPsel;

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
  auto excludePDG  = ParseCSVPDG(excludePDGCSV);

  std::cout<<std::fixed<<std::setprecision(3)
           <<"[INFO] BH2 thr="<<thrBH2<<" MeV, HTOF thr="<<thrHTOF<<" MeV\n";
  if(!fullExclude.empty()){
    std::cout<<"[INFO] FULL EXCLUDE tiles: ";
    bool first=true; for(int t:fullExclude){ if(!first) std::cout<<","; first=false; std::cout<<t; }
    std::cout<<"\n";
  }
  if(!excludePDG.empty()){
    std::cout<<"[INFO] EXCLUDE PDG: ";
    bool first=true; for(int p:excludePDG){ if(!first) std::cout<<","; first=false; std::cout<<p; }
    std::cout<<"\n";
  }
  std::cout<<"[INFO] primaryOriginCode = "<<primaryOriginCode<<"\n";

  // Sections BH2 ranges
  const int S1_lo = bh2_lo, S1_hi = bh2_hi; // Full Beam
  const int S2_lo = 4,      S2_hi = 9;      // Narrow 1
  const int S3_lo = 5,      S3_hi = 10;     // Narrow 2

  // MP-only exclusion sets (멀티≥2 판정에만 적용)
  const std::set<int> MP_EXCL_1_4 = {1,2,3,4};
  const std::set<int> MP_EXCL_0_5 = {0,1,2,3,4,5};
  const std::set<int> MP_EXCL_0   = {}; // none

  // ---- Section 0 (Baseline) ----
  auto C0 = ProcessRange_HTOFselect(T, S1_lo, S1_hi, thrBH2, thrHTOF,
                                    hasBH2edep, hasHTOFedep, fullExclude, MP_EXCL_0,
                                    primaryOriginCode, excludePDG);
  PrintSection0(C0, S1_lo, S1_hi);

  // ---- Sections 1–3 (MP-excl {1,2,3,4}) ----
  auto C1 = ProcessRange_HTOFselect(T, S1_lo, S1_hi, thrBH2, thrHTOF,
                                    hasBH2edep, hasHTOFedep, fullExclude, MP_EXCL_1_4,
                                    primaryOriginCode, excludePDG);
  auto C2 = ProcessRange_HTOFselect(T, S2_lo, S2_hi, thrBH2, thrHTOF,
                                    hasBH2edep, hasHTOFedep, fullExclude, MP_EXCL_1_4,
                                    primaryOriginCode, excludePDG);
  auto C3 = ProcessRange_HTOFselect(T, S3_lo, S3_hi, thrBH2, thrHTOF,
                                    hasBH2edep, hasHTOFedep, fullExclude, MP_EXCL_1_4,
                                    primaryOriginCode, excludePDG);

  // ---- Sections 4–6 (MP-excl {0..5}) ----
  auto C4 = ProcessRange_HTOFselect(T, S1_lo, S1_hi, thrBH2, thrHTOF,
                                    hasBH2edep, hasHTOFedep, fullExclude, MP_EXCL_0_5,
                                    primaryOriginCode, excludePDG);
  auto C5 = ProcessRange_HTOFselect(T, S2_lo, S2_hi, thrBH2, thrHTOF,
                                    hasBH2edep, hasHTOFedep, fullExclude, MP_EXCL_0_5,
                                    primaryOriginCode, excludePDG);
  auto C6 = ProcessRange_HTOFselect(T, S3_lo, S3_hi, thrBH2, thrHTOF,
                                    hasBH2edep, hasHTOFedep, fullExclude, MP_EXCL_0_5,
                                    primaryOriginCode, excludePDG);

  // ---- Print ----
  PrintSection("Section 1 (Full Beam) MP-excl{1,2,3,4}", C1, S1_lo, S1_hi, MP_EXCL_1_4, TOPN);
  PrintSection("Section 2 (Narrow 1) MP-excl{1,2,3,4}",  C2, S2_lo, S2_hi, MP_EXCL_1_4, TOPN);
  PrintSection("Section 3 (Narrow 2) MP-excl{1,2,3,4}",  C3, S3_lo, S3_hi, MP_EXCL_1_4, TOPN);

  PrintSection("Section 4 (Full Beam) MP-excl{0,1,2,3,4,5}", C4, S1_lo, S1_hi, MP_EXCL_0_5, TOPN);
  PrintSection("Section 5 (Narrow 1) MP-excl{0,1,2,3,4,5}",  C5, S2_lo, S2_hi, MP_EXCL_0_5, TOPN);
  PrintSection("Section 6 (Narrow 2) MP-excl{0,1,2,3,4,5}",  C6, S3_lo, S3_hi, MP_EXCL_0_5, TOPN);
}

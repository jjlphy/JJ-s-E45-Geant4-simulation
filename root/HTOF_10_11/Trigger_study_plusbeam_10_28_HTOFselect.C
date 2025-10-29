// -*- C++ -*-
// Trigger_study_plusbeam_10_28_v3.C  (2025-10-29 for jaejin; OR-pair veto + combos + root-cause)
//
// PLUS-beam (π+) version with OR-based BeamVeto over adjacent HTOF pairs.
//
// Veto variants (tiles near 17–21):
//   tight      : pair(17,18) OR pair(18,19)
//   fit        : tight       OR pair(19,20)
//   wide       : fit         OR pair(16,17)
//   ultra-wide : wide        OR pair(20,21)
//
// Sections:
//   (Sec1) BH2 4–10
//   (Sec2) BH2 4–9
//   (Sec3) BH2 5–10
//
// Trigger logic per section:
//   Beam      = (BH2 in [lo,hi])
//   Trig1     = Beam && (HTOF multiplicity >= 2)          // after exclusions
//   Trig2(*)  = Trig1 && !(BeamVeto(*))                   // tiles from "full-valid"
//
// Denominator for % : Beam (= #events with BH2 in [lo,hi])
//
// Options:
//   - excludeTilesCSV : e.g. "20,21" → 해당 타일을 멀티/비토 모두에서 제외
//   - excludePDGCSV   : e.g. "11,-11,22" → 해당 PDG 입자를 분석에서 제외(EM 제거 테스트)
//
// Build/Run:
//   root -l
//   .L Trigger_study_plusbeam_10_28_v3.C+
//   Trigger_study_plusbeam_10_28_v3("E45_Nov_beamplus_098.root","g4hyptpc",
//                                   4,10, 0.10,2.0, 5.0,10.0, "", "");
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

namespace TSP2_ORX {

static const int    kNBH2Seg   = 15;     // 0..14
static const double kBH2_x0    = -10.0;  // [mm]
static const double kBH2_pitch = 14.0;   // [mm]
static const double kBH2_total = kNBH2Seg * kBH2_pitch;
static const int    kNHTOF     = 34;     // 0..33

// PDG helpers
static inline bool IsEM(int pdg){ return (pdg==11 || pdg==-11 || pdg==22); }
static inline bool IsPiPlus(int pdg){ return (pdg==+211); }

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
static inline double pct(long long a, long long b){
  return (b>0)? (100.0*double(a)/double(b)):0.0;
}

// parse "a,b,c" → {a,b,c}
static std::set<int> ParseCSVInt(const char* csv){
  std::set<int> s; if(!csv) return s;
  std::stringstream ss(csv); std::string tok;
  while(std::getline(ss, tok, ',')){ if(tok.empty()) continue; s.insert(std::stoi(tok)); }
  return s;
}
// parse "a,b,c" → PDG set
static std::set<int> ParseCSVPDG(const char* csv){ return ParseCSVInt(csv); }

static inline bool PairHit(const std::set<int>& tiles, int a, int b){
  return tiles.count(a) && tiles.count(b);
}

// Veto variants (PLUS-beam, OR over adjacent pairs)
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

// set → comma key "17,21,22"
static std::string KeyFromSet(const std::set<int>& tiles){
  if(tiles.empty()) return std::string("-");
  std::ostringstream os; bool first=true;
  for(int t: tiles){ if(!first) os<<","; first=false; os<<t; }
  return os.str();
}

struct Counts {
  long long N_total=0;
  long long N_beam =0;
  long long N_trig1=0;
  long long N_trig2_tight=0, N_trig2_fit=0, N_trig2_wide=0, N_trig2_ultra=0;
  long long N_veto_tight=0,  N_veto_fit=0,  N_veto_wide=0,  N_veto_ultra=0;

  // presence (diagnostic)
  long long have16=0, have17=0, have18=0, have19=0, have20=0, have21=0;

  // survivors combos
  std::map<std::string,long long> Comb_tight, Comb_fit, Comb_wide, Comb_ultra;

  // root-cause tallies per variant
  long long RC_allPiPri_tight=0, RC_hasEM_tight=0, RC_hasSec_tight=0;
  long long RC_allPiPri_fit=0,   RC_hasEM_fit=0,   RC_hasSec_fit=0;
  long long RC_allPiPri_wide=0,  RC_hasEM_wide=0,  RC_hasSec_wide=0;
  long long RC_allPiPri_ultra=0, RC_hasEM_ultra=0, RC_hasSec_ultra=0;

  // ΣEdep over tiles (per-variant) min/max
  double sumE_min_tight=+1e99, sumE_max_tight=-1e99;
  double sumE_min_fit=+1e99,   sumE_max_fit=-1e99;
  double sumE_min_wide=+1e99,  sumE_max_wide=-1e99;
  double sumE_min_ultra=+1e99, sumE_max_ultra=-1e99;
};

struct CauseFlags {
  bool allPiPrimary=true;  // 시작은 true, 위배되면 false
  bool hasEM=false;
  bool hasSecondary=false;
};

// try to infer primary-origin code (non-negative UniqueID mode)
static int InferPrimaryOriginCode(TTree* T){
  std::vector<TParticle>* HTOF=nullptr; T->SetBranchAddress("HTOF",&HTOF);
  std::map<int,long long> freq;
  const Long64_t N = std::min<Long64_t>(T->GetEntries(), 20000); // sample up to 20k
  for(Long64_t ie=0; ie<N; ++ie){
    T->GetEntry(ie);
    if(!HTOF) continue;
    for(const auto& p : *HTOF){
      int uid = p.GetUniqueID();
      if(uid>=0) freq[uid]++;
    }
  }
  if(freq.empty()) return 0; // fallback
  int best=-1; long long cnt=-1;
  for(const auto& kv:freq){ if(kv.second>cnt){ cnt=kv.second; best=kv.first; } }
  if(best<0) best=0;
  return best;
}

static Counts ProcessRange(TTree* T,
                           int bh2_lo, int bh2_hi,
                           double thrBH2, double thrHTOF,
                           bool hasBH2edep, bool hasHTOFedep,
                           const std::set<int>& excludeTiles,
                           const std::set<int>& excludePDG,
                           int primaryOriginCode)
{
  std::vector<TParticle>* BH2=nullptr;   T->SetBranchAddress("BH2",&BH2);
  std::vector<TParticle>* HTOF=nullptr;  T->SetBranchAddress("HTOF",&HTOF);
  std::vector<double>* BH2_edep=nullptr; if(hasBH2edep)  T->SetBranchAddress("BH2_edep",&BH2_edep);
  std::vector<double>* HTOF_edep=nullptr;if(hasHTOFedep) T->SetBranchAddress("HTOF_edep",&HTOF_edep);

  Counts C;
  const Long64_t N=T->GetEntries();
  for(Long64_t ie=0; ie<N; ++ie){
    T->GetEntry(ie); C.N_total++;

    // ---- BH2 valid
    std::map<int,double> bh2E;
    if(BH2){
      for(size_t i=0;i<BH2->size();++i){
        const TParticle& p = BH2->at(i);
        if(excludePDG.count(p.GetPdgCode())) continue;
        int sid = MapBH2_WorldToSeg(p.Vx(),p.Vy(),p.Vz());
        if(0<=sid && sid<kNBH2Seg){
          const double ed = (hasBH2edep && BH2_edep && i<BH2_edep->size()) ? BH2_edep->at(i) : p.GetWeight();
          bh2E[sid]+=ed;
        }
      }
    }
    std::set<int> bh2Valid;
    for(const auto& kv:bh2E) if(kv.second>=thrBH2) bh2Valid.insert(kv.first);
    bool beam=false; for(int s:bh2Valid){ if(bh2_lo<=s && s<=bh2_hi){ beam=true; break; } }
    if(beam) C.N_beam++;

    // ---- HTOF full-valid tiles + per-event cause flags
    std::map<int,double> htofE;
    CauseFlags cf; // starts allPiPrimary=true
    double sumE_thisEvent = 0.0;

    if(HTOF){
      for(size_t i=0;i<HTOF->size();++i){
        const TParticle& p = HTOF->at(i);
        const int pdg = p.GetPdgCode();
        if(excludePDG.count(pdg)) continue;                 // PDG 배제
        const int tid = p.GetStatusCode();
        if(!(0<=tid && tid<kNHTOF)) continue;
        if(excludeTiles.count(tid)) continue;               // 타일 배제

        const double ed = (hasHTOFedep && HTOF_edep && i<HTOF_edep->size()) ? HTOF_edep->at(i) : p.GetWeight();
        htofE[tid]+=ed;

        // cause flags (per-hit update):
        if(!IsPiPlus(pdg) || p.GetUniqueID()!=primaryOriginCode) cf.allPiPrimary=false;
        if(IsEM(pdg)) cf.hasEM=true;
        if(p.GetUniqueID()!=primaryOriginCode) cf.hasSecondary=true;

        sumE_thisEvent += ed;
      }
    }

    std::set<int> tiles;
    for(const auto& kv:htofE) if(kv.second>=thrHTOF) tiles.insert(kv.first);

    // presence diag
    if(tiles.count(16)) C.have16++;
    if(tiles.count(17)) C.have17++;
    if(tiles.count(18)) C.have18++;
    if(tiles.count(19)) C.have19++;
    if(tiles.count(20)) C.have20++;
    if(tiles.count(21)) C.have21++;

    const int mult = (int)tiles.size();

    if(beam && mult>=2){
      C.N_trig1++;
      const auto v = EvalVeto_ORpairs(tiles);
      const std::string key = KeyFromSet(tiles);

      auto upd = [&](bool survived, long long& n_surv, long long& n_veto,
                     std::map<std::string,long long>& comb,
                     long long& rc_allPiPri, long long& rc_hasEM, long long& rc_hasSec,
                     double& smin, double& smax){
        if(survived){
          n_surv++;
          comb[key]++;
          // root-cause
          if(cf.allPiPrimary) rc_allPiPri++;
          if(cf.hasEM)       rc_hasEM++;
          if(cf.hasSecondary)rc_hasSec++;
          // sumE range
          if(sumE_thisEvent<smin) smin=sumE_thisEvent;
          if(sumE_thisEvent>smax) smax=sumE_thisEvent;
        }else{
          n_veto++;
        }
      };

      upd(!v.tight, C.N_trig2_tight, C.N_veto_tight, C.Comb_tight,
          C.RC_allPiPri_tight, C.RC_hasEM_tight, C.RC_hasSec_tight,
          C.sumE_min_tight, C.sumE_max_tight);

      upd(!v.fit,   C.N_trig2_fit,   C.N_veto_fit,   C.Comb_fit,
          C.RC_allPiPri_fit,   C.RC_hasEM_fit,   C.RC_hasSec_fit,
          C.sumE_min_fit,   C.sumE_max_fit);

      upd(!v.wide,  C.N_trig2_wide,  C.N_veto_wide,  C.Comb_wide,
          C.RC_allPiPri_wide,  C.RC_hasEM_wide,  C.RC_hasSec_wide,
          C.sumE_min_wide,  C.sumE_max_wide);

      upd(!v.ultra, C.N_trig2_ultra, C.N_veto_ultra, C.Comb_ultra,
          C.RC_allPiPri_ultra, C.RC_hasEM_ultra, C.RC_hasSec_ultra,
          C.sumE_min_ultra, C.sumE_max_ultra);
    }
  }
  return C;
}

static void PrintCombos(const char* head,
                        const std::map<std::string,long long>& M)
{
  std::vector<std::pair<std::string,long long>> v(M.begin(), M.end());
  std::sort(v.begin(), v.end(),
            [](auto& a, auto& b){
              if(a.second!=b.second) return a.second>b.second;
              return a.first<b.first;
            });
  std::cout<<head<<"\n";
  if(v.empty()){ std::cout<<"  (none)\n"; return; }
  for(const auto& kv: v){
    std::cout<<"  comb("<<kv.first<<"): "<<kv.second<<"\n";
  }
}

static void PrintSection(const char* title, const Counts& C, int lo, int hi){
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

  std::cout<<"[Tile presence among Beam events (diagnostic)]\n";
  std::cout<<"  have16="<<C.have16<<", have17="<<C.have17<<", have18="<<C.have18
           <<", have19="<<C.have19<<", have20="<<C.have20<<", have21="<<C.have21<<"\n";

  // survivors combos
  auto PrintCombos = [](const char* head, const std::map<std::string,long long>& M){
    std::vector<std::pair<std::string,long long>> v(M.begin(), M.end());
    std::sort(v.begin(), v.end(),
              [](auto& a, auto& b){
                if(a.second!=b.second) return a.second>b.second;
                return a.first<b.first;
              });
    std::cout<<head<<"\n";
    if(v.empty()){ std::cout<<"  (none)\n"; return; }
    for(const auto& kv: v) std::cout<<"  comb("<<kv.first<<"): "<<kv.second<<"\n";
  };
  PrintCombos("== Survived combinations (tight) ==", C.Comb_tight);
  PrintCombos("== Survived combinations (fit) ==",   C.Comb_fit);
  PrintCombos("== Survived combinations (wide) ==",  C.Comb_wide);
  PrintCombos("== Survived combinations (ultra) ==", C.Comb_ultra);

  // root-cause tallies
  auto PrintRC = [&](const char* tag,
                     long long a,long long b,long long c,
                     double smin,double smax){
    std::cout<<"[Root-cause tallies | "<<tag<<"]\n";
    std::cout<<"  ALL_PI_PRIMARY : "<<a<<"\n";
    std::cout<<"  HAS_EM         : "<<b<<"\n";
    std::cout<<"  HAS_SECONDARY  : "<<c<<"\n";
    if(smin<=smax){
      std::cout<<"  ΣEdep over tiles  : min="<<std::setprecision(3)<<std::fixed
               <<smin<<" MeV, max="<<smax<<" MeV\n";
    }else{
      std::cout<<"  ΣEdep over tiles  : (no survivors)\n";
    }
  };
  PrintRC("tight", C.RC_allPiPri_tight, C.RC_hasEM_tight, C.RC_hasSec_tight,
          C.sumE_min_tight, C.sumE_max_tight);
  PrintRC("fit  ", C.RC_allPiPri_fit,   C.RC_hasEM_fit,   C.RC_hasSec_fit,
          C.sumE_min_fit,   C.sumE_max_fit);
  PrintRC("wide ", C.RC_allPiPri_wide,  C.RC_hasEM_wide,  C.RC_hasSec_wide,
          C.sumE_min_wide,  C.sumE_max_wide);
  PrintRC("ultra", C.RC_allPiPri_ultra, C.RC_hasEM_ultra, C.RC_hasSec_ultra,
          C.sumE_min_ultra, C.sumE_max_ultra);
}

} // namespace TSP2_ORX


void Trigger_study_plusbeam_10_28_v3(const char* filename,
                                     const char* treename="g4hyptpc",
                                     int bh2_lo=4, int bh2_hi=10,
                                     double mipFrac=0.10,
                                     double mipMeVperCm=2.0,
                                     double BH2_thickness_mm=5.0,
                                     double HTOF_thickness_mm=10.0,
                                     const char* excludeTilesCSV="",   // e.g. "20,21"
                                     const char* excludePDGCSV="" )    // e.g. "11,-11,22"
{
  using namespace TSP2_ORX;

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
  const auto excludeTiles = ParseCSVInt(excludeTilesCSV);
  const auto excludePDG   = ParseCSVPDG(excludePDGCSV);

  // infer primary-origin code (fallback to 0)
  const int primaryOriginCode = InferPrimaryOriginCode(T);
  std::cout<<std::fixed<<std::setprecision(3)
           <<"[INFO] BH2 thr="<<thrBH2<<" MeV, HTOF thr="<<thrHTOF<<" MeV\n"
           <<"[INFO] primaryOriginCode = "<<primaryOriginCode<<"\n";
  if(!excludeTiles.empty()){
    std::cout<<"[INFO] Excluding HTOF tiles: ";
    bool first=true; for(int t:excludeTiles){ if(!first) std::cout<<","; first=false; std::cout<<t; }
    std::cout<<"\n";
  }
  if(!excludePDG.empty()){
    std::cout<<"[INFO] Excluding PDGs: ";
    bool first=true; for(int p:excludePDG){ if(!first) std::cout<<","; first=false; std::cout<<p; }
    std::cout<<"\n";
  }

  // sections
  const int s1_lo = bh2_lo, s1_hi = bh2_hi; // Full
  const int s2_lo = 4,      s2_hi = 9;      // Narrow 1
  const int s3_lo = 5,      s3_hi = 10;     // Narrow 2

  auto C1 = ProcessRange(T, s1_lo, s1_hi, thrBH2, thrHTOF,
                         hasBH2edep, hasHTOFedep,
                         excludeTiles, excludePDG, primaryOriginCode);
  auto C2 = ProcessRange(T, s2_lo, s2_hi, thrBH2, thrHTOF,
                         hasBH2edep, hasHTOFedep,
                         excludeTiles, excludePDG, primaryOriginCode);
  auto C3 = ProcessRange(T, s3_lo, s3_hi, thrBH2, thrHTOF,
                         hasBH2edep, hasHTOFedep,
                         excludeTiles, excludePDG, primaryOriginCode);

  PrintSection("Section 1 (Full Beam)",   C1, s1_lo, s1_hi);
  PrintSection("Section 2 (Narrow 1)",    C2, s2_lo, s2_hi);
  PrintSection("Section 3 (Narrow 2)",    C3, s3_lo, s3_hi);
}

// -*- C++ -*-
// Trigger_study_plusreaction_10_28_mpOnly.C
// (2025-10-29 for jaejin; plus-beam OR-pair veto, reaction-only, MP-only exclusions)
//
// 목적: π+ (plus-beam) 리액션 하전입자(Origin tag)만으로 HTOF 기반 트리거/빔비토 요약
//
// BeamVeto (OR over adjacent HTOF pairs; copy-no):
//   tight      : (17,18) OR (18,19)
//   fit        : tight    OR (19,20)
//   wide       : fit      OR (16,17)
//   ultra-wide : wide     OR (20,21)
//
// Sections:
//   Sec1: BH2 4–10  , MP-only excl {1,2,3,4}
//   Sec2: BH2 4–9   , MP-only excl {1,2,3,4}
//   Sec3: BH2 5–10  , MP-only excl {1,2,3,4}
//   Sec4: BH2 4–10  , MP-only excl {0,1,2,3,4,5}
//   Sec5: BH2 4–9   , MP-only excl {0,1,2,3,4,5}
//   Sec6: BH2 5–10  , MP-only excl {0,1,2,3,4,5}
//
// 정의 / 분모:
//   Beam = BH2 in-range (옵션: 타겟터치 필터 통과)   ⇒ %는 모두 Beam 기준
//   Overkill(%) = 100 × (BeamVeto / Beam)
//
// Reaction-only:
//   HTOF의 vector<TParticle>에서 UniqueID==ORI_FORCED2PI_CHILD 만 사용.
//   (초기 2000 evt에서 해당 태그가 안 보이면 경고 후 “모든 하전” 폴백)
//
// 사용 예:
//   root -l
//   .L Trigger_study_plusreaction_10_28_mpOnly.C+
//   Trigger_study_plusreaction_10_28_mpOnly("../rootfile/E45_Nov_pipluspipn_105.root","g4hyptpc",4,10, 0.10,2.0, 5.0,10.0,/*excludeTilesCSV*/"",/*requireTargetHit*/true, /*targetBranch*/"Target", /*targetThickness_mm*/50.0);
//                                            4,10, 0.10,2.0, 5.0,10.0,
//                                            4,10, 0.10,2.0, 5.0,10.0,/*excludeTilesCSV*/"",/*requireTargetHit*/true, /*targetBranch*/"Target", /*targetThickness_mm*/50.0);
//*requireTargetHit*/true, /*targetBranch*/"Target", /*targetThickness_mm*/50.0);
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
#include <cmath>
#include <sstream>
#include <algorithm>

namespace TSPR {

static const int    kNBH2Seg   = 15;     // 0..14
static const double kBH2_x0    = -10.0;  // [mm]
static const double kBH2_pitch = 14.0;   // [mm]
static const double kBH2_total = kNBH2Seg * kBH2_pitch;
static const int    kNHTOF     = 34;     // 0..33

// 프로젝트의 Origin 코드와 동기화 필요. 여기서는 다음과 같이 가정:
static const int ORI_FORCED2PI_CHILD = 1;   // 리액션 자식(강제 2π)

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

static inline bool IsChargedByPDG(int pdg){
  if(pdg==0) return true; // SD에서 neutral 미기록 가정
  const int a=std::abs(pdg);
  return (a==211 || a==321 || a==2212 || a==11 || a==13);
}

static inline double pct(long long a, long long b){ return (b>0)? (100.0*double(a)/double(b)) : 0.0; }

// OR-기반 비토
struct VetoFired { bool tight=false, fit=false, wide=false, ultra=false; };
static inline bool HasPair(const std::set<int>& S, int a, int b){ return S.count(a)&&S.count(b); }
static inline VetoFired EvalVetoOR(const std::set<int>& tiles){
  const bool p1617 = HasPair(tiles,16,17);
  const bool p1718 = HasPair(tiles,17,18);
  const bool p1819 = HasPair(tiles,18,19);
  const bool p1920 = HasPair(tiles,19,20);
  const bool p2021 = HasPair(tiles,20,21);
  VetoFired v;
  v.tight = (p1718 || p1819);
  v.fit   = (v.tight || p1920);
  v.wide  = (v.fit   || p1617);
  v.ultra = (v.wide  || p2021);
  return v;
}

// Target 브랜치 선택 도우미
static const char* PickTargetBranch(TTree* T, const char* preferred){
  if(preferred && T->GetBranch(preferred)) return preferred;
  if(T->GetBranch("Target"))    return "Target";
  if(T->GetBranch("LH2"))       return "LH2";
  if(T->GetBranch("TargetLH2")) return "TargetLH2";
  return nullptr;
}

struct Counts {
  long long N_total=0;     // 읽은 총 이벤트
  long long N_used =0;     // 타겟 필터 통과(옵션) 후 남은 이벤트
  long long N_beam =0;     // BH2 in [lo,hi]
  long long N_trig1=0;     // Beam && (MP>=2)  [MP는 mp-only 제외 반영]
  long long N_veto_tight=0, N_veto_fit=0, N_veto_wide=0, N_veto_ultra=0; // fired
};

// 공통 유틸: CSV → set<int>
static std::set<int> ParseCSVInt(const char* csv){
  std::set<int> s; if(!csv) return s;
  std::stringstream ss(csv); std::string tok;
  while(std::getline(ss, tok, ',')){ if(tok.empty()) continue; s.insert(std::stoi(tok)); }
  return s;
}

static Counts ProcessRange(TTree* T,
                           int bh2_lo, int bh2_hi,
                           double thrBH2, double thrHTOF, double thrTarget,
                           bool hasBH2edep, bool hasHTOFedep,
                           bool originSeen,
                           bool requireTargetHit, const char* targetBranch,
                           const std::set<int>& fullExclude,     // 멀티/비토 공통 제외
                           const std::set<int>& mpOnlyExclude)   // 멀티 전용 제외
{
  std::vector<TParticle>* BH2=nullptr;   T->SetBranchAddress("BH2",&BH2);
  std::vector<TParticle>* HTOF=nullptr;  T->SetBranchAddress("HTOF",&HTOF);
  std::vector<double>* BH2_edep=nullptr; if(hasBH2edep)  T->SetBranchAddress("BH2_edep",&BH2_edep);
  std::vector<double>* HTOF_edep=nullptr;if(hasHTOFedep) T->SetBranchAddress("HTOF_edep",&HTOF_edep);

  // 타겟 필터 (vector<TParticle> + *_edep 가정)
  std::vector<TParticle>* TGT=nullptr; std::vector<double>* TGT_edep=nullptr;
  if(requireTargetHit && targetBranch){
    T->SetBranchAddress(targetBranch, &TGT);
    std::string edn = std::string(targetBranch) + "_edep";
    if(T->GetBranch(edn.c_str())) T->SetBranchAddress(edn.c_str(), &TGT_edep);
  }

  // 선택적 tgt_touch_flag(있으면 우선 사용)
  const bool hasTgtFlag = (T->GetBranch("tgt_touch_flag") != nullptr);
  int tgt_touch_flag=1; if(hasTgtFlag) T->SetBranchAddress("tgt_touch_flag",&tgt_touch_flag);

  Counts C;
  const Long64_t N=T->GetEntries();
  for(Long64_t ie=0; ie<N; ++ie){
    T->GetEntry(ie); C.N_total++;

    // --- Target filter ---
    if(requireTargetHit){
      bool pass=false;
      if(hasTgtFlag){ pass = (tgt_touch_flag==1); }
      else if(targetBranch){
        double Etot=0.0;
        if(TGT){
          for(size_t i=0;i<TGT->size();++i){
            const double ed=(TGT_edep && i<TGT_edep->size()) ? TGT_edep->at(i) : TGT->at(i).GetWeight();
            Etot += ed;
          }
        }
        pass = (Etot >= thrTarget);
      }else{
        pass = true; // 브랜치가 없으면 필터 미적용
      }
      if(!pass) continue;
    }
    C.N_used++;

    // --- BH2 유효 세그먼트 ---
    std::map<int,double> bh2E;
    if(BH2){
      for(size_t i=0;i<BH2->size();++i){
        const TParticle& p=BH2->at(i);
        int sid=MapBH2_WorldToSeg(p.Vx(),p.Vy(),p.Vz());
        if(0<=sid && sid<kNBH2Seg){
          const double ed=(hasBH2edep && BH2_edep && i<BH2_edep->size()) ? BH2_edep->at(i) : p.GetWeight();
          bh2E[sid]+=ed;
        }
      }
    }
    std::set<int> bh2Valid;
    for(const auto& kv:bh2E) if(kv.second>=thrBH2) bh2Valid.insert(kv.first);

    bool beam=false; for(int s:bh2Valid){ if(bh2_lo<=s && s<=bh2_hi){ beam=true; break; } }
    if(!beam) continue;
    C.N_beam++;

    // --- HTOF: full-valid (비토/표시), mp-valid (멀티 전용) ---
    std::map<int,double> htofE_full, htofE_mp;
    if(HTOF){
      for(size_t i=0;i<HTOF->size();++i){
        const TParticle& p=HTOF->at(i);
        const int tid=p.GetStatusCode();
        if(!(0<=tid && tid<kNHTOF)) continue;

        // Reaction-only origin
        if(originSeen){
          if(p.GetUniqueID()!=ORI_FORCED2PI_CHILD) continue;
          if(!IsChargedByPDG(p.GetPdgCode())) continue;
        }else{
          // Fallback: 모든 하전
          if(!IsChargedByPDG(p.GetPdgCode())) continue;
        }

        if(fullExclude.count(tid)) continue; // 완전 제외

        const double ed=(hasHTOFedep && HTOF_edep && i<HTOF_edep->size()) ? HTOF_edep->at(i) : p.GetWeight();
        htofE_full[tid]+=ed;
        if(!mpOnlyExclude.count(tid)) htofE_mp[tid]+=ed; // mp-only 제외 반영
      }
    }
    std::set<int> tiles_full, tiles_mp;
    for(const auto& kv:htofE_full) if(kv.second>=thrHTOF) tiles_full.insert(kv.first);
    for(const auto& kv:htofE_mp)   if(kv.second>=thrHTOF) tiles_mp.insert(kv.first);

    const int mult = (int)tiles_mp.size();
    if(mult<2) continue;
    C.N_trig1++;

    // --- 비토는 full-valid 세트로 평가 ---
    const auto v = EvalVetoOR(tiles_full);
    if(v.tight) C.N_veto_tight++;
    if(v.fit)   C.N_veto_fit++;
    if(v.wide)  C.N_veto_wide++;
    if(v.ultra) C.N_veto_ultra++;
  }
  return C;
}

static void PrintSection(const char* title, const Counts& C, int lo, int hi,
                         const std::set<int>& mpOnlyExclude)
{
  std::cout<<"\n== "<<title<<" | BH2 "<<lo<<"–"<<hi<<" ==\n";
  std::cout<<"Total events (read)           : "<<C.N_total<<"\n";
  std::cout<<"Target-used (after filter)    : "<<C.N_used <<"\n";
  std::cout<<"Beam  (BH2 in-range)          : "<<C.N_beam<<"\n";
  std::cout<<"Trig1 (Beam && MP>=2)         : "<<C.N_trig1<<"  ("<<std::fixed<<std::setprecision(3)<<pct(C.N_trig1, C.N_beam)<<" % of Beam)\n";

  if(!mpOnlyExclude.empty()){
    std::cout<<"[MP-only exclusion for multiplicity] ";
    bool first=true; for(int t:mpOnlyExclude){ if(!first) std::cout<<","; first=false; std::cout<<t; }
    std::cout<<"\n";
  }

  auto pr=[&](const char* name, long long Nveto){
    const double overkill = pct(Nveto, C.N_beam);
    std::cout<<"BeamVeto "<<std::setw(9)<<name<<": "
             <<Nveto<<"  ("<<std::fixed<<std::setprecision(3)<<overkill<<" % of Beam)  "
             <<"| Overkill = "<<std::setprecision(3)<<overkill<<" %\n";
  };
  pr("tight", C.N_veto_tight);
  pr("fit",   C.N_veto_fit);
  pr("wide",  C.N_veto_wide);
  pr("ultra", C.N_veto_ultra);
}

} // namespace TSPR


void Trigger_study_plusreaction_10_28_mpOnly(const char* filename,
                                             const char* treename="g4hyptpc",
                                             int bh2_lo=4, int bh2_hi=10,
                                             double mipFrac=0.10,
                                             double mipMeVperCm=2.0,
                                             double BH2_thickness_mm=5.0,
                                             double HTOF_thickness_mm=10.0,
                                             const char* excludeTilesCSV="",
                                             bool   requireTargetHit=true,
                                             const char* targetBranch="Target",
                                             double targetThickness_mm=50.0)
{
  using namespace TSPR;

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

  // thresholds
  const double thrBH2   = MIPThrMeV(mipFrac,mipMeVperCm,BH2_thickness_mm);
  const double thrHTOF  = MIPThrMeV(mipFrac,mipMeVperCm,HTOF_thickness_mm);
  const double thrTgt   = MIPThrMeV(mipFrac,mipMeVperCm,targetThickness_mm);

  // full exclude: 멀티/비토 공통 제외 (필요 없으면 "")
  auto fullExclude = ParseCSVInt(excludeTilesCSV);

  // target branch 결정
  const char* tgtBr = requireTargetHit ? PickTargetBranch(T, targetBranch) : nullptr;
  if(requireTargetHit && !tgtBr && !T->GetBranch("tgt_touch_flag")){
    std::cout<<"[WARN] Target branch/flag not found; target filter disabled.\n";
  }

  // origin tag 빠른 탐색
  bool originSeen=false;
  {
    std::vector<TParticle>* HTOF=nullptr;  T->SetBranchAddress("HTOF",&HTOF);
    const Long64_t Npeek=std::min<Long64_t>(T->GetEntries(), 2000);
    for(Long64_t ie=0; ie<Npeek && !originSeen; ++ie){
      T->GetEntry(ie);
      if(!HTOF) continue;
      for(size_t i=0;i<HTOF->size();++i){
        if(HTOF->at(i).GetUniqueID()==ORI_FORCED2PI_CHILD){ originSeen=true; break; }
      }
    }
    if(!originSeen){
      std::cerr<<"[WARN] No HTOF hit has UniqueID==kForced2PiChild in first "<<Npeek
               <<" events. Falling back to 'ALL charged hits'.\n";
    }
  }

  std::cout<<std::fixed<<std::setprecision(3)
           <<"[INFO] BH2 thr="<<thrBH2<<" MeV, HTOF thr="<<thrHTOF<<" MeV"
           <<", Target thr="<<thrTgt<<" MeV | "
           <<"BH2 range (Sec1) "<<bh2_lo<<"-"<<bh2_hi
           <<", Reaction-origin "<<(originSeen? "DETECTED":"MISSING→FALLBACK")<<"\n";
  if(!fullExclude.empty()){
    std::cout<<"[INFO] FULL EXCLUDE tiles (both mult & veto): ";
    bool first=true; for(int t:fullExclude){ if(!first) std::cout<<","; first=false; std::cout<<t; }
    std::cout<<"\n";
  }

  // Sections
  const int S1_lo=bh2_lo, S1_hi=bh2_hi;
  const int S2_lo=4,      S2_hi=9;
  const int S3_lo=5,      S3_hi=10;

  // MP-only exclusion sets
  const std::set<int> MP_EXCL_1_3 = {1,2,3,4};
  const std::set<int> MP_EXCL_4_6 = {0,1,2,3,4,5};

  // Run 1–3 (mp-only {1,2,3,4})
  auto C1 = ProcessRange(T, S1_lo, S1_hi, thrBH2, thrHTOF, thrTgt,
                         hasBH2edep, hasHTOFedep,
                         originSeen,
                         requireTargetHit, tgtBr,
                         fullExclude, MP_EXCL_1_3);
  auto C2 = ProcessRange(T, S2_lo, S2_hi, thrBH2, thrHTOF, thrTgt,
                         hasBH2edep, hasHTOFedep,
                         originSeen,
                         requireTargetHit, tgtBr,
                         fullExclude, MP_EXCL_1_3);
  auto C3 = ProcessRange(T, S3_lo, S3_hi, thrBH2, thrHTOF, thrTgt,
                         hasBH2edep, hasHTOFedep,
                         originSeen,
                         requireTargetHit, tgtBr,
                         fullExclude, MP_EXCL_1_3);

  // Run 4–6 (mp-only {0,1,2,3,4,5})
  auto C4 = ProcessRange(T, S1_lo, S1_hi, thrBH2, thrHTOF, thrTgt,
                         hasBH2edep, hasHTOFedep,
                         originSeen,
                         requireTargetHit, tgtBr,
                         fullExclude, MP_EXCL_4_6);
  auto C5 = ProcessRange(T, S2_lo, S2_hi, thrBH2, thrHTOF, thrTgt,
                         hasBH2edep, hasHTOFedep,
                         originSeen,
                         requireTargetHit, tgtBr,
                         fullExclude, MP_EXCL_4_6);
  auto C6 = ProcessRange(T, S3_lo, S3_hi, thrBH2, thrHTOF, thrTgt,
                         hasBH2edep, hasHTOFedep,
                         originSeen,
                         requireTargetHit, tgtBr,
                         fullExclude, MP_EXCL_4_6);

  // Print
  PrintSection("Section 1 (Full Beam) MP-excl{1,2,3,4}",      C1, S1_lo, S1_hi, MP_EXCL_1_3);
  PrintSection("Section 2 (Narrow 1)  MP-excl{1,2,3,4}",      C2, S2_lo, S2_hi, MP_EXCL_1_3);
  PrintSection("Section 3 (Narrow 2)  MP-excl{1,2,3,4}",      C3, S3_lo, S3_hi, MP_EXCL_1_3);
  PrintSection("Section 4 (Full Beam) MP-excl{0,1,2,3,4,5}",  C4, S1_lo, S1_hi, MP_EXCL_4_6);
  PrintSection("Section 5 (Narrow 1)  MP-excl{0,1,2,3,4,5}",  C5, S2_lo, S2_hi, MP_EXCL_4_6);
  PrintSection("Section 6 (Narrow 2)  MP-excl{0,1,2,3,4,5}",  C6, S3_lo, S3_hi, MP_EXCL_4_6);
}

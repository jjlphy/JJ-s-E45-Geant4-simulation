// -*- C++ -*-
// Trigger_study_plusreaction_10_28.C  (2025-10-28 for jaejin)
//
// Purpose: π+ (plus-beam) 설정에서 "리액션 하전입자(강제 2π 자식)"만으로
//          HTOF 기반 트리거/빔비토 동작을 섹션별로 요약.
// Sections:
//   (Sec1) BH2 4–10
//   (Sec2) BH2 4–9
//   (Sec3) BH2 5–10
//
// Logic:
//   Beam      = (BH2 in [lo,hi])
//   Trig1     = Beam && (HTOF multiplicity >= 2)         // Σedep>=thr
//   BeamVeto tight (OR): (17,18) or (18,19)
//            fit   (OR): tight or (19,20)
//            wide  (OR): fit   or (16,17)
//            ultra (OR): wide  or (20,21)
//   **Reaction-only**: HTOF hits with UniqueID==kForced2PiChild 만 사용.
//   (단, 처음 2000 evt에서 origin 태그를 못 보면 경고 후 "모든 하전" 폴백)
//
// %의 분모: Beam  (BH2 in-range)
//
// 사용 예:
//   root -l
//   .L Trigger_study_plusreaction_10_28.C+
//   Trigger_study_plusreaction_10_28("../rootfile/E45_Nov_pipluspipn_105.root","g4hyptpc");
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
#include <algorithm>

namespace TSPR {

static const int    kNBH2Seg   = 15;     // 0..14
static const double kBH2_x0    = -10.0;  // [mm]
static const double kBH2_pitch = 14.0;   // [mm]
static const double kBH2_total = kNBH2Seg * kBH2_pitch;
static const int    kNHTOF     = 34;     // 0..33

// EOrigin codes (프로젝트 정의와 동기화 필요)
// 여기서는 kForced2PiChild=1 로 가정
static const int ORI_FORCED2PI_CHILD = 1;

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

static inline bool IsChargedByPDG(int pdg){
  if(pdg==0) return true; // SD에서 neutral 미기록 가정
  const int a=std::abs(pdg);
  return (a==211 || a==321 || a==2212 || a==11 || a==13);
}

static inline double pct(long long a, long long b){
  return (b>0)? (100.0*double(a)/double(b)) : 0.0;
}

// π+ (plus-beam) OR-기반 Veto: pair가 "하나라도" 맞으면 veto fire
struct VetoFired {
  bool tight=false, fit=false, wide=false, ultra=false;
};

static inline bool HasPair(const std::set<int>& S, int a, int b){
  return S.count(a) && S.count(b);
}

static inline VetoFired EvalVetoOR(const std::set<int>& tiles){
  const bool p1718 = HasPair(tiles,17,18);
  const bool p1819 = HasPair(tiles,18,19);
  const bool p1920 = HasPair(tiles,19,20);
  const bool p1617 = HasPair(tiles,16,17);
  const bool p2021 = HasPair(tiles,20,21);

  VetoFired v;
  v.tight = (p1718 || p1819);
  v.fit   = (v.tight || p1920);
  v.wide  = (v.fit   || p1617);
  v.ultra = (v.wide  || p2021);
  return v;
}

struct Counts {
  long long N_total=0;     // total read
  long long N_tgt_used=0;  // if tgt_touch_flag exists: only flag==1 used; else ==N_total
  long long N_beam=0;      // BH2 in [lo,hi]
  long long N_trig1=0;     // Beam && MP>=2
  // veto fired (NOT survivors). Overkill 계산에 사용
  long long N_veto_tight=0, N_veto_fit=0, N_veto_wide=0, N_veto_ultra=0;
};

static Counts ProcessRange(TTree* T,
                           int bh2_lo, int bh2_hi,
                           double thrBH2, double thrHTOF,
                           bool hasBH2edep, bool hasHTOFedep,
                           bool originSeen, bool hasTgtFlag)
{
  std::vector<TParticle>* BH2=nullptr;   T->SetBranchAddress("BH2",&BH2);
  std::vector<TParticle>* HTOF=nullptr;  T->SetBranchAddress("HTOF",&HTOF);
  std::vector<double>* BH2_edep=nullptr; if(hasBH2edep)  T->SetBranchAddress("BH2_edep",&BH2_edep);
  std::vector<double>* HTOF_edep=nullptr;if(hasHTOFedep) T->SetBranchAddress("HTOF_edep",&HTOF_edep);
  int tgt_touch_flag=1; if(hasTgtFlag) T->SetBranchAddress("tgt_touch_flag",&tgt_touch_flag);

  Counts C;
  const Long64_t N=T->GetEntries();
  for(Long64_t ie=0; ie<N; ++ie){
    T->GetEntry(ie); C.N_total++;

    if(hasTgtFlag && tgt_touch_flag!=1) continue; // reaction only
    C.N_tgt_used++;

    // ---- BH2 유효 세그먼트
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

    bool beam=false;
    for(int s:bh2Valid){ if(bh2_lo<=s && s<=bh2_hi){ beam=true; break; } }
    if(!beam) continue;
    C.N_beam++;

    // ---- HTOF (reaction-only origin filter) & VALID tiles
    std::map<int,double> htofE;
    if(HTOF){
      for(size_t i=0;i<HTOF->size();++i){
        const TParticle& p=HTOF->at(i);
        const int tid=p.GetStatusCode();
        if(!(0<=tid && tid<kNHTOF)) continue;

        if(originSeen){
          if(p.GetUniqueID()!=ORI_FORCED2PI_CHILD) continue; // reaction charged only
          // (추가 보강) 만약 neutral이 섞일 가능성 의심 시 PDG 하전 체크
          if(!IsChargedByPDG(p.GetPdgCode())) continue;
        }else{
          // 폴백: 모든 하전(원빔 포함 가능)
          if(!IsChargedByPDG(p.GetPdgCode())) continue;
        }

        const double ed=(hasHTOFedep && HTOF_edep && i<HTOF_edep->size()) ? HTOF_edep->at(i) : p.GetWeight();
        htofE[tid]+=ed;
      }
    }
    std::set<int> tiles;
    for(const auto& kv:htofE) if(kv.second>=thrHTOF) tiles.insert(kv.first);
    const int mult=(int)tiles.size();

    // ---- Trig1 & Veto(OR) fired
    if(mult>=2){
      C.N_trig1++;
      const auto v=EvalVetoOR(tiles);
      if(v.tight) C.N_veto_tight++;
      if(v.fit)   C.N_veto_fit++;
      if(v.wide)  C.N_veto_wide++;
      if(v.ultra) C.N_veto_ultra++;
    }
  }
  return C;
}

static void PrintSection(const char* title, const Counts& C, int lo, int hi){
  std::cout<<"\n== "<<title<<" | BH2 "<<lo<<"–"<<hi<<" ==\n";
  std::cout<<"Total events (read)           : "<<C.N_total<<"\n";
  std::cout<<"Target-touched (used)         : "<<C.N_tgt_used<<"\n";
  std::cout<<"Beam  (BH2 in-range)          : "<<C.N_beam<<"\n";
  std::cout<<"Trig1 (Beam && MP>=2)         : "<<C.N_trig1<<"  ("<<std::fixed<<std::setprecision(3)<<pct(C.N_trig1, C.N_beam)<<" % of Beam)\n";

  auto pr=[&](const char* name, long long Nveto){
    const double frac = pct(Nveto, C.N_beam);
    // Overkill 정의 = 100*BeamVeto/Beam
    const double overkill = frac;
    std::cout<<"BeamVeto "<<std::setw(9)<<name<<": "
             <<Nveto<<"  ("<<std::fixed<<std::setprecision(3)<<frac<<" % of Beam)  "
             <<"| Overkill = "<<std::setprecision(3)<<overkill<<" %\n";
  };
  pr("tight", C.N_veto_tight);
  pr("fit",   C.N_veto_fit);
  pr("wide",  C.N_veto_wide);
  pr("ultra", C.N_veto_ultra);
}

} // namespace TSPR


void Trigger_study_plusreaction_10_28(const char* filename,
                                      const char* treename="g4hyptpc",
                                      int bh2_lo=4, int bh2_hi=10,
                                      double mipFrac=0.10,
                                      double mipMeVperCm=2.0,
                                      double BH2_thickness_mm=5.0,
                                      double HTOF_thickness_mm=10.0)
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
  const bool hasTgtFlag  = (T->GetBranch("tgt_touch_flag") != nullptr);

  // thresholds
  const double thrBH2  = MIPThrMeV(mipFrac,mipMeVperCm,BH2_thickness_mm);
  const double thrHTOF = MIPThrMeV(mipFrac,mipMeVperCm,HTOF_thickness_mm);

  // origin tag 존재 여부 빠른 감지
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
               <<" events. Falling back to 'ALL charged hits' (beam may be included).\n";
    }
  }

  std::cout<<std::fixed<<std::setprecision(3)
           <<"[INFO] BH2 thr="<<thrBH2<<" MeV, HTOF thr="<<thrHTOF<<" MeV | "
           <<"BH2 range (Sec1) "<<bh2_lo<<"-"<<bh2_hi
           <<", Reaction-only origin tag "<<(originSeen? "DETECTED":"MISSING→FALLBACK")<<"\n";

  // Sections
  const int s1_lo=bh2_lo, s1_hi=bh2_hi; // 입력값 (보통 4–10)
  const int s2_lo=4, s2_hi=9;
  const int s3_lo=5, s3_hi=10;

  auto C1 = ProcessRange(T, s1_lo, s1_hi, thrBH2, thrHTOF, hasBH2edep, hasHTOFedep, originSeen, hasTgtFlag);
  auto C2 = ProcessRange(T, s2_lo, s2_hi, thrBH2, thrHTOF, hasBH2edep, hasHTOFedep, originSeen, hasTgtFlag);
  auto C3 = ProcessRange(T, s3_lo, s3_hi, thrBH2, thrHTOF, hasBH2edep, hasHTOFedep, originSeen, hasTgtFlag);

  PrintSection("Section 1 (Full Beam)", C1, s1_lo, s1_hi);
  PrintSection("Section 2 (Narrow 1)",  C2, s2_lo, s2_hi);
  PrintSection("Section 3 (Narrow 2)",  C3, s3_lo, s3_hi);
}

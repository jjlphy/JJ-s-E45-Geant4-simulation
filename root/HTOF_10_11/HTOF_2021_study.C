// HTOF_2021_study.C
// 목적:
//  (기본 BH2 4–10, 변경 가능) 아래 항목들을 터미널에 출력
//   1) N(BH2 in [lo,hi])
//   2) N(BH2 in [lo,hi] && HTOF multiplicity >= 2)
//   3) (20,21) 조합을 포함한 이벤트 수  + 비율: (1)분모, (2)분모
//   4) (22,23) 조합을 포함한 이벤트 수  + 비율: (1)분모, (2)분모
//   5) {0..5} 중 ≥1 포함 && {20,21,22,23} 중 "정확히 1개만" 포함한 이벤트 수 + 비율(1),(2)
//
// 사용법:
//   root -l
//   .L HTOF_2021_study.C+
//   HTOF_2021_study("../rootfile/E45_fix_Beam_098.root","g4hyptpc",4,10,0.1,2.0,5.0,10.0);

#include "TFile.h"
#include "TTree.h"
#include "TParticle.h"
#include "TSystem.h"
#include "TInterpreter.h"
#include <vector>
#include <map>
#include <set>
#include <iostream>
#include <iomanip>
#include <cmath>

namespace H21 {

// geometry (E72-like)
static const int    kNBH2Seg   = 15;     // 0..14
static const double kBH2_x0    = -10.0;  // [mm]
static const double kBH2_pitch = 14.0;   // [mm]
static const double kBH2_total = kNBH2Seg * kBH2_pitch;
static const int    kNHTOF     = 34;     // 0..33

struct Dict {
  Dict(){ gSystem->Load("libPhysics");
          gInterpreter->GenerateDictionary("vector<TParticle>","TParticle.h;vector"); }
} _dict;

static inline double MIPThrMeV(double frac,double dEdx,double thk_mm){
  return frac * dEdx * (thk_mm*0.1); // mm→cm
}

static int MapBH2_WorldToSeg(double x, double, double){
  const double xloc = x - kBH2_x0;
  const double xmin = -0.5*kBH2_total, xmax = 0.5*kBH2_total;
  if (xloc < xmin || xloc >= xmax) return -1;
  int idx = int(std::floor((xloc - xmin)/kBH2_pitch));
  if(idx<0) idx=0; if(idx>=kNBH2Seg) idx=kNBH2Seg-1;
  return idx;
}

static inline bool InBH2Range(const std::set<int>& S, int lo, int hi){
  for(int s:S) if(lo<=s && s<=hi) return true;
  return false;
}

} // namespace H21

void HTOF_2021_study(const char* filename,
                     const char* treename="g4hyptpc",
                     int bh2_lo=4, int bh2_hi=10,
                     double mipFrac=0.1,
                     double mipMeVperCm=2.0,
                     double BH2_thickness_mm=5.0,
                     double HTOF_thickness_mm=10.0)
{
  using namespace H21;

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
  std::vector<double>* BH2_edep=nullptr;  if(hasBH2edep)  T->SetBranchAddress("BH2_edep",&BH2_edep);
  std::vector<double>* HTOF_edep=nullptr; if(hasHTOFedep) T->SetBranchAddress("HTOF_edep",&HTOF_edep);

  const double thrBH2  = MIPThrMeV(mipFrac,mipMeVperCm,BH2_thickness_mm);   // e.g. 0.10 MeV
  const double thrHTOF = MIPThrMeV(mipFrac,mipMeVperCm,HTOF_thickness_mm);  // e.g. 0.20 MeV
  std::cout<<std::fixed<<std::setprecision(3)
           <<"[INFO] thrBH2="<<thrBH2<<" MeV, thrHTOF="<<thrHTOF<<" MeV | "
           <<"BH2 range "<<bh2_lo<<"-"<<bh2_hi<<"\n";
  std::cout<<"[INFO] IDs: HTOF=StatusCode, BH2=coord-map; edep: *_edep→Weight; valid if sum>=thr.\n";

  // counters
  long long N_bh2=0;        // 1)
  long long N_bh2_m2=0;     // 2) multiplicity>=2
  long long N_2021=0;       // 3) includes tiles {20,21} (둘 다 포함)
  long long N_2223=0;       // 4) includes tiles {22,23} (둘 다 포함)
  long long N_edgeOnly1=0;  // 5) {0..5}∩HTOF ≠ ∅  &&  exactly-one-of{20,21,22,23}

  const Long64_t N=T->GetEntries();
  for(Long64_t ie=0; ie<N; ++ie){
    T->GetEntry(ie);

    // ---- BH2 valid set (sum edep ≥ thr)
    std::map<int,double> bh2E;
    if(BH2){
      for(size_t i=0;i<BH2->size();++i){
        const TParticle& p=BH2->at(i);
        int sid = MapBH2_WorldToSeg(p.Vx(),p.Vy(),p.Vz());
        if(0<=sid && sid<kNBH2Seg){
          const double ed = (hasBH2edep && BH2_edep && i<BH2_edep->size()) ? BH2_edep->at(i)
                                                                            : p.GetWeight();
          bh2E[sid] += ed;
        }
      }
    }
    std::set<int> bh2Valid;
    for(auto& kv:bh2E) if(kv.second>=thrBH2) bh2Valid.insert(kv.first);

    // ---- HTOF valid set (ID=StatusCode, sum edep ≥ thr)
    std::map<int,double> htofE;
    if(HTOF){
      for(size_t i=0;i<HTOF->size();++i){
        const TParticle& p=HTOF->at(i);
        int tid = p.GetStatusCode();
        if(0<=tid && tid<kNHTOF){
          const double ed = (hasHTOFedep && HTOF_edep && i<HTOF_edep->size()) ? HTOF_edep->at(i)
                                                                               : p.GetWeight();
          htofE[tid] += ed;
        }
      }
    }
    std::set<int> htofValid;
    for(auto& kv:htofE) if(kv.second>=thrHTOF) htofValid.insert(kv.first);
    const int mult = (int)htofValid.size();

    // --- (1) BH2 in [lo,hi]
    const bool inRange = InBH2Range(bh2Valid, bh2_lo, bh2_hi);
    if(!inRange) continue;
    N_bh2++;

    // --- (2) multiplicity >= 2
    if(mult>=2) N_bh2_m2++;

    // --- 타일 플래그
    const bool has20 = htofValid.count(20);
    const bool has21 = htofValid.count(21);
    const bool has22 = htofValid.count(22);
    const bool has23 = htofValid.count(23);

    // (3) includes both 20 and 21 (다른 타일 더 있어도 포함)
    if(has20 && has21) N_2021++;

    // (4) includes both 22 and 23
    if(has22 && has23) N_2223++;

    // (5) {0..5}∩HTOF ≠ ∅ && exactly-one-of {20,21,22,23}
    bool has0to5=false;
    for(int t=0;t<=5;++t){ if(htofValid.count(t)){ has0to5=true; break; } }
    int nEdge = (has20?1:0) + (has21?1:0) + (has22?1:0) + (has23?1:0);
    if(has0to5 && nEdge==1) N_edgeOnly1++;
  }

  auto pct = [](long long a, long long b)->double{ return (b>0)? (100.0*(double)a/(double)b):0.0; };

  // ---- print summary ----
  std::cout<<"\n========== HTOF_2021_study (BH2 "<<bh2_lo<<"-"<<bh2_hi<<") ==========\n";
  std::cout<<"(1) N[BH2 in range]                                 : "<<N_bh2<<"\n";
  std::cout<<"(2) N[BH2 in range && HTOF mult>=2]                : "<<N_bh2_m2
           <<"  ("<<std::setprecision(3)<<pct(N_bh2_m2,N_bh2)<<" % of (1))\n";

  std::cout<<"(3) includes (20 & 21)                              : "<<N_2021
           <<"  | vs (1): "<<pct(N_2021,N_bh2)<<" %"
           <<"  | vs (2): "<<pct(N_2021,N_bh2_m2)<<" %\n";

  std::cout<<"(4) includes (22 & 23)                              : "<<N_2223
           <<"  | vs (1): "<<pct(N_2223,N_bh2)<<" %"
           <<"  | vs (2): "<<pct(N_2223,N_bh2_m2)<<" %\n";

  std::cout<<"(5) has any of {0..5} AND exactly-one of {20,21,22,23} : "<<N_edgeOnly1
           <<"  | vs (1): "<<pct(N_edgeOnly1,N_bh2)<<" %"
           <<"  | vs (2): "<<pct(N_edgeOnly1,N_bh2_m2)<<" %\n";
  std::cout<<"============================================================\n";
}

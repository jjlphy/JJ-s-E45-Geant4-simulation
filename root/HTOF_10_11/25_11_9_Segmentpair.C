// -*- C++ -*-
// HTOF_pair_scan_14_27.C  (2025-11-09)
// - 입력: E45.root, g4hyptpc, BH2 윈도우(예: 4–10), 임계치 파라미터
// - 출력: Beam 사건(=BH2 in-range) 중, HTOF 인접 페어 (14,15)~(26,27)
//         가 동시에 유효(hit)인 이벤트 개수 (no-excl vs excl{2..5})
//
// 사용법:
//   root -l
//   .L HTOF_pair_scan_14_27.C+
//   HTOF_pair_scan_14_27("E45.root","g4hyptpc", 4,10, 0.10,2.0, 5.0,10.0);
//
// 메모:
//  - BH2/HTOF 브랜치: vector<TParticle>
//  - *_edep 브랜치가 있으면 우선 사용, 없으면 TParticle::Weight 사용

#include "TFile.h"
#include "TTree.h"
#include "TParticle.h"
#include "TSystem.h"
#include "TInterpreter.h"
#include <vector>
#include <map>
#include <set>
#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>

namespace HPS { // HTOF Pair Scan

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

} // namespace HPS

void HTOF_pair_scan_14_27(const char* filename,
                          const char* treename="g4hyptpc",
                          int bh2_lo=4, int bh2_hi=10,
                          double mipFrac=0.10,
                          double mipMeVperCm=2.0,
                          double BH2_thickness_mm=5.0,
                          double HTOF_thickness_mm=10.0)
{
  using namespace HPS;

  // open & tree
  TFile* f=TFile::Open(filename,"READ");
  if(!f||f->IsZombie()){ std::cerr<<"[ERR] open "<<filename<<" failed\n"; return; }
  TTree* T=(TTree*)f->Get(treename);
  if(!T){ std::cerr<<"[ERR] tree "<<treename<<" not found\n"; return; }
  if(!T->GetBranch("BH2") || !T->GetBranch("HTOF")){
    std::cerr<<"[ERR] need BH2 & HTOF (vector<TParticle>) branches\n"; return;
  }
  const bool hasBH2edep  = (T->GetBranch("BH2_edep")  != nullptr);
  const bool hasHTOFedep = (T->GetBranch("HTOF_edep") != nullptr);

  std::vector<TParticle>* BH2=nullptr;   T->SetBranchAddress("BH2",&BH2);
  std::vector<TParticle>* HTOF=nullptr;  T->SetBranchAddress("HTOF",&HTOF);
  std::vector<double>* BH2_edep=nullptr; if(hasBH2edep)  T->SetBranchAddress("BH2_edep",&BH2_edep);
  std::vector<double>* HTOF_edep=nullptr;if(hasHTOFedep) T->SetBranchAddress("HTOF_edep",&HTOF_edep);

  // thresholds
  const double thrBH2  = MIPThrMeV(mipFrac,mipMeVperCm,BH2_thickness_mm);
  const double thrHTOF = MIPThrMeV(mipFrac,mipMeVperCm,HTOF_thickness_mm);

  // 페어 목록: (14,15) .. (26,27)  → 13개
  std::vector<std::pair<int,int>> pairs;
  for(int a=14; a<=26; ++a) pairs.emplace_back(a, a+1);

  // 카운터
  long long N_total=0, N_beam=0;
  std::vector<long long> cnt_noexcl(pairs.size(), 0);
  std::vector<long long> cnt_excl(pairs.size(),   0);

  const std::set<int> EXCL_2_5 = {2,3,4,5}; // MP/유효성 판정에서 제외할 타일(요구사항)

  const Long64_t N=T->GetEntries();
  for(Long64_t ie=0; ie<N; ++ie){
    T->GetEntry(ie); N_total++;

    // ---- BH2: in-range? ----
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
    if(!beam) continue;
    N_beam++;

    // ---- HTOF 누적 → 유효 타일 집합 (no-excl / excl{2..5}) ----
    std::map<int,double> acc_noexcl, acc_excl;
    if(HTOF){
      for(size_t i=0;i<HTOF->size();++i){
        const TParticle& p=HTOF->at(i);
        const int tid = p.GetStatusCode();
        if(!(0<=tid && tid<kNHTOF)) continue;
        const double ed=(hasHTOFedep && HTOF_edep && i<HTOF_edep->size())? HTOF_edep->at(i) : p.GetWeight();
        // no-excl: 전 타일 누적
        acc_noexcl[tid] += ed;
        // excl: 2..5 제외 후 누적
        if(!EXCL_2_5.count(tid)) acc_excl[tid] += ed;
      }
    }
    std::set<int> tiles_noexcl, tiles_excl;
    for(const auto& kv:acc_noexcl) if(kv.second>=thrHTOF) tiles_noexcl.insert(kv.first);
    for(const auto& kv:acc_excl)   if(kv.second>=thrHTOF) tiles_excl.insert(kv.first);

    // ---- 페어 히트 판정 & 카운트 ----
    for(size_t ip=0; ip<pairs.size(); ++ip){
      const int a = pairs[ip].first;
      const int b = pairs[ip].second;
      if(PairHit(tiles_noexcl,a,b)) cnt_noexcl[ip]++;
      if(PairHit(tiles_excl,  a,b)) cnt_excl[ip]++;
    }
  }

  // ---- 출력 ----
  std::cout<<std::fixed<<std::setprecision(3);
  std::cout<<"[INFO] BH2 thr="<<thrBH2<<" MeV, HTOF thr="<<thrHTOF<<" MeV\n";
  std::cout<<"[INFO] BH2 window: "<<bh2_lo<<"–"<<bh2_hi<<"\n";
  std::cout<<"Total events : "<<N_total<<"\n";
  std::cout<<"Beam (BH2 in-range) : "<<N_beam<<"\n\n";

  auto pct = [&](long long x)->double{
    return (N_beam>0)? (100.0*double(x)/double(N_beam)) : 0.0;
  };

  std::cout<<"== HTOF adjacent-pair hits among Beam ==\n";
  std::cout<<" Pair    no-excl(count, %)      excl{2–5}(count, %)\n";
  for(size_t ip=0; ip<pairs.size(); ++ip){
    const auto& pr = pairs[ip];
    std::cout<<" ("<<std::setw(2)<<pr.first<<","<<std::setw(2)<<pr.second<<") : "
             <<std::setw(8)<<cnt_noexcl[ip]<<" ("<<std::setw(6)<<pct(cnt_noexcl[ip])<<"%)"
             <<"    "
             <<std::setw(8)<<cnt_excl[ip]  <<" ("<<std::setw(6)<<pct(cnt_excl[ip])  <<"%)\n";
  }
}

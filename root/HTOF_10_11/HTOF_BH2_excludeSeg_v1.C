// HTOF_BH2_excludeSeg_v1.C
// <Section3-1> BH2(3~10) & HTOF>=1 & !HTOF{1,2,3,4}
// <Section4-1> BH2(4~9)  & HTOF>=1 & !HTOF{1,2,3,4}
//
// 사용법:
//   root -l
//   .L HTOF_BH2_excludeSeg_v1.C+
//   HTOF_BH2_excludeSeg_v1("../rootfile/E45_Beam1.root","g4hyptpc");

#include "TFile.h"
#include "TTree.h"
#include "TParticle.h"
#include "TSystem.h"
#include "TInterpreter.h"
#include <vector>
#include <set>
#include <map>
#include <iostream>
#include <iomanip>
#include <cmath>

namespace HBSX {

static const int kNBH2Seg = 15;
static const int kNHTOF = 34;

struct DictGuard {
  DictGuard() {
    gSystem->Load("libPhysics");
    gInterpreter->GenerateDictionary("vector<TParticle>", "TParticle.h;vector");
  }
} _dict_guard_;

static bool HasBranch(TTree* tr, const char* bname) {
  return tr && tr->GetBranch(bname);
}

static inline double pct(double a, double b) { return (b>0)?100.0*a/b:0.0; }

static int MapBH2_WorldToSeg(double x, double y, double z) {
  const double x0 = -10.0;
  const double pitch = 14.0;
  const double total = kNBH2Seg*pitch;
  const double xmin = -0.5*total;
  const double rel = (x - x0) - xmin;
  int idx = int(std::floor(rel/pitch));
  if (idx<0 || idx>=kNBH2Seg) return -1;
  return idx;
}

static bool BH2HitsInRange(const std::set<int>& bh2Seg, int lo, int hi){
  for(int s : bh2Seg){ if(lo<=s && s<=hi) return true; }
  return false;
}

static bool HTOFExcludes(std::set<int>& htofSeg){
  for(int ex : {1,2,3,4}) if(htofSeg.count(ex)) return false;
  return true;
}

struct PairSummary {
  std::map<std::pair<int,int>, long long> count;
  long long others = 0;
  PairSummary(){ for(int i=18;i<=24;++i) count[{i,i+1}] = 0; }

  void Fill(const std::set<int>& tiles){
    if(tiles.size()!=2) return;
    auto it = tiles.begin(); int a=*it; ++it; int b=*it;
    if(a>b) std::swap(a,b);
    auto f=count.find({a,b});
    if(f!=count.end()) f->second++;
    else others++;
  }

  void Print(const char* title, long long base){
    std::cout<<"\n<"<<title<<" HTOF==2 pair summary>\n";
    for(auto &kv: count){
      std::cout<<"  pair("<<kv.first.first<<","<<kv.first.second<<") : "
               <<kv.second<<"  ("<<std::fixed<<std::setprecision(3)
               <<pct(kv.second, base)<<" % of [BH2 & HTOF>=2])\n";
    }
    std::cout<<"  others : "<<others<<"\n";
  }
};

} // namespace

void HTOF_BH2_excludeSeg_v1(const char* filename="E45.root",
                            const char* treename="g4hyptpc")
{
  using namespace HBSX;
  TFile* f=TFile::Open(filename,"READ");
  if(!f||f->IsZombie()){std::cerr<<"[ERR] cannot open "<<filename<<"\n";return;}
  TTree* T=(TTree*)f->Get(treename);
  if(!T){std::cerr<<"[ERR] tree not found\n";return;}

  std::vector<TParticle>* BH2=nullptr;
  std::vector<TParticle>* HTOF=nullptr;
  std::vector<int>* HTOF_copyNo=nullptr;

  const bool hasHTOFvec  = HasBranch(T,"HTOF");
  const bool hasHTOFcopy = HasBranch(T,"HTOF_copyNo");
  T->SetBranchAddress("BH2",&BH2);
  if(hasHTOFvec)  T->SetBranchAddress("HTOF",&HTOF);
  if(hasHTOFcopy) T->SetBranchAddress("HTOF_copyNo",&HTOF_copyNo);

  // Counters
  long long N_total=0;
  long long S31_BH2sel=0, S31_ge1=0, S31_ge2=0, S31_eq2=0;
  long long S41_BH2sel=0, S41_ge1=0, S41_ge2=0, S41_eq2=0;
  PairSummary P31, P41;

  const Long64_t N=T->GetEntries();
  for(Long64_t i=0;i<N;++i){
    T->GetEntry(i);
    N_total++;

    // ---- BH2 segment set ----
    std::set<int> bh2Seg;
    if(BH2){
      for(auto &p:*BH2){
        int s=MapBH2_WorldToSeg(p.Vx(),p.Vy(),p.Vz());
        if(s>=0&&s<kNBH2Seg) bh2Seg.insert(s);
      }
    }
    const bool pass3_10 = BH2HitsInRange(bh2Seg,3,10);
    const bool pass4_9  = BH2HitsInRange(bh2Seg,4,9);
    if(!pass3_10 && !pass4_9) continue; // skip rest if neither

    // ---- HTOF segment set ----
    std::set<int> htofSeg;
    if(hasHTOFcopy && HTOF_copyNo){
      for(int cn:*HTOF_copyNo) if(cn>=0&&cn<kNHTOF) htofSeg.insert(cn);
    }else if(hasHTOFvec && HTOF){
      for(auto &p:*HTOF){
        int cn=p.GetStatusCode();
        if(cn>=0&&cn<kNHTOF) htofSeg.insert(cn);
      }
    }
    if(htofSeg.empty()) continue;
    if(!HTOFExcludes(htofSeg)) continue; // skip if hits 1,2,3,4

    const int mult=htofSeg.size();

    // ---- Section3-1 ----
    if(pass3_10){
      S31_BH2sel++;
      if(mult>=1) S31_ge1++;
      if(mult>=2) S31_ge2++;
      if(mult==2){ S31_eq2++; P31.Fill(htofSeg); }
    }

    // ---- Section4-1 ----
    if(pass4_9){
      S41_BH2sel++;
      if(mult>=1) S41_ge1++;
      if(mult>=2) S41_ge2++;
      if(mult==2){ S41_eq2++; P41.Fill(htofSeg); }
    }
  }

  // ---- 결과 출력 ----
  std::cout<<"\n============================\n";
  std::cout<<"<Section3-1> BH2(3-10) & HTOF>=1 & !HTOF(1,2,3,4)\n";
  std::cout<<"  Total events: "<<N_total<<"\n";
  std::cout<<"  BH2-selected : "<<S31_BH2sel<<"\n";
  std::cout<<"  BH2&HTOF>=1  : "<<S31_ge1<<"\n";
  std::cout<<"  BH2&HTOF>=2  : "<<S31_ge2<<" ("<<pct(S31_ge2,S31_BH2sel)<<" % of BH2)\n";
  std::cout<<"  BH2&HTOF==2  : "<<S31_eq2<<"\n";
  P31.Print("Section3-1", S31_ge2);

  std::cout<<"\n----------------------------\n";
  std::cout<<"<Section4-1> BH2(4-9) & HTOF>=1 & !HTOF(1,2,3,4)\n";
  std::cout<<"  Total events: "<<N_total<<"\n";
  std::cout<<"  BH2-selected : "<<S41_BH2sel<<"\n";
  std::cout<<"  BH2&HTOF>=1  : "<<S41_ge1<<"\n";
  std::cout<<"  BH2&HTOF>=2  : "<<S41_ge2<<" ("<<pct(S41_ge2,S41_BH2sel)<<" % of BH2)\n";
  std::cout<<"  BH2&HTOF==2  : "<<S41_eq2<<"\n";
  P41.Print("Section4-1", S41_ge2);
  std::cout<<"============================\n";
}

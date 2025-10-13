// HTOF_searching_others.C  (합본: OTHERS ∪ (HTOF mult>2), BH2 1D만 그림, HTOF 조합 콘솔 출력)
// 목적: BH2 [lo,hi] 구간 통과 이벤트 중 (A) HTOF==2 others, (B) HTOF>2 를 합쳐서 선택.
//       선택 이벤트에 대해 (1) BH2 hit pattern (1D)만 그리기, (2) 각 이벤트의 HTOF 조합을 터미널에 출력.
// 정책(기존 논리 유지):
//   * HTOF 타일 ID = TParticle::StatusCode()
//   * BH2 세그 ID  = 좌표→세그 매핑(MapBH2_WorldToSeg)
//   * edep = *_edep 브랜치 우선, 없으면 TParticle::Weight()
//   * 컷: 0.1 MIP (기본), dE/dx=2 MeV/cm, BH2=5 mm(0.10 MeV), HTOF=10 mm(0.20 MeV)

#include "TFile.h"
#include "TTree.h"
#include "TParticle.h"
#include "TSystem.h"
#include "TInterpreter.h"
#include "TH1I.h"
#include "TCanvas.h"
#include "TString.h"
#include <vector>
#include <map>
#include <set>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>

namespace HSO {

// geometry (E72-like)
static const int    kNBH2Seg   = 15;     // 0..14
static const double kBH2_x0    = -10.0;  // [mm]
static const double kBH2_pitch = 14.0;   // [mm]
static const double kBH2_total = kNBH2Seg * kBH2_pitch;
static const int    kNHTOF     = 34;     // 0..33

// ROOT dict for vector<TParticle>
struct DictGuard {
  DictGuard(){
    gSystem->Load("libPhysics");
    gInterpreter->GenerateDictionary("vector<TParticle>","TParticle.h;vector");
  }
};
static DictGuard _dg;

static inline double MIPThrMeV(double frac,double dEdx,double thk_mm){
  return frac * dEdx * (thk_mm*0.1); // mm → cm
}

// 월드좌표 → BH2 세그(0..14), 실패 -1
static int MapBH2_WorldToSeg(double x, double, double){
  const double xloc = x - kBH2_x0;
  const double xmin = -0.5*kBH2_total;
  const double xmax =  0.5*kBH2_total;
  if (xloc < xmin || xloc >= xmax) return -1;
  int idx = int(std::floor((xloc - xmin)/kBH2_pitch));
  if(idx<0) idx=0; if(idx>=kNBH2Seg) idx=kNBH2Seg-1;
  return idx;
}

// others 판단: HTOF==2 & 인접(18–25) 7쌍이 아니면 others
static inline bool IsOthersPair(const std::set<int>& tiles){
  if((int)tiles.size()!=2) return false;
  auto it=tiles.begin(); int a=*it; ++it; int b=*it;
  if(a>b) std::swap(a,b);
  bool isAdj = (b==a+1) || (a==0 && b==(kNHTOF-1));
  bool isAdj7 = isAdj && (a>=18 && b<=25);
  return !isAdj7;
}

// set<int> → "(a,b,c)" 문자열
static std::string SetToTuple(const std::set<int>& S){
  std::string out="(";
  bool first=true;
  for(int v:S){ if(!first) out += ","; first=false; out += std::to_string(v); }
  out += ")";
  return out;
}

} // namespace HSO


// ================= PUBLIC =================
void HTOF_searching_others(const char* filename,
                           const char* treename="g4hyptpc",
                           int bh2_lo=4, int bh2_hi=10,
                           double mipFrac=0.1,
                           double mipMeVperCm=2.0,
                           double BH2_thickness_mm=5.0,
                           double HTOF_thickness_mm=10.0,
                           bool save=true,
                           const char* tag="BH2_4_10")
{
  using namespace HSO;

  // open
  TFile* f=TFile::Open(filename,"READ");
  if(!f||f->IsZombie()){ std::cerr<<"[ERR] open "<<filename<<" failed\n"; return; }
  TTree* T=(TTree*)f->Get(treename);
  if(!T){ std::cerr<<"[ERR] tree "<<treename<<" not found\n"; return; }

  if(!T->GetBranch("BH2") || !T->GetBranch("HTOF")){
    std::cerr<<"[ERR] need BH2 & HTOF branches (vector<TParticle>)\n"; return;
  }
  const bool hasBH2edep  = (T->GetBranch("BH2_edep")  != nullptr);
  const bool hasHTOFedep = (T->GetBranch("HTOF_edep") != nullptr);

  // branches
  std::vector<TParticle>* BH2=nullptr;   T->SetBranchAddress("BH2",&BH2);
  std::vector<TParticle>* HTOF=nullptr;  T->SetBranchAddress("HTOF",&HTOF);
  std::vector<double>* BH2_edep=nullptr;  if(hasBH2edep)  T->SetBranchAddress("BH2_edep",&BH2_edep);
  std::vector<double>* HTOF_edep=nullptr; if(hasHTOFedep) T->SetBranchAddress("HTOF_edep",&HTOF_edep);

  // thresholds
  const double thrBH2  = MIPThrMeV(mipFrac,mipMeVperCm,BH2_thickness_mm);   // 0.10 MeV
  const double thrHTOF = MIPThrMeV(mipFrac,mipMeVperCm,HTOF_thickness_mm);  // 0.20 MeV
  std::cout<<std::fixed<<std::setprecision(3)
           <<"[INFO] thrBH2="<<thrBH2<<" MeV, thrHTOF="<<thrHTOF<<" MeV | "
           <<"BH2 seg range "<<bh2_lo<<"-"<<bh2_hi<<"\n";
  std::cout<<"[INFO] ID: HTOF=StatusCode, BH2=coordinate mapping. edep=";
  if(hasBH2edep||hasHTOFedep) std::cout<<"*_edep or Weight\n"; else std::cout<<"Weight\n";

  // histogram: BH2 hit pattern (선택 이벤트만)
  TH1I* hBH2 = new TH1I("hBH2_sel", "BH2 hit pattern (OTHERS ∪ mult>2);BH2 seg;Events",
                        kNBH2Seg,-0.5,kNBH2Seg-0.5);
  hBH2->SetDirectory(nullptr);

  auto inBH2Range = [=](const std::set<int>& S){
    for(int s:S){ if(bh2_lo<=s && s<=bh2_hi) return true; }
    return false;
  };

  Long64_t nSelected=0;
  const Long64_t N=T->GetEntries();
  for(Long64_t ie=0; ie<N; ++ie){
    T->GetEntry(ie);

    // --- BH2 valid set (∑edep with fallback as in existing logic) ---
    std::map<int,double> bh2E; std::set<int> bh2Any;
    if(BH2){
      for(size_t i=0;i<BH2->size();++i){
        const TParticle& p=BH2->at(i);
        int sid = MapBH2_WorldToSeg(p.Vx(),p.Vy(),p.Vz());
        if(0<=sid && sid<kNBH2Seg){
          double ed = (hasBH2edep && BH2_edep && i<BH2_edep->size()) ? BH2_edep->at(i)
                                                                     : p.GetWeight();
          bh2E[sid]+=ed; bh2Any.insert(sid);
        }
      }
    }
    std::set<int> bh2Valid;
    bool anyE=false; for(auto& kv:bh2E){ if(kv.second>0){ anyE=true; break; } }
    if(anyE){ for(auto& kv:bh2E){ if(kv.second>=thrBH2) bh2Valid.insert(kv.first); } }
    else     { bh2Valid = bh2Any; }

    if(bh2Valid.empty() || !inBH2Range(bh2Valid)) continue; // BH2 구간 조건

    // --- HTOF valid set (ID=StatusCode, ∑edep with fallback as in existing logic) ---
    std::map<int,double> htofE; std::set<int> htofAny;
    if(HTOF){
      for(size_t i=0;i<HTOF->size();++i){
        const TParticle& p=HTOF->at(i);
        int tid = p.GetStatusCode();
        if(0<=tid && tid<kNHTOF){
          double ed = (hasHTOFedep && HTOF_edep && i<HTOF_edep->size()) ? HTOF_edep->at(i)
                                                                         : p.GetWeight();
          htofE[tid]+=ed; htofAny.insert(tid);
        }
      }
    }
    std::set<int> htofValid;
    bool anyEH=false; for(auto& kv:htofE){ if(kv.second>0){ anyEH=true; break; } }
    if(anyEH){ for(auto& kv:htofE){ if(kv.second>=thrHTOF) htofValid.insert(kv.first); } }
    else     { htofValid = htofAny; }

    const int m = (int)htofValid.size();
    if(m<2) continue;

    // --- 합본 선택: (m>=3) OR (m==2 && others) ---
    bool isOthers = (m==2) && IsOthersPair(htofValid);
    bool isGT2    = (m>=3);
    if(!(isOthers || isGT2)) continue;

    // 콘솔 출력: HTOF 조합
    std::cout << "HTOF " << SetToTuple(htofValid) << "\n";

    // BH2 1D 채우기 (유니크 세그 기준, 이벤트당 한 번씩 증가)
    for(int s:bh2Valid) hBH2->Fill(s);
    ++nSelected;
  }

  std::cout<<"[DONE] Selected events (OTHERS ∪ mult>2): "<<nSelected<<"\n";

  // --- draw & save ---
  TCanvas* c1=new TCanvas("cBH2_sel","BH2 hit pattern (OTHERS ∪ mult>2)",800,600);
  hBH2->Draw("hist");
  if(save){
    TString t(tag);
    c1->SaveAs("BH2pattern_selected_"+t+".png");
  }
}

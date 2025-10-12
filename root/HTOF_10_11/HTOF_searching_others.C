// HTOF_searching_others.C
// 목적: BH2 [lo,hi] 구간 통과 이벤트 중 (A) HTOF==2 others, (B) HTOF>2 를 선별하여
//       (1) BH2 hit pattern (1D), (2) HTOF hit pattern (1D), (3) BH2 vs HTOF (2D) 히스토그램 생성/저장.
// 정책(기존 논리):
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
#include "TH2I.h"
#include "TCanvas.h"
#include "TString.h"
#include <vector>
#include <map>
#include <set>
#include <iostream>
#include <iomanip>
#include <cmath>

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

// set<int> 유니크 아이템으로 1D, 2D 히스토 채우기
static inline void FillPatterns(const std::set<int>& bh2Seg,
                                const std::set<int>& htofSeg,
                                TH1I* hBH2, TH1I* hHTOF, TH2I* h2){
  if(hBH2)  for(int s:bh2Seg)  hBH2->Fill(s);
  if(hHTOF) for(int t:htofSeg) hHTOF->Fill(t);
  if(h2)    for(int s:bh2Seg) for(int t:htofSeg) h2->Fill(s,t);
}

} // namespace HSO


// ================= PUBLIC =================
void HTOF_searching_others(const char* filename,
                           const char* treename="g4hyptpc",
                           int bh2_lo=3, int bh2_hi=11,
                           double mipFrac=0.1,
                           double mipMeVperCm=2.0,
                           double BH2_thickness_mm=5.0,
                           double HTOF_thickness_mm=10.0,
                           bool save=true,
                           const char* tag="BH2_3_11")
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

  // histograms (OTHERS)
  TH1I* hBH2_O = new TH1I("hBH2_others",  "BH2 hit pattern (OTHERS);BH2 seg;Events", kNBH2Seg,-0.5,kNBH2Seg-0.5);
  TH1I* hHTOF_O= new TH1I("hHTOF_others", "HTOF hit pattern (OTHERS);HTOF tile;Events", kNHTOF,-0.5,kNHTOF-0.5);
  TH2I* h2_O   = new TH2I("hBH2vHTOF_others","BH2 vs HTOF (OTHERS);BH2 seg;HTOF tile",
                          kNBH2Seg,-0.5,kNBH2Seg-0.5, kNHTOF,-0.5,kNHTOF-0.5);

  // histograms (HTOF>2)
  TH1I* hBH2_G = new TH1I("hBH2_gt2",  "BH2 hit pattern (HTOF>2);BH2 seg;Events", kNBH2Seg,-0.5,kNBH2Seg-0.5);
  TH1I* hHTOF_G= new TH1I("hHTOF_gt2", "HTOF hit pattern (HTOF>2);HTOF tile;Events", kNHTOF,-0.5,kNHTOF-0.5);
  TH2I* h2_G   = new TH2I("hBH2vHTOF_gt2","BH2 vs HTOF (HTOF>2);BH2 seg;HTOF tile",
                          kNBH2Seg,-0.5,kNBH2Seg-0.5, kNHTOF,-0.5,kNHTOF-0.5);

  hBH2_O->SetDirectory(nullptr); hHTOF_O->SetDirectory(nullptr); h2_O->SetDirectory(nullptr);
  hBH2_G->SetDirectory(nullptr); hHTOF_G->SetDirectory(nullptr); h2_G->SetDirectory(nullptr);

  auto inBH2Range = [=](const std::set<int>& S){
    for(int s:S){ if(bh2_lo<=s && s<=bh2_hi) return true; }
    return false;
  };

  Long64_t nOthers=0, nGT2=0;
  const Long64_t N=T->GetEntries();
  for(Long64_t ie=0; ie<N; ++ie){
    T->GetEntry(ie);

    // --- BH2 valid set ---
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

    // --- HTOF valid set (ID=StatusCode) ---
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

    if(m==2){
      if(IsOthersPair(htofValid)){
        FillPatterns(bh2Valid, htofValid, hBH2_O, hHTOF_O, h2_O);
        ++nOthers;
      }
    }else{ // m>=3
      FillPatterns(bh2Valid, htofValid, hBH2_G, hHTOF_G, h2_G);
      ++nGT2;
    }
  }

  std::cout<<"[DONE] Selected events: OTHERS="<<nOthers<<", HTOF>2="<<nGT2<<"\n";

  // --- draw & save ---
  TCanvas* c1=new TCanvas("cOthers_BH2","OTHERS - BH2 hit pattern",800,600);
  hBH2_O->Draw("hist");
  TCanvas* c2=new TCanvas("cOthers_HTOF","OTHERS - HTOF hit pattern",800,600);
  hHTOF_O->Draw("hist");
  TCanvas* c3=new TCanvas("cOthers_2D","OTHERS - BH2 vs HTOF",900,700);
  h2_O->Draw("COLZ");

  TCanvas* c4=new TCanvas("cGT2_BH2","HTOF>2 - BH2 hit pattern",800,600);
  hBH2_G->Draw("hist");
  TCanvas* c5=new TCanvas("cGT2_HTOF","HTOF>2 - HTOF hit pattern",800,600);
  hHTOF_G->Draw("hist");
  TCanvas* c6=new TCanvas("cGT2_2D","HTOF>2 - BH2 vs HTOF",900,700);
  h2_G->Draw("COLZ");

  if(save){
    TString t(tag);
    c1->SaveAs("others_BH2pattern_"+t+".png");
    c2->SaveAs("others_HTOFpattern_"+t+".png");
    c3->SaveAs("others_BH2vHTOF_"+t+".png");
    c4->SaveAs("gt2_BH2pattern_"+t+".png");
    c5->SaveAs("gt2_HTOFpattern_"+t+".png");
    c6->SaveAs("gt2_BH2vHTOF_"+t+".png");
  }
}

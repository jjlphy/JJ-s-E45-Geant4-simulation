// HTOF_BH2_summary_v2.C
// (기존 논리판: HTOF ID=StatusCode, BH2 ID=좌표→세그, e-dep=Weight 사용)
// 사용법:
//   root -l
//   .L HTOF_BH2_summary_v2.C+
//   HTOF_BH2_summary_v2("../rootfile/E45_fix_Beam_098.root","g4hyptpc","3-11,4-9",true);
//
// 기능 요약:
//   - Section1: BH2(any), 추가 섹션은 "a-b,a-b,..." 문자열로 지정
//   - e-dep 컷(기본 0.1 MIP): 폴리스티렌 dE/dx=2 MeV/cm, BH2=5 mm, HTOF=10 mm
//   - e-dep 소스 우선순위: *_edep 브랜치 → TParticle::Weight()
//   - ID 규약(기존 논리 고정):
//       * HTOF tile ID = TParticle::StatusCode()
//       * BH2 seg ID   = 좌표→세그 매핑(MapBH2_WorldToSeg)
//   - HTOF==2 페어: 특수 7쌍(18–19…24–25) + others 집계, 인접(i-(i+1)) 히스토그램

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1I.h"
#include "TSystem.h"
#include "TInterpreter.h"
#include "TParticle.h"
#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"
#include <vector>
#include <set>
#include <map>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <algorithm>

namespace HBS {

// -------- geometry (E72-like) --------
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

static bool HasBranch(TTree* tr, const char* b){ return tr && tr->GetBranch(b); }

static int MapBH2_WorldToSeg(double x, double, double){
  const double xloc = x - kBH2_x0;
  const double xmin = -0.5*kBH2_total;
  const double xmax =  0.5*kBH2_total;
  if (xloc < xmin || xloc >= xmax) return -1;
  int idx = int(std::floor( (xloc - xmin)/kBH2_pitch ));
  if(idx<0) idx=0; if(idx>=kNBH2Seg) idx=kNBH2Seg-1;
  return idx;
}

static inline double pct(double a,double b){ return (b>0)?100.0*a/b:0.0; }

// ---- Pair summary (special 18..24) ----
struct PairCount {
  std::map<std::pair<int,int>, long long> counts;
  long long others=0;
  PairCount(){ for(int i=18;i<=24;++i) counts[{i,i+1}] = 0; }
  void Fill(const std::set<int>& tiles){
    if((int)tiles.size()!=2) return;
    auto it=tiles.begin(); int a=*it; ++it; int b=*it;
    if(a>b) std::swap(a,b);
    auto f = counts.find({a,b});
    if(f!=counts.end()) f->second++;
    else others++;
  }
  void PrintSummary(const char* tag, long long N_ge2){
    std::cout << "  ["<<tag<<"] HTOF==2 pair summary:\n";
    long long total=others; for(auto& kv:counts) total+=kv.second;
    for(auto& kv:counts){
      std::cout<<"    pair("<<kv.first.first<<","<<kv.first.second<<") : "
               <<kv.second<<"  ("<<std::fixed<<std::setprecision(3)
               <<pct(kv.second,N_ge2)<<" % of [BH2 & HTOF>=2])\n";
    }
    std::cout<<"    others : "<<others<<"\n";
    std::cout<<"    total(HTOF==2) counted in summary : "<<total<<"\n";
  }
};

struct StatsBox {
  long long N_total=0, N_BH2sel=0, N_ge1=0, N_ge2=0, N_eq2=0;
  PairCount pairSummary;
  TH1I* hHTOF=nullptr; TH1I* hPairAdj=nullptr;

  void InitHist(const TString& key,const TString& title){
    hHTOF=new TH1I(key,title,kNHTOF,-0.5,kNHTOF-0.5); hHTOF->SetDirectory(nullptr);
  }
  void InitPairHist(const TString& key,const TString& title){
    hPairAdj=new TH1I(key,title,kNHTOF,-0.5,kNHTOF-0.5); hPairAdj->SetDirectory(nullptr);
    for(int i=0;i<kNHTOF;++i){ int j=(i+1)%kNHTOF; hPairAdj->GetXaxis()->SetBinLabel(i+1,Form("%d-%d",i,j)); }
  }
  void FillPairHistIfAdjacent(const std::set<int>& tiles){
    if(!hPairAdj || (int)tiles.size()!=2) return;
    auto it=tiles.begin(); int a=*it; ++it; int b=*it; if(a>b) std::swap(a,b);
    if(b==a+1) hPairAdj->Fill(a);
    else if(a==0 && b==kNHTOF-1) hPairAdj->Fill(kNHTOF-1);
  }
  void Print(const char* title){
    std::cout<<"\n<"<<title<<">\n";
    std::cout<<"  Total events                            : "<<N_total<<"\n";
    std::cout<<"  BH2-selected events                     : "<<N_BH2sel
             <<"  ("<<std::fixed<<std::setprecision(3)<<pct(N_BH2sel,N_total)<<" % of Total)\n";
    std::cout<<"  BH2 & HTOF>=1                           : "<<N_ge1
             <<"  ("<<pct(N_ge1,N_total)<<" % of Total, "<<pct(N_ge1,N_BH2sel)<<" % of BH2)\n";
    std::cout<<"  BH2 & HTOF>=2                           : "<<N_ge2
             <<"  ("<<pct(N_ge2,N_total)<<" % of Total, "<<pct(N_ge2,N_BH2sel)<<" % of BH2, "<<pct(N_ge2,N_ge1)<<" % of [BH2 & HTOF>=1])\n";
    std::cout<<"  BH2 & HTOF==2                           : "<<N_eq2
             <<"  ("<<pct(N_eq2,N_total)<<" % of Total, "<<pct(N_eq2,N_BH2sel)<<" % of BH2, "<<pct(N_eq2,N_ge1)<<" % of [BH2 & HTOF>=1])\n";
    pairSummary.PrintSummary(title,N_ge2);
  }
};

static bool BH2HitsInRange(const std::set<int>& bh2Seg,int lo,int hi){
  for(int s:bh2Seg){ if(lo<=s && s<=hi) return true; } return false;
}

// ---- sections parser ----
struct SectionSpec{ bool isAny=false; int lo=0,hi=0; TString tag,key; };
static inline TString TrimWS(const TString& s){ TString t=s; t.ReplaceAll(" ",""); t.ReplaceAll("\t",""); return t; }

static std::vector<SectionSpec> ParseRanges(const char* csvIn){
  std::vector<SectionSpec> out;
  SectionSpec any; any.isAny=true; any.tag="BH2(any seg)"; any.key="BH2_any"; out.push_back(any);
  if(!csvIn) return out;
  TString csv=csvIn; csv=csv.Strip(TString::kBoth); if(csv.Length()==0) return out;
  TObjArray* toks=csv.Tokenize(","); if(!toks) return out;
  for(int i=0;i<toks->GetEntriesFast();++i){
    auto* o=toks->At(i); if(!o) continue;
    TString tok=((TObjString*)o)->GetString(); tok=TrimWS(tok); if(tok.Length()==0) continue;
    Ssiz_t d=tok.Index("-"); if(d==kNPOS){ std::cerr<<"[WARN] bad range token: "<<tok<<"\n"; continue; }
    int lo=TString(tok(0,d)).Atoi(); int hi=TString(tok(d+1,tok.Length()-d-1)).Atoi(); if(lo>hi) std::swap(lo,hi);
    if(lo<0) lo=0; if(hi>=kNBH2Seg) hi=kNBH2Seg-1;
    SectionSpec sp; sp.isAny=false; sp.lo=lo; sp.hi=hi; sp.tag.Form("BH2 seg %d-%d",lo,hi); sp.key.Form("BH2_%d_%d",lo,hi);
    out.push_back(sp);
  }
  delete toks; return out;
}

// ---- thresholds ----
static inline double MIPThrMeV(double frac,double dEdx_MeVperCm,double thick_mm){
  return frac*dEdx_MeVperCm*(thick_mm*0.1);
}

} // namespace HBS


// ================== PUBLIC ==================
void HTOF_BH2_summary_v2(const char* filename="E45.root",
                         const char* treename="g4hyptpc",
                         const char* ranges="4-9,4-11,3-8",
                         bool save=true,
                         double mipFrac=0.1,
                         double mipMeVperCm=2.0,
                         double BH2_thickness_mm=5.0,
                         double HTOF_thickness_mm=10.0)
{
  using namespace HBS;

  TFile* f=TFile::Open(filename,"READ");
  if(!f||f->IsZombie()){ std::cerr<<"[ERR] open "<<filename<<" failed\n"; return; }
  TTree* T=(TTree*)f->Get(treename);
  if(!T){ std::cerr<<"[ERR] tree "<<treename<<" not found\n"; return; }

  if(!HasBranch(T,"BH2")||!HasBranch(T,"HTOF")){
    std::cerr<<"[ERR] need branches 'BH2' and 'HTOF' (vector<TParticle>)\n"; return;
  }
  const bool hasBH2edep = HasBranch(T,"BH2_edep");
  const bool hasHTOFedep= HasBranch(T,"HTOF_edep");

  std::vector<TParticle>* BH2=nullptr; std::vector<TParticle>* HTOF=nullptr;
  std::vector<double>*   BH2_edep=nullptr; std::vector<double>* HTOF_edep=nullptr;

  T->SetBranchAddress("BH2",&BH2);
  T->SetBranchAddress("HTOF",&HTOF);
  if(hasBH2edep)  T->SetBranchAddress("BH2_edep",&BH2_edep);
  if(hasHTOFedep) T->SetBranchAddress("HTOF_edep",&HTOF_edep);

  const double thrBH2  = MIPThrMeV(mipFrac,mipMeVperCm,BH2_thickness_mm);   // 0.10 MeV
  const double thrHTOF = MIPThrMeV(mipFrac,mipMeVperCm,HTOF_thickness_mm);  // 0.20 MeV

  std::cout<<std::fixed<<std::setprecision(3);
  std::cout<<"[INFO] Energy cuts: mipFrac="<<mipFrac
           <<", dE/dx="<<mipMeVperCm<<" MeV/cm"
           <<", BH2_th="<<BH2_thickness_mm<<" mm (thr="<<thrBH2<<" MeV)"
           <<", HTOF_th="<<HTOF_thickness_mm<<" mm (thr="<<thrHTOF<<" MeV)\n";
  if(!(hasBH2edep||hasHTOFedep))
    std::cout<<"[INFO] Using embedded TParticle::Weight() as edep (no *_edep branches).\n";
  std::cout<<"[INFO] ID policy: HTOF=StatusCode, BH2=coordinate-mapping.\n";

  auto specs = ParseRanges(ranges);

  struct SectionRun{ SectionSpec spec; StatsBox box; TCanvas* c1=nullptr; TCanvas* c2=nullptr; };
  std::vector<SectionRun> runs; runs.reserve(specs.size());
  for(auto& sp:specs){
    SectionRun r; r.spec=sp;
    TString hn,ht, hp,tp;
    hn.Form("hHTOF_%s",sp.key.Data());
    ht = sp.isAny? "HTOF pattern | BH2(any seg);HTOF tile;Events"
                 : Form("HTOF pattern | %s;HTOF tile;Events",sp.tag.Data());
    r.box.InitHist(hn,ht);
    hp.Form("hPairs_%s",sp.key.Data());
    tp = sp.isAny? "Adjacent HTOF pairs (i-(i+1)) | BH2(any seg);pair;Events"
                  : Form("Adjacent HTOF pairs (i-(i+1)) | %s;pair;Events",sp.tag.Data());
    r.box.InitPairHist(hp,tp);
    runs.push_back(r);
  }

  const Long64_t N=T->GetEntries();
  for(Long64_t ie=0; ie<N; ++ie){
    T->GetEntry(ie);
    for(auto& r:runs) r.box.N_total++;

    // --- BH2: 좌표→세그, edep = *_edep or Weight ---
    std::map<int,double> bh2E; std::set<int> bh2Any;
    if(BH2){
      for(size_t i=0;i<BH2->size();++i){
        const TParticle& p=BH2->at(i);
        int sid = MapBH2_WorldToSeg(p.Vx(),p.Vy(),p.Vz());
        if(0<=sid && sid<kNBH2Seg){
          double edep = (hasBH2edep && BH2_edep && i<BH2_edep->size()) ? BH2_edep->at(i)
                                                                        : p.GetWeight();
          bh2E[sid]+=edep; bh2Any.insert(sid);
        }
      }
    }
    std::set<int> bh2Valid;
    bool anyE=false; for(auto& kv:bh2E){ if(kv.second>0){ anyE=true; break; } }
    if(anyE){ for(auto& kv:bh2E){ if(kv.second>=thrBH2) bh2Valid.insert(kv.first); } }
    else     { bh2Valid = bh2Any; }
    const bool passBH2_any = !bh2Valid.empty();

    // --- HTOF: ID = StatusCode, edep = *_edep or Weight ---
    std::map<int,double> htofE; std::set<int> htofAny;
    if(HTOF){
      for(size_t i=0;i<HTOF->size();++i){
        const TParticle& p=HTOF->at(i);
        int tid = p.GetStatusCode();
        if(0<=tid && tid<kNHTOF){
          double edep = (hasHTOFedep && HTOF_edep && i<HTOF_edep->size()) ? HTOF_edep->at(i)
                                                                           : p.GetWeight();
          htofE[tid]+=edep; htofAny.insert(tid);
        }
      }
    }
    std::set<int> htofValid;
    bool anyEH=false; for(auto& kv:htofE){ if(kv.second>0){ anyEH=true; break; } }
    if(anyEH){ for(auto& kv:htofE){ if(kv.second>=thrHTOF) htofValid.insert(kv.first); } }
    else     { htofValid = htofAny; }
    const int htMult=(int)htofValid.size();

    // --- Sections ---
    for(auto& r:runs){
      bool pass = r.spec.isAny ? passBH2_any
                               : BH2HitsInRange(bh2Valid,r.spec.lo,r.spec.hi);
      if(!pass) continue;
      r.box.N_BH2sel++;
      if(htMult>=1) r.box.N_ge1++;
      if(htMult>=2) r.box.N_ge2++;
      if(htMult==2){
        r.box.N_eq2++;
        r.box.pairSummary.Fill(htofValid);
        r.box.FillPairHistIfAdjacent(htofValid);
      }
      for(int t:htofValid) r.box.hHTOF->Fill(t);
    }
  }

  // --- Print & Draw ---
  for(auto& r:runs){
    TString title = r.spec.isAny ? "BH2(any seg)" : r.spec.tag;
    r.box.Print(title);
    TString cn1; cn1.Form("cHTOF_%s",r.spec.key.Data());
    r.c1=new TCanvas(cn1,title,800,600); r.box.hHTOF->Draw("hist");
    TString cn2; cn2.Form("cPairs_%s",r.spec.key.Data());
    r.c2=new TCanvas(cn2,TString("Adjacent pairs | ")+title,900,600);
    r.box.hPairAdj->LabelsOption("v","X"); r.box.hPairAdj->Draw("hist");
    if(save){
      TString f1; f1.Form("HTOF_pattern_%s.png",r.spec.key.Data());
      TString f2; f2.Form("HTOF_pairs_%s.png",  r.spec.key.Data());
      r.c1->SaveAs(f1); r.c2->SaveAs(f2);
    }
  }
}

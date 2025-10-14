// 2pi_study_v2.C
// 파일에 origin/pdg 브랜치가 없어도 PRM/SEC/TGT로 분석하는 버전
// 사용법:
//   root -l
//   .L 2pi_study_v2.C+
//   two_pi_study_v2("../rootfile/E45_fix_piplusn_098.root","g4hyptpc",4,10,true);

#include "TFile.h"
#include "TTree.h"
#include "TParticle.h"
#include "TSystem.h"
#include "TInterpreter.h"
#include "TCanvas.h"
#include "TH1I.h"
#include "TH2I.h"
#include <vector>
#include <set>
#include <map>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>

namespace TPI2 {
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

static inline double MIPThrMeV(double frac,double dEdx,double thk_mm){
  return frac * dEdx * (thk_mm*0.1); // mm→cm
}

// x(world) → BH2 seg id
static int MapBH2_WorldToSeg(double x, double, double){
  const double xloc = x - kBH2_x0;
  const double xmin = -0.5*kBH2_total, xmax = 0.5*kBH2_total;
  if (xloc < xmin || xloc >= xmax) return -1;
  int idx = int(std::floor((xloc - xmin)/kBH2_pitch));
  if(idx<0) idx=0; if(idx>=kNBH2Seg) idx=kNBH2Seg-1;
  return idx;
}

static inline double pct(long long a,long long b){ return (b>0)? 100.0*double(a)/double(b) : 0.0; }
} // ns

void two_pi_study_v2(const char* filename,
                     const char* treename="g4hyptpc",
                     int bh2_lo=4, int bh2_hi=10,
                     bool savePlots=true,
                     double mipFrac=0.1,
                     double mipMeVperCm=2.0,
                     double BH2_thickness_mm=5.0,
                     double HTOF_thickness_mm=10.0,
                     double TGT_anyEdep_thr=0.0) // TGT 히트판정(합edep ≥ thr). 기본 0: 존재만으로 인정
{
  using namespace TPI2;

  TFile* f=TFile::Open(filename,"READ");
  if(!f||f->IsZombie()){ std::cerr<<"[ERR] open "<<filename<<" failed\n"; return; }
  TTree* T=(TTree*)f->Get(treename);
  if(!T){ std::cerr<<"[ERR] tree "<<treename<<" not found\n"; return; }

  // 필요한 브랜치 체크
  auto need = [&](const char* b)->bool{ return T->GetBranch(b)!=nullptr; };
  if(!need("BH2") || !need("HTOF") || !need("PRM") || !need("SEC") || !need("TGT")){
    std::cerr<<"[ERR] need BH2, HTOF, PRM, SEC, TGT (vector<TParticle>)\n";
    return;
  }
  const bool hasBH2edep  = (T->GetBranch("BH2_edep")  != nullptr);
  const bool hasHTOFedep = (T->GetBranch("HTOF_edep") != nullptr);
  const bool hasTGTedep  = (T->GetBranch("TGT_edep")  != nullptr);

  // 브랜치 포인터
  std::vector<TParticle>* BH2=nullptr;   T->SetBranchAddress("BH2",&BH2);
  std::vector<TParticle>* HTOF=nullptr;  T->SetBranchAddress("HTOF",&HTOF);
  std::vector<TParticle>* PRM=nullptr;   T->SetBranchAddress("PRM",&PRM);
  std::vector<TParticle>* SECv=nullptr;  T->SetBranchAddress("SEC",&SECv);
  std::vector<TParticle>* TGT=nullptr;   T->SetBranchAddress("TGT",&TGT);
  std::vector<double>* BH2_edep=nullptr; if(hasBH2edep)  T->SetBranchAddress("BH2_edep",&BH2_edep);
  std::vector<double>* HTOF_edep=nullptr;if(hasHTOFedep) T->SetBranchAddress("HTOF_edep",&HTOF_edep);
  std::vector<double>* TGT_edep=nullptr; if(hasTGTedep)  T->SetBranchAddress("TGT_edep",&TGT_edep);

  const double thrBH2  = MIPThrMeV(mipFrac,mipMeVperCm,BH2_thickness_mm);   // ex) 0.10 MeV
  const double thrHTOF = MIPThrMeV(mipFrac,mipMeVperCm,HTOF_thickness_mm);  // ex) 0.20 MeV

  std::cout<<std::fixed<<std::setprecision(3);
  std::cout<<"[INFO] BH2 thr="<<thrBH2<<" MeV, HTOF thr="<<thrHTOF
           <<" MeV | BH2 range = "<<bh2_lo<<"-"<<bh2_hi<<"\n";
  std::cout<<"[INFO] TGT-hit policy: sum(edep) >= "<<TGT_anyEdep_thr<<" MeV (0 → presence-only)\n";

  // 카운터
  Long64_t N_total=0;
  Long64_t N_BH2_4_10_all=0;

  Long64_t N_noReaction=0; // TGT-hit 없음
  Long64_t N_BH2_4_10_in_noReaction=0;
  Long64_t N_BH2_4_10_reactionOnly=0;
  Long64_t N_trig1_in_noReaction=0;

  // (요청 4번) 히스토그램: "no-reaction에서 원빔 제외"라고 적혀있지만, 그 집합은 반응입자가 없으므로
  // 실제 유의미한 것은 "reaction-only(=TGT-hit 있음) & 원빔 제외"가 맞음.
  // 두 가지를 모두 만들어서 비교 가능하게 제공.
  TH1I *hNpiMinus_noReact=nullptr, *hNpiPlus_noReact=nullptr;
  TH2I *h2_pim_vs_piplus_noReact=nullptr, *h2_pim_vs_p_noReact=nullptr;

  TH1I *hNpiMinus_react=nullptr, *hNpiPlus_react=nullptr;
  TH2I *h2_pim_vs_piplus_react=nullptr, *h2_pim_vs_p_react=nullptr;

  auto makeH = [&](const char* suf){
    auto h1a=new TH1I(TString::Format("hNpiMinus_%s",suf),
                      TString::Format("N(pi-) per event | %s & no-beam;N(pi-);Events",suf),20,-0.5,19.5);
    auto h1b=new TH1I(TString::Format("hNpiPlus_%s",suf),
                      TString::Format("N(pi+) per event | %s & no-beam;N(pi+);Events",suf),20,-0.5,19.5);
    auto h2a=new TH2I(TString::Format("h2_Npim_vs_Piplus_%s",suf),
                      TString::Format("N(pi-) vs N(pi+) | %s & no-beam;N(pi+);N(pi-)",suf),
                      20,-0.5,19.5,20,-0.5,19.5);
    auto h2b=new TH2I(TString::Format("h2_Npim_vs_P_%s",suf),
                      TString::Format("N(pi-) vs N(p) | %s & no-beam;N(p);N(pi-)",suf),
                      20,-0.5,19.5,20,-0.5,19.5);
    h1a->SetDirectory(nullptr); h1b->SetDirectory(nullptr);
    h2a->SetDirectory(nullptr); h2b->SetDirectory(nullptr);
    return std::make_tuple(h1a,h1b,h2a,h2b);
  };

  std::tie(hNpiMinus_noReact,hNpiPlus_noReact,h2_pim_vs_piplus_noReact,h2_pim_vs_p_noReact) = makeH("noReaction");
  std::tie(hNpiMinus_react,  hNpiPlus_react,  h2_pim_vs_piplus_react,  h2_pim_vs_p_react)   = makeH("reaction");

  const Long64_t N=T->GetEntries();
  for(Long64_t ie=0; ie<N; ++ie){
    T->GetEntry(ie); N_total++;

    // --- BH2 유효 세그
    std::map<int,double> bh2E;
    for(size_t i=0; BH2 && i<BH2->size(); ++i){
      const TParticle& p=BH2->at(i);
      int sid = MapBH2_WorldToSeg(p.Vx(),p.Vy(),p.Vz());
      if(0<=sid && sid<kNBH2Seg){
        const double ed = (hasBH2edep && BH2_edep && i<BH2_edep->size()) ? BH2_edep->at(i) : p.GetWeight();
        bh2E[sid]+=ed;
      }
    }
    std::set<int> bh2Valid; for(auto& kv:bh2E) if(kv.second>=thrBH2) bh2Valid.insert(kv.first);
    const bool bh2_in_4_10 = std::any_of(bh2Valid.begin(), bh2Valid.end(),
                                         [&](int s){ return (bh2_lo<=s && s<=bh2_hi); });
    if(bh2_in_4_10) N_BH2_4_10_all++;

    // --- HTOF 유효 타일
    std::map<int,double> htofE;
    for(size_t i=0; HTOF && i<HTOF->size(); ++i){
      const TParticle& p=HTOF->at(i);
      int tid = p.GetStatusCode();
      if(0<=tid && tid<kNHTOF){
        const double ed = (hasHTOFedep && HTOF_edep && i<HTOF_edep->size()) ? HTOF_edep->at(i) : p.GetWeight();
        htofE[tid]+=ed;
      }
    }
    std::set<int> htofValid; for(auto& kv:htofE) if(kv.second>=thrHTOF) htofValid.insert(kv.first);
    const int htofMult = (int)htofValid.size();

    // --- TGT 히트 존재 여부(합edep >= thr). thr=0이면 존재만으로 true
    double tgtE=0.0;
    for(size_t i=0; TGT && i<TGT->size(); ++i){
      const double ed = (hasTGTedep && TGT_edep && i<TGT_edep->size()) ? TGT_edep->at(i) : TGT->at(i).GetWeight();
      tgtE += ed;
    }
    const bool hasTGT = (TGT && (!TGT->empty())) && (tgtE >= TGT_anyEdep_thr);

    // --- 반응입자 카운트(원빔 제외 → PRM의 pi- 제외, SEC에서만 카운트)
    auto fillCounts = [](const std::vector<TParticle>* V, int& nPim, int& nPip, int& nP){
      for(size_t i=0; V && i<V->size(); ++i){
        const int pdg = V->at(i).GetPdgCode();
        if(pdg==-211) nPim++;
        else if(pdg==211) nPip++;
        else if(pdg==2212) nP++;
      }
    };

    // (2) no-reaction = !hasTGT
    if(!hasTGT){
      N_noReaction++;
      if(bh2_in_4_10) N_BH2_4_10_in_noReaction++;
      if(bh2_in_4_10 && htofMult>=2) N_trig1_in_noReaction++;

      // 요청 4번 원문대로면 여기서 "reaction만"을 그리라 했지만, no-reaction 집합엔 반응입자가 없음.
      // 그래도 확인용으로 SEC 카운트를 그려보면 대부분 0일 것임.
      int nPim=0,nPip=0,nP=0;
      fillCounts(SECv,nPim,nPip,nP); // SEC만(원빔제외 효과)
      hNpiMinus_noReact->Fill(nPim);
      hNpiPlus_noReact ->Fill(nPip);
      h2_pim_vs_piplus_noReact->Fill(nPip,nPim);
      h2_pim_vs_p_noReact     ->Fill(nP,nPim);
    } else {
      // (참고) reaction-only에서의 분포(아마 이게 네가 원하는 그림)
      int nPim=0,nPip=0,nP=0;
      fillCounts(SECv,nPim,nPip,nP); // SEC만 사용 → 원빔(PRM의 pi-) 제외
      hNpiMinus_react->Fill(nPim);
      hNpiPlus_react ->Fill(nPip);
      h2_pim_vs_piplus_react->Fill(nPip,nPim);
      h2_pim_vs_p_react     ->Fill(nP,nPim);
    }
  } // loop

  N_BH2_4_10_reactionOnly = N_BH2_4_10_all - N_BH2_4_10_in_noReaction;

  // 요약 출력
  std::cout<<"\n===== 2pi_study_v2 summary (BH2 "<<bh2_lo<<"-"<<bh2_hi<<") =====\n";
  std::cout<<"Total events                                 : "<<N_total<<"\n";
  std::cout<<"BH2[4-10] hit events (all)                   : "<<N_BH2_4_10_all
           <<"  ("<<pct(N_BH2_4_10_all,N_total)<<" % of Total)\n";
  std::cout<<"No-reaction (no TGT hit)                     : "<<N_noReaction
           <<"  ("<<pct(N_noReaction,N_total)<<" % of Total)\n";
  std::cout<<"BH2[4-10] after excluding no-reaction        : "<<N_BH2_4_10_reactionOnly
           <<"  ("<<pct(N_BH2_4_10_reactionOnly, N_total - N_noReaction)<<" % of [Total - NoReaction])\n";
  std::cout<<"Trig1 in no-reaction (BH2[4-10] & m>=2)      : "<<N_trig1_in_noReaction
           <<"  ("<<pct(N_trig1_in_noReaction, N_noReaction)<<" % of NoReaction)\n";

  // 그림 저장
  TCanvas* c1=new TCanvas("cNpi_noReact","N(pi), no-reaction set",1000,500);
  c1->Divide(2,1);
  c1->cd(1); hNpiMinus_noReact->Draw("hist");
  c1->cd(2); hNpiPlus_noReact ->Draw("hist");

  TCanvas* c2=new TCanvas("cPairs_noReact","pairs, no-reaction set",1000,500);
  c2->Divide(2,1);
  c2->cd(1); h2_pim_vs_piplus_noReact->Draw("colz");
  c2->cd(2); h2_pim_vs_p_noReact     ->Draw("colz");

  TCanvas* c3=new TCanvas("cNpi_react","N(pi), reaction-only",1000,500);
  c3->Divide(2,1);
  c3->cd(1); hNpiMinus_react->Draw("hist");
  c3->cd(2); hNpiPlus_react ->Draw("hist");

  TCanvas* c4=new TCanvas("cPairs_react","pairs, reaction-only",1000,500);
  c4->Divide(2,1);
  c4->cd(1); h2_pim_vs_piplus_react->Draw("colz");
  c4->cd(2); h2_pim_vs_p_react     ->Draw("colz");

  if(savePlots){
    c1->SaveAs("2pi_noReaction_Npi.png");
    c2->SaveAs("2pi_noReaction_pairs.png");
    c3->SaveAs("2pi_reaction_Npi.png");
    c4->SaveAs("2pi_reaction_pairs.png");
  }
}

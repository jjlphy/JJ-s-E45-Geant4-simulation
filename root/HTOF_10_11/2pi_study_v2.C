// 2pi_study_v2.C
// - origin/pdg 전용 브랜치 없이도 PRM/SEC/TGT로 분석
// - 강제 2π 플래그/채널(tgt_touch_flag, forced2pi_flag, forced2pi_channel) 활용
// - vector<TParticle> 딕셔너리를 함수 시작 시 보장(맥/ACLiC 안전)
//
// 사용:
//   root -l
//   .L 2pi_study_v2.C+
//   two_pi_study_v2("../rootfile/E45_with_SCH.root","g4hyptpc",4,10,true);

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TParticle.h"
#include "TSystem.h"
#include "TInterpreter.h"
#include "TCanvas.h"
#include "TH1I.h"
#include "TH2I.h"
#include "TString.h"

#include <vector>
#include <set>
#include <map>
#include <tuple>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>

// ----------------------------------------------------------------------------
// 딕셔너리 보장 (반드시 SetBranchAddress 전에 호출)
// ----------------------------------------------------------------------------
static inline void EnsureVecTParticleDict(){
  gSystem->Load("libPhysics");
  gInterpreter->GenerateDictionary("vector<TParticle>","TParticle.h;vector");
}

// ----------------------------------------------------------------------------
// 기하/임계값 유틸
// ----------------------------------------------------------------------------
namespace TPI2 {
static const int    kNBH2Seg   = 15;     // seg id: 0..14
static const double kBH2_x0    = -10.0;  // [mm] BH2 원점 보정(필요시 조정)
static const double kBH2_pitch = 14.0;   // [mm]
static const double kBH2_total = kNBH2Seg * kBH2_pitch;
static const int    kNHTOF     = 34;     // tile id: 0..33

static inline double MIPThrMeV(double frac,double dEdx,double thk_mm){
  return frac * dEdx * (thk_mm*0.1); // mm→cm
}

// x(world) → BH2 seg id (단순 1D 매핑; 필요시 보정)
static inline int MapBH2_WorldToSeg(double x, double y, double z){
  (void)y; (void)z;
  const double xloc = x - kBH2_x0;
  const double xmin = -0.5*kBH2_total, xmax = 0.5*kBH2_total;
  if (xloc < xmin || xloc >= xmax) return -1;
  int idx = int(std::floor((xloc - xmin)/kBH2_pitch));
  if(idx<0) idx=0; if(idx>=kNBH2Seg) idx=kNBH2Seg-1;
  return idx;
}

static inline double pct(long long a,long long b){ return (b>0)? 100.0*double(a)/double(b) : 0.0; }
} // ns

// ----------------------------------------------------------------------------
// 메인
// ----------------------------------------------------------------------------
void two_pi_study_v2(const char* filename,
                     const char* treename="g4hyptpc",
                     int bh2_lo=4, int bh2_hi=10,
                     bool savePlots=true,
                     // 임계값 파라미터(필요시 조정)
                     double mipFrac=0.1,
                     double mipMeVperCm=2.0,
                     double BH2_thickness_mm=5.0,
                     double HTOF_thickness_mm=10.0,
                     // TGT 판정 정책: <0 → tgt_touch_flag 사용, >=0 → TGT_edep 합≥thr
                     double TGT_anyEdep_thr=-1.0)
{
  using namespace TPI2;

  EnsureVecTParticleDict(); // ★ 반드시 가장 먼저

  // 파일/트리
  TFile* f=TFile::Open(filename,"READ");
  if(!f||f->IsZombie()){ std::cerr<<"[ERR] open "<<filename<<" failed\n"; return; }
  TTree* T=(TTree*)f->Get(treename);
  if(!T){ std::cerr<<"[ERR] tree "<<treename<<" not found\n"; return; }

  // 필요한 컬렉션(필수)
  auto need = [&](const char* b)->bool{ return T->GetBranch(b)!=nullptr; };
  if(!need("BH2") || !need("HTOF") || !need("PRM") || !need("SEC") || !need("TGT")){
    std::cerr<<"[ERR] need branches: BH2, HTOF, PRM, SEC, TGT (vector<TParticle>)\n";
    return;
  }

  // 에너지 브랜치 존재 여부
  const bool hasBH2edep  = (T->GetBranch("BH2_edep")  != nullptr);
  const bool hasHTOFedep = (T->GetBranch("HTOF_edep") != nullptr);
  const bool hasTGTedep  = (T->GetBranch("TGT_edep")  != nullptr);

  // 새 truth 플래그 (없으면 graceful fallback)
  const bool hasForced2piFlag    = (T->GetBranch("forced2pi_flag")    != nullptr);
  const bool hasForced2piChannel = (T->GetBranch("forced2pi_channel") != nullptr);
  const bool hasTgtTouchFlag     = (T->GetBranch("tgt_touch_flag")    != nullptr);

  // 브랜치 주소
  std::vector<TParticle>* BH2=nullptr;   T->SetBranchAddress("BH2",&BH2);
  std::vector<TParticle>* HTOF=nullptr;  T->SetBranchAddress("HTOF",&HTOF);
  std::vector<TParticle>* PRM=nullptr;   T->SetBranchAddress("PRM",&PRM);
  std::vector<TParticle>* SEC=nullptr;   T->SetBranchAddress("SEC",&SEC);
  std::vector<TParticle>* TGT=nullptr;   T->SetBranchAddress("TGT",&TGT);

  std::vector<double>* BH2_edep=nullptr; if(hasBH2edep)  T->SetBranchAddress("BH2_edep",&BH2_edep);
  std::vector<double>* HTOF_edep=nullptr;if(hasHTOFedep) T->SetBranchAddress("HTOF_edep",&HTOF_edep);
  std::vector<double>* TGT_edp =nullptr; if(hasTGTedep)  T->SetBranchAddress("TGT_edep",&TGT_edp);

  int forced2pi_flag=0, forced2pi_channel=-1, tgt_touch_flag=0;
  if(hasForced2piFlag)    T->SetBranchAddress("forced2pi_flag",&forced2pi_flag);
  if(hasForced2piChannel) T->SetBranchAddress("forced2pi_channel",&forced2pi_channel);
  if(hasTgtTouchFlag)     T->SetBranchAddress("tgt_touch_flag",&tgt_touch_flag);

  // 딕셔너리/바인딩 sanity check
  T->GetEntry(0);
  if(!SEC){ std::cerr<<"[FATAL] SEC not attached (vector<TParticle> dict missing?)\n"; return; }

  const double thrBH2  = MIPThrMeV(mipFrac,mipMeVperCm,BH2_thickness_mm);   // ex) 0.10 MeV
  const double thrHTOF = MIPThrMeV(mipFrac,mipMeVperCm,HTOF_thickness_mm);  // ex) 0.20 MeV

  std::cout<<std::fixed<<std::setprecision(3);
  std::cout<<"[INFO] BH2 thr="<<thrBH2<<" MeV, HTOF thr="<<thrHTOF
           <<" MeV | BH2 range = "<<bh2_lo<<"-"<<bh2_hi<<"\n";
  if(hasTgtTouchFlag && TGT_anyEdep_thr<0)
    std::cout<<"[INFO] TGT-hit policy: use tgt_touch_flag(1=hit)\n";
  else
    std::cout<<"[INFO] TGT-hit policy: sum(TGT edep) >= "<<std::max(0.0,TGT_anyEdep_thr)<<" MeV\n";

  // 카운터
  Long64_t N_total=0, N_BH2_4_10_all=0;
  Long64_t N_noReaction=0, N_BH2_4_10_in_noReaction=0, N_BH2_4_10_reactionOnly=0, N_trig1_in_noReaction=0;

  Long64_t N_ch_total[2]={0,0}, N_ch_bh2_in[2]={0,0}; // per-channel

  // 히스토그램
  auto makeH = [](const char* suf){
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

  TH1I *hNpiMinus_noReact,*hNpiPlus_noReact; TH2I *h2_pim_vs_piplus_noReact,*h2_pim_vs_p_noReact;
  TH1I *hNpiMinus_react,*hNpiPlus_react;     TH2I *h2_pim_vs_piplus_react,*h2_pim_vs_p_react;

  std::tie(hNpiMinus_noReact,hNpiPlus_noReact,h2_pim_vs_piplus_noReact,h2_pim_vs_p_noReact) = makeH("noReaction");
  std::tie(hNpiMinus_react,  hNpiPlus_react,  h2_pim_vs_piplus_react,  h2_pim_vs_p_react)   = makeH("reaction");

  TH1I *hNpiMinus_react_ch0,*hNpiPlus_react_ch0; TH2I *h2_pim_vs_piplus_react_ch0,*h2_pim_vs_p_react_ch0;
  TH1I *hNpiMinus_react_ch1,*hNpiPlus_react_ch1; TH2I *h2_pim_vs_piplus_react_ch1,*h2_pim_vs_p_react_ch1;

  std::tie(hNpiMinus_react_ch0,hNpiPlus_react_ch0,h2_pim_vs_piplus_react_ch0,h2_pim_vs_p_react_ch0) = makeH("reaction_ch0");
  std::tie(hNpiMinus_react_ch1,hNpiPlus_react_ch1,h2_pim_vs_piplus_react_ch1,h2_pim_vs_p_react_ch1) = makeH("reaction_ch1");

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

    // --- HTOF 유효 타일(multiplicity m)
    std::map<int,double> htofE;
    for(size_t i=0; HTOF && i<HTOF->size(); ++i){
      const TParticle& p=HTOF->at(i);
      int tid = p.GetStatusCode(); // 타일 id를 StatusCode에 저장해둠
      if(0<=tid && tid<kNHTOF){
        const double ed = (hasHTOFedep && HTOF_edep && i<HTOF_edep->size()) ? HTOF_edep->at(i) : p.GetWeight();
        htofE[tid]+=ed;
      }
    }
    std::set<int> htofValid; for(auto& kv:htofE) if(kv.second>=thrHTOF) htofValid.insert(kv.first);
    const int htofMult = (int)htofValid.size();

    // --- TGT hit 판정 (flag 우선, 없으면 edep 합)
    bool hasTGT=false;
    if(hasTgtTouchFlag && TGT_anyEdep_thr<0){
      hasTGT = (tgt_touch_flag!=0);
    }else{
      double tgtE=0.0;
      if(hasTGTedep && TGT_edp){
        for(size_t i=0; TGT && i<TGT->size() && i<TGT_edp->size(); ++i) tgtE += TGT_edp->at(i);
      }else{
        // edep 브랜치가 없으면 weight 합(0이면 존재만으로 판정)
        for(size_t i=0; TGT && i<TGT->size(); ++i) tgtE += std::max(0.0, (double)TGT->at(i).GetWeight());
      }
      const double thr = std::max(0.0, TGT_anyEdep_thr);
      hasTGT = (TGT && (!TGT->empty())) && (tgtE >= thr);
    }

    auto countSEC = [](const std::vector<TParticle>* V, int& nPim, int& nPip, int& nP){
      nPim=nPip=nP=0;
      for(size_t i=0; V && i<V->size(); ++i){
        const int pdg = V->at(i).GetPdgCode();
        if(pdg==-211) nPim++;
        else if(pdg==211) nPip++;
        else if(pdg==2212) nP++;
      }
    };

    if(!hasTGT){
      N_noReaction++;
      if(bh2_in_4_10) N_BH2_4_10_in_noReaction++;
      if(bh2_in_4_10 && htofMult>=2) N_trig1_in_noReaction++;

      int nPim=0,nPip=0,nP=0;
      countSEC(SEC,nPim,nPip,nP); // SEC만 사용 → 원빔(PRM) 제외
      hNpiMinus_noReact->Fill(nPim);
      hNpiPlus_noReact ->Fill(nPip);
      h2_pim_vs_piplus_noReact->Fill(nPip,nPim);
      h2_pim_vs_p_noReact     ->Fill(nP,nPim);
    } else {
      int nPim=0,nPip=0,nP=0;
      countSEC(SEC,nPim,nPip,nP);
      hNpiMinus_react->Fill(nPim);
      hNpiPlus_react ->Fill(nPip);
      h2_pim_vs_piplus_react->Fill(nPip,nPim);
      h2_pim_vs_p_react     ->Fill(nP,nPim);

      // per-channel (강제2π만 집계)
      if(hasForced2piFlag && hasForced2piChannel && forced2pi_flag==1 &&
         (forced2pi_channel==0 || forced2pi_channel==1))
      {
        const int ch = forced2pi_channel;
        N_ch_total[ch]++; if(bh2_in_4_10) N_ch_bh2_in[ch]++;
        if(ch==0){
          hNpiMinus_react_ch0->Fill(nPim); hNpiPlus_react_ch0->Fill(nPip);
          h2_pim_vs_piplus_react_ch0->Fill(nPip,nPim);
          h2_pim_vs_p_react_ch0     ->Fill(nP,nPim);
        }else{
          hNpiMinus_react_ch1->Fill(nPim); hNpiPlus_react_ch1->Fill(nPip);
          h2_pim_vs_piplus_react_ch1->Fill(nPip,nPim);
          h2_pim_vs_p_react_ch1     ->Fill(nP,nPim);
        }
      }
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

  if(T->GetBranch("forced2pi_flag") && T->GetBranch("forced2pi_channel")){
    std::cout<<"\n-- Per-channel (forced2pi_flag==1) --\n";
    for(int ch=0; ch<2; ++ch){
      std::cout<<"  ch"<<ch<<": total="<<N_ch_total[ch]
               <<",  BH2-in-range="<<N_ch_bh2_in[ch]
               <<"  ("<<pct(N_ch_bh2_in[ch], std::max<Long64_t>(1,N_ch_total[ch]))<<" %)\n";
    }
  }

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

  TCanvas* c5=new TCanvas("cNpi_react_byCh","N(pi) | reaction-only by channel",1000,500);
  c5->Divide(2,1);
  c5->cd(1); hNpiMinus_react_ch0->Draw("hist"); c5->cd(1)->SetGrid();
  c5->cd(2); hNpiPlus_react_ch0 ->Draw("hist"); c5->cd(2)->SetGrid();

  TCanvas* c6=new TCanvas("cPairs_react_byCh","pairs | reaction-only by channel",1000,500);
  c6->Divide(2,1);
  c6->cd(1); h2_pim_vs_piplus_react_ch0->Draw("colz");
  c6->cd(2); h2_pim_vs_p_react_ch1     ->Draw("colz"); // 의도적으로 ch1만 표시해 비교 가능

  if(savePlots){
    c1->SaveAs("2pi_noReaction_Npi.png");
    c2->SaveAs("2pi_noReaction_pairs.png");
    c3->SaveAs("2pi_reaction_Npi.png");
    c4->SaveAs("2pi_reaction_pairs.png");
    c5->SaveAs("2pi_reaction_Npi_byChannel.png");
    c6->SaveAs("2pi_reaction_pairs_byChannel.png");
  }
}

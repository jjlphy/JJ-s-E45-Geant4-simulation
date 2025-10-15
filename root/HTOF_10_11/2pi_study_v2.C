// 2pi_study_v2.C
// 파일에 origin/pdg 브랜치가 없어도 PRM/SEC/TGT로 분석하는 버전(+ HTOF 2D hit pattern by copy_no)
// 사용법:
//   root -l
//   .L 2pi_study_v2.C+
//   two_pi_study_v2("../rootfile/E45_with_SCH.root","g4hyptpc",4,10,true);

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

// ====== Geometry / mapping (E72-like) ======
static const int    kNBH2Seg   = 15;     // 0..14
static const double kBH2_x0    = -10.0;  // [mm]
static const double kBH2_pitch = 14.0;   // [mm]
static const double kBH2_total = kNBH2Seg * kBH2_pitch;
static const int    kNHTOF     = 34;     // 0..33

// ROOT dict for vector<TParticle>
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

// ============================================================================
// 메인 매크로
// ============================================================================
void two_pi_study_v2(const char* filename,
                     const char* treename="g4hyptpc",
                     int bh2_lo=4, int bh2_hi=10,
                     bool savePlots=true,
                     double mipFrac=0.1,
                     double mipMeVperCm=2.0,
                     double BH2_thickness_mm=5.0,
                     double HTOF_thickness_mm=10.0,
                     double TGT_anyEdep_thr=0.0) // 0: 존재만으로 true, >0: 합edep 컷
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

  const bool hasBH2edep   = need("BH2_edep");
  const bool hasHTOFedep  = need("HTOF_edep");
  const bool hasHTOFcopyN = need("HTOF_copyNo");         // ★ copy_no 사용 여부

  // 브랜치 포인터
  std::vector<TParticle>* BH2=nullptr;    T->SetBranchAddress("BH2",&BH2);
  std::vector<TParticle>* HTOF=nullptr;   T->SetBranchAddress("HTOF",&HTOF);
  std::vector<TParticle>* PRM=nullptr;    T->SetBranchAddress("PRM",&PRM);
  std::vector<TParticle>* SECv=nullptr;   T->SetBranchAddress("SEC",&SECv);
  std::vector<TParticle>* TGT=nullptr;    T->SetBranchAddress("TGT",&TGT);
  std::vector<double>* BH2_edep=nullptr;  if(hasBH2edep)   T->SetBranchAddress("BH2_edep",&BH2_edep);
  std::vector<double>* HTOF_edep=nullptr; if(hasHTOFedep)  T->SetBranchAddress("HTOF_edep",&HTOF_edep);
  std::vector<int>*    HTOF_copyNo=nullptr; if(hasHTOFcopyN) T->SetBranchAddress("HTOF_copyNo",&HTOF_copyNo);

  // event-level flags (있을 때만 사용)
  int forced2pi_flag=0, forced2pi_ch=-1, tgt_touch_flag=0;
  const bool hasF2Pflag = need("forced2pi_flag");
  const bool hasF2Pchan = need("forced2pi_channel");
  const bool hasTgtFlag = need("tgt_touch_flag");
  if(hasF2Pflag) T->SetBranchAddress("forced2pi_flag",&forced2pi_flag);
  if(hasF2Pchan) T->SetBranchAddress("forced2pi_channel",&forced2pi_ch);
  if(hasTgtFlag) T->SetBranchAddress("tgt_touch_flag",&tgt_touch_flag);

  // thresholds
  const double thrBH2  = MIPThrMeV(mipFrac,mipMeVperCm,BH2_thickness_mm);   // ex) 0.10 MeV
  const double thrHTOF = MIPThrMeV(mipFrac,mipMeVperCm,HTOF_thickness_mm);  // ex) 0.20 MeV

  std::cout<<std::fixed<<std::setprecision(3);
  std::cout<<"[INFO] BH2 thr="<<thrBH2<<" MeV, HTOF thr="<<thrHTOF
           <<" MeV | BH2 range = "<<bh2_lo<<"-"<<bh2_hi<<"\n";
  std::cout<<"[INFO] TGT-hit policy: "
           <<(TGT_anyEdep_thr>0.0? Form("sum(edep) >= %.3f MeV",TGT_anyEdep_thr)
                                 : "presence-only (>=1 hit)")<<"\n";

  // ---------------------- 통계 카운터 ----------------------
  Long64_t N_total=0;
  Long64_t N_BH2_4_10_all=0;

  Long64_t N_noReaction=0; // TGT-hit 없음
  Long64_t N_BH2_4_10_in_noReaction=0;
  Long64_t N_BH2_4_10_reactionOnly=0;
  Long64_t N_trig1_in_noReaction=0; // (BH2 in range & HTOF m>=2)

  // channel별 통계 (forced2pi_flag==1일 때만 집계)
  Long64_t ch_total[2]={0,0};
  Long64_t ch_bh2in [2]={0,0};

  // ---------------------- 히스토그램 ----------------------
  auto makeH1 = [&](const char* name,const char* title){
    auto h=new TH1I(name,title,20,-0.5,19.5); h->SetDirectory(nullptr); return h;
  };
  auto makeH2 = [&](const char* name,const char* title,
                    const char* xlab,const char* ylab){
    auto h=new TH2I(name,title,20,-0.5,19.5,20,-0.5,19.5);
    h->SetDirectory(nullptr); h->GetXaxis()->SetTitle(xlab);
    h->GetYaxis()->SetTitle(ylab); return h;
  };

  TH1I *hNpiMinus_noReact=makeH1("hNpiMinus_noReaction","N(pi-) per event | noReaction & no-beam;N(pi-);Events");
  TH1I *hNpiPlus_noReact =makeH1("hNpiPlus_noReaction","N(pi+) per event | noReaction & no-beam;N(pi+);Events");
  TH2I *h2_pim_vs_piplus_noReact = makeH2("h2_Npim_vs_Piplus_noReaction","N(pi-) vs N(pi+) | noReaction & no-beam","N(pi+)","N(pi-)");
  TH2I *h2_pim_vs_p_noReact      = makeH2("h2_Npim_vs_P_noReaction","N(pi-) vs N(p) | noReaction & no-beam","N(p)","N(pi-)");

  TH1I *hNpiMinus_react=makeH1("hNpiMinus_reaction","N(pi-) per event | reaction & no-beam;N(pi-);Events");
  TH1I *hNpiPlus_react =makeH1("hNpiPlus_reaction","N(pi+) per event | reaction & no-beam;N(pi+);Events");
  TH2I *h2_pim_vs_piplus_react = makeH2("h2_Npim_vs_Piplus_reaction","N(pi-) vs N(pi+) | reaction & no-beam","N(pi+)","N(pi-)");
  TH2I *h2_pim_vs_p_react      = makeH2("h2_Npim_vs_P_reaction","N(pi-) vs N(p) | reaction & no-beam","N(p)","N(pi-)");

  // per-channel(0/1) 반응 전용
  TH1I *hNpiMinus_react_ch0=makeH1("hNpiMinus_reaction_ch0","N(pi-) per event | reaction_ch0 & no-beam;N(pi-);Events");
  TH1I *hNpiPlus_react_ch0 =makeH1("hNpiPlus_reaction_ch0","N(pi+) per event | reaction_ch0 & no-beam;N(pi+);Events");
  TH1I *hNpiMinus_react_ch1=makeH1("hNpiMinus_reaction_ch1","N(pi-) per event | reaction_ch1 & no-beam;N(pi-);Events");
  TH1I *hNpiPlus_react_ch1 =makeH1("hNpiPlus_reaction_ch1","N(pi+) per event | reaction_ch1 & no-beam;N(pi+);Events");
  TH2I *h2_pim_vs_piplus_react_ch0 = makeH2("h2_Npim_vs_Piplus_reaction_ch0","N(pi-) vs N(pi+) | reaction_ch0 & no-beam","N(pi+)","N(pi-)");
  TH2I *h2_pim_vs_p_react_ch1      = makeH2("h2_Npim_vs_P_reaction_ch1","N(pi-) vs N(p) | reaction_ch1 & no-beam","N(p)","N(pi-)");

  // ===== HTOF 2D hit pattern by copy_no (0..33 × 0..33) =====
  TH2I* h2_HTOF_hitpattern = new TH2I("h2_HTOF_hitpattern",
    "HTOF 2D hit pattern (tileID by copy_no);HTOF tile i (0..33);HTOF tile j (0..33)",
    34,-0.5,33.5, 34,-0.5,33.5);
  h2_HTOF_hitpattern->SetDirectory(nullptr);
  const bool includeDiagonal = false; // i==j 포함 여부

  // ---------------------- 이벤트 루프 ----------------------
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

    // --- HTOF 유효 타일 집합 (copy_no 우선, 없으면 statusCode)
    std::set<int> htofValid;
    if(hasHTOFcopyN && HTOF_copyNo && HTOF && HTOF_copyNo->size()==HTOF->size()){
      // copy_no별 에너지 합 → thrHTOF 이상이면 유효 타일
      std::map<int,double> eByCopyNo;
      for(size_t i=0;i<HTOF->size();++i){
        const int cn = HTOF_copyNo->at(i);
        if(cn<0 || cn>=kNHTOF) continue;
        const double ed = (hasHTOFedep && HTOF_edep && i<HTOF_edep->size()) ? HTOF_edep->at(i)
                                                                             : HTOF->at(i).GetWeight();
        eByCopyNo[cn] += ed;
      }
      for(auto& kv: eByCopyNo) if(kv.second>=thrHTOF) htofValid.insert(kv.first);
    } else {
      // 백업: statusCode별 에너지 합 → thrHTOF 이상이면 유효 타일
      std::map<int,double> eByStatus;
      for(size_t i=0; HTOF && i<HTOF->size(); ++i){
        const int tid = HTOF->at(i).GetStatusCode();
        if(tid<0 || tid>=kNHTOF) continue;
        const double ed = (hasHTOFedep && HTOF_edep && i<HTOF_edep->size()) ? HTOF_edep->at(i)
                                                                             : HTOF->at(i).GetWeight();
        eByStatus[tid] += ed;
      }
      for(auto& kv: eByStatus) if(kv.second>=thrHTOF) htofValid.insert(kv.first);
    }
    const int htofMult = (int)htofValid.size();

    // --- TGT 히트 존재 여부(합edep >= thr). thr=0이면 존재만으로 true
    bool hasTGT=false;
    if(TGT_anyEdep_thr<=0.0){
      hasTGT = (TGT && !TGT->empty());
    } else {
      double tgtE=0.0;
      for(size_t i=0; TGT && i<TGT->size(); ++i){
        const double ed = TGT->at(i).GetWeight(); // TGT_edep 없다고 가정(있으면 바꿔도 됨)
        tgtE += ed;
      }
      hasTGT = (TGT && (!TGT->empty())) && (tgtE >= TGT_anyEdep_thr);
    }

    // --- 반응입자 카운트(원빔 제외 → SEC만 카운트)
    auto countSEC = [](const std::vector<TParticle>* V, int& nPim, int& nPip, int& nP){
      nPim=nPip=nP=0;
      for(size_t i=0; V && i<V->size(); ++i){
        const int pdg = V->at(i).GetPdgCode();
        if(pdg==-211) nPim++;
        else if(pdg== 211) nPip++;
        else if(pdg==2212) nP++;
      }
    };

    // --- HTOF 2D hit pattern (by copy_no) : 모든 순서쌍(i,j) 채우기
    if(!htofValid.empty()){
      std::vector<int> tiles(htofValid.begin(), htofValid.end());
      for(size_t a=0; a<tiles.size(); ++a){
        for(size_t b=0; b<tiles.size(); ++b){
          if(!includeDiagonal && a==b) continue;
          h2_HTOF_hitpattern->Fill(tiles[a], tiles[b]);
        }
      }
    }

    // --- BH2 범위 조건 만족시 통계
    if(!bh2_in_4_10) continue;

    if(!hasTGT){
      N_noReaction++;
      if(bh2_in_4_10) N_BH2_4_10_in_noReaction++;
      if(bh2_in_4_10 && htofMult>=2) N_trig1_in_noReaction++;

      int nPim=0,nPip=0,nP=0;
      countSEC(SECv,nPim,nPip,nP);
      hNpiMinus_noReact->Fill(nPim);
      hNpiPlus_noReact ->Fill(nPip);
      h2_pim_vs_piplus_noReact->Fill(nPip,nPim);
      h2_pim_vs_p_noReact     ->Fill(nP,nPim);

    } else {
      // reaction-only 분포
      int nPim=0,nPip=0,nP=0;
      countSEC(SECv,nPim,nPip,nP);
      hNpiMinus_react->Fill(nPim);
      hNpiPlus_react ->Fill(nPip);
      h2_pim_vs_piplus_react->Fill(nPip,nPim);
      h2_pim_vs_p_react     ->Fill(nP,nPim);

      // per-channel 집계(강제 2pi 플래그가 있을 때만)
      if(hasF2Pflag && hasF2Pchan && forced2pi_flag==1 && (forced2pi_ch==0 || forced2pi_ch==1)){
        ch_total[forced2pi_ch]++; if(bh2_in_4_10) ch_bh2in[forced2pi_ch]++;
        if(forced2pi_ch==0){
          hNpiMinus_react_ch0->Fill(nPim);
          hNpiPlus_react_ch0 ->Fill(nPip);
          h2_pim_vs_piplus_react_ch0->Fill(nPip,nPim);
        }else{
          hNpiMinus_react_ch1->Fill(nPim);
          hNpiPlus_react_ch1 ->Fill(nPip);
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

  if(hasF2Pflag && hasF2Pchan){
    std::cout<<"\n-- Per-channel (forced2pi_flag==1) --\n";
    std::cout<<"  ch0: total="<<ch_total[0]<<",  BH2-in-range="<<ch_bh2in[0]
             <<"  ("<<pct(ch_bh2in[0], (ch_total[0]? ch_total[0]:1))<<" %)\n";
    std::cout<<"  ch1: total="<<ch_total[1]<<",  BH2-in-range="<<ch_bh2in[1]
             <<"  ("<<pct(ch_bh2in[1], (ch_total[1]? ch_total[1]:1))<<" %)\n";
  }

  // --------- 그림 저장 ---------
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

  // per-channel
  TCanvas* c5=new TCanvas("cNpi_react_byChannel","N(pi) reaction-only by channel",1000,500);
  c5->Divide(2,1);
  c5->cd(1); hNpiMinus_react_ch0->Draw("hist");
  c5->cd(2); hNpiPlus_react_ch0 ->Draw("hist");

  TCanvas* c6=new TCanvas("cPairs_react_byChannel","pairs reaction-only by channel",1000,500);
  c6->Divide(2,1);
  c6->cd(1); h2_pim_vs_piplus_react_ch0->Draw("colz");
  c6->cd(2); h2_pim_vs_p_react_ch1     ->Draw("colz");

  // HTOF 2D hit pattern
  TCanvas* cHTOFmap = new TCanvas("cHTOFmap","HTOF 2D hit pattern",650,600);
  h2_HTOF_hitpattern->Draw("colz");

  if(savePlots){
    c1->SaveAs("2pi_noReaction_Npi.png");
    c2->SaveAs("2pi_noReaction_pairs.png");
    c3->SaveAs("2pi_reaction_Npi.png");
    c4->SaveAs("2pi_reaction_pairs.png");
    c5->SaveAs("2pi_reaction_Npi_byChannel.png");
    c6->SaveAs("2pi_reaction_pairs_byChannel.png");
    cHTOFmap->SaveAs("2pi_HTOF_hitpattern.png");  // ★ copy_no 기반 2D 맵
  }
}

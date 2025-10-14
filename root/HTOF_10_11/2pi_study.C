// 2pi_study.C
// E45/E72 forced-2pi 분석 유틸 (BH2・HTOF 정책은 HTOF_* 참고본과 동일)
//
// 사용 예:
//   root -l
//   .L 2pi_study.C+
//   two_pi_study("../rootfile/E45_fix_Beam_098.root", "g4hyptpc", 4, 10, true);
//
// 출력:
//   - 터미널: 총 이벤트, BH2[4-10], no-reaction 개수, (no-reaction 제외) BH2[4-10],
//             (no-reaction 집합에서) Trig1 = BH2[4-10] ⊗ HTOF mult≥2
//   - 그림: (no-reaction 집합 & origin!=1)에서
//           hNpiMinus, hNpiPlus, h2_NpiMinus_vs_NpiPlus, h2_NpiMinus_vs_Np
//
// 임계값(기본): 0.1 MIP, dE/dx=2 MeV/cm, BH2=5 mm, HTOF=10 mm
// ID정책: HTOF tile = TParticle::StatusCode(), BH2 seg = 좌표→세그 매핑(MapBH2_WorldToSeg)
// edep 정책: *_edep 브랜치 우선, 없으면 TParticle::Weight(); 합계가 thr 이상일 때만 유효 히트

#include "TFile.h"
#include "TTree.h"
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
#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <type_traits>

namespace TPI {

// --- geometry (E72-like) ---
static const int    kNBH2Seg   = 15;     // 0..14
static const double kBH2_x0    = -10.0;  // [mm]
static const double kBH2_pitch = 14.0;   // [mm]
static const double kBH2_total = kNBH2Seg * kBH2_pitch;
static const int    kNHTOF     = 34;     // 0..33

// dictionary for vector<TParticle>
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

static inline bool HasBranch(TTree* tr, const char* b){ return tr && tr->GetBranch(b); }

// helper: try multiple branch names
template<class T>
static bool AttachBranch(TTree* T, const std::vector<const char*>& names, T** ptrOut, TString& chosen){
  for(const char* nm : names){
    if(HasBranch(T,nm)){ T->SetBranchAddress(nm, ptrOut); chosen = nm; return true; }
  }
  return false;
}

static inline double pct(long long a,long long b){ return (b>0)? 100.0*double(a)/double(b) : 0.0; }

} // namespace TPI


void two_pi_study(const char* filename,
                  const char* treename="g4hyptpc",
                  int bh2_lo=4, int bh2_hi=10,
                  bool savePlots=true,
                  double mipFrac=0.1,
                  double mipMeVperCm=2.0,
                  double BH2_thickness_mm=5.0,
                  double HTOF_thickness_mm=10.0)
{
  using namespace TPI;

  // --- open
  TFile* f=TFile::Open(filename,"READ");
  if(!f||f->IsZombie()){ std::cerr<<"[ERR] open "<<filename<<" failed\n"; return; }
  TTree* T=(TTree*)f->Get(treename);
  if(!T){ std::cerr<<"[ERR] tree "<<treename<<" not found\n"; return; }

  // --- required detector branches
  if(!HasBranch(T,"BH2") || !HasBranch(T,"HTOF")){
    std::cerr<<"[ERR] need branches 'BH2' and 'HTOF' (vector<TParticle>)\n"; return;
  }
  const bool hasBH2edep  = HasBranch(T,"BH2_edep");
  const bool hasHTOFedep = HasBranch(T,"HTOF_edep");

  std::vector<TParticle>* BH2=nullptr;  T->SetBranchAddress("BH2",&BH2);
  std::vector<TParticle>* HTOF=nullptr; T->SetBranchAddress("HTOF",&HTOF);
  std::vector<double>* BH2_edep=nullptr;  if(hasBH2edep)  T->SetBranchAddress("BH2_edep",&BH2_edep);
  std::vector<double>* HTOF_edep=nullptr; if(hasHTOFedep) T->SetBranchAddress("HTOF_edep",&HTOF_edep);

  // --- optional truth branches for reaction & beam tagging
  // origin: 0=unknown, 1=kPrimaryBeam, 2=kForced2PiChild(반응 자식)
  std::vector<int>*   vOrigin=nullptr; TString bOrigin;
  std::vector<int>*   vPDG   =nullptr; TString bPDG;

  bool hasOrigin = AttachBranch(T, {"track_origin","trk_origin","origin"}, &vOrigin, bOrigin);
  bool hasPDG    = AttachBranch(T, {"track_pdg","trk_pdg","pdg"},         &vPDG,    bPDG);

  if(!hasOrigin){
    std::cerr<<"[WARN] origin branch not found (track_origin/trk_origin/origin). "
                "No-reaction set will be empty; counts still printed.\n";
  }
  if(!hasPDG){
    std::cerr<<"[WARN] PDG branch not found (track_pdg/trk_pdg/pdg). "
                "Pi histograms will be skipped.\n";
  }

  // --- thresholds
  const double thrBH2  = MIPThrMeV(mipFrac,mipMeVperCm,BH2_thickness_mm);   // e.g. 0.10 MeV
  const double thrHTOF = MIPThrMeV(mipFrac,mipMeVperCm,HTOF_thickness_mm);  // e.g. 0.20 MeV

  std::cout<<std::fixed<<std::setprecision(3);
  std::cout<<"[INFO] BH2 thr="<<thrBH2<<" MeV, HTOF thr="<<thrHTOF<<" MeV | "
           <<"BH2 range = "<<bh2_lo<<"-"<<bh2_hi<<"\n";
  std::cout<<"[INFO] ID policy: HTOF=StatusCode, BH2=coord-map; edep: *_edep→Weight; valid if sum>=thr.\n";

  // --- counters
  Long64_t N_total=0;
  Long64_t N_BH2_4_10_all=0;

  Long64_t N_noReaction=0;
  Long64_t N_BH2_4_10_in_noReaction=0;
  Long64_t N_BH2_4_10_reactionOnly=0; // derived later

  Long64_t N_trig1_in_noReaction=0;   // (no-reaction) & [BH2 in 4-10] & [HTOF mult>=2]

  // --- histograms (no-reaction & origin!=1 에서)
  TH1I* hNpiMinus=nullptr;
  TH1I* hNpiPlus =nullptr;
  TH2I* h2_Npim_vs_Piplus=nullptr;
  TH2I* h2_Npim_vs_P   =nullptr;

  if(hasPDG){
    hNpiMinus = new TH1I("hNpiMinus_noReaction_noBeam","N(pi-) per event | no-reaction & origin!=1;N(pi-);Events",20,-0.5,19.5);
    hNpiPlus  = new TH1I("hNpiPlus_noReaction_noBeam", "N(pi+) per event | no-reaction & origin!=1;N(pi+);Events",20,-0.5,19.5);
    h2_Npim_vs_Piplus = new TH2I("h2_Npim_vs_Piplus_noReaction_noBeam","N(pi-) vs N(pi+) | no-reaction & origin!=1;N(pi+);N(pi-)",
                                 20,-0.5,19.5, 20,-0.5,19.5);
    h2_Npim_vs_P      = new TH2I("h2_Npim_vs_P_noReaction_noBeam","N(pi-) vs N(p) | no-reaction & origin!=1;N(p);N(pi-)",
                                 20,-0.5,19.5, 20,-0.5,19.5);
    hNpiMinus->SetDirectory(nullptr);
    hNpiPlus ->SetDirectory(nullptr);
    h2_Npim_vs_Piplus->SetDirectory(nullptr);
    h2_Npim_vs_P     ->SetDirectory(nullptr);
  }

  const Long64_t N=T->GetEntries();
  for(Long64_t ie=0; ie<N; ++ie){
    T->GetEntry(ie); N_total++;

    // ---- BH2 유효 세그 집합 (edep 합 ≥ thr)
    std::map<int,double> bh2E;
    if(BH2){
      for(size_t i=0;i<BH2->size();++i){
        const TParticle& p=BH2->at(i);
        int sid = MapBH2_WorldToSeg(p.Vx(),p.Vy(),p.Vz());
        if(0<=sid && sid<kNBH2Seg){
          const double ed = (hasBH2edep && BH2_edep && i<BH2_edep->size()) ? BH2_edep->at(i)
                                                                            : p.GetWeight();
          bh2E[sid]+=ed;
        }
      }
    }
    std::set<int> bh2Valid;
    for(auto& kv:bh2E) if(kv.second>=thrBH2) bh2Valid.insert(kv.first);

    const bool bh2_in_4_10 = std::any_of(bh2Valid.begin(), bh2Valid.end(),
                                         [&](int s){ return (bh2_lo<=s && s<=bh2_hi); });
    if(bh2_in_4_10) N_BH2_4_10_all++;

    // ---- HTOF 유효 타일 집합
    std::map<int,double> htofE;
    if(HTOF){
      for(size_t i=0;i<HTOF->size();++i){
        const TParticle& p=HTOF->at(i);
        int tid = p.GetStatusCode();
        if(0<=tid && tid<kNHTOF){
          const double ed = (hasHTOFedep && HTOF_edep && i<HTOF_edep->size()) ? HTOF_edep->at(i)
                                                                               : p.GetWeight();
          htofE[tid]+=ed;
        }
      }
    }
    std::set<int> htofValid;
    for(auto& kv:htofE) if(kv.second>=thrHTOF) htofValid.insert(kv.first);
    const int htofMult = (int)htofValid.size();

    // ---- reaction 판정(origin==2 존재?) & 원빔 제외용(origin==1)
    bool hasForced2Pi=false;
    if(hasOrigin && vOrigin){
      for(int o : *vOrigin){ if(o==2){ hasForced2Pi=true; break; } }
    }

    // (2) no-reaction 집합(=타겟 미도달로 2π 생성 없음)
    bool isNoReaction = hasOrigin ? (!hasForced2Pi) : false;
    if(isNoReaction){
      N_noReaction++;
      if(bh2_in_4_10) N_BH2_4_10_in_noReaction++;

      // (3) no-reaction 기준 Trig1
      if(bh2_in_4_10 && htofMult>=2) N_trig1_in_noReaction++;

      // (4) no-reaction & 원빔 제외(origin!=1)에서 PDG 히스토그램
      if(hasPDG && vPDG && vOrigin){
        int nPim=0, nPip=0, nP=0;
        const size_t ntrk = std::min(vPDG->size(), vOrigin->size());
        for(size_t it=0; it<ntrk; ++it){
          if(vOrigin->at(it)==1) continue; // exclude primary beam
          const int pdg = vPDG->at(it);
          if(pdg==-211) nPim++;
          else if(pdg== 211) nPip++;
          else if(pdg==2212) nP++;
        }
        if(hNpiMinus){ hNpiMinus->Fill(nPim); }
        if(hNpiPlus ){ hNpiPlus ->Fill(nPip); }
        if(h2_Npim_vs_Piplus){ h2_Npim_vs_Piplus->Fill(nPip, nPim); }
        if(h2_Npim_vs_P     ){ h2_Npim_vs_P     ->Fill(nP,   nPim); }
      }
    }
  } // end loop

  N_BH2_4_10_reactionOnly = N_BH2_4_10_all - N_BH2_4_10_in_noReaction;

  // ---- print summary ----
  std::cout<<"\n===== 2pi_study summary (BH2 "<<bh2_lo<<"-"<<bh2_hi<<") =====\n";
  std::cout<<"Total events                                 : "<<N_total<<"\n";
  std::cout<<"BH2[4-10] hit events (all)                   : "<<N_BH2_4_10_all
           <<"  ("<<pct(N_BH2_4_10_all,N_total)<<" % of Total)\n";
  std::cout<<"No-reaction events (no origin==2)            : "<<N_noReaction
           <<"  ("<<pct(N_noReaction,N_total)<<" % of Total)\n";
  std::cout<<"BH2[4-10] after excluding no-reaction        : "<<N_BH2_4_10_reactionOnly
           <<"  ("<<pct(N_BH2_4_10_reactionOnly, N_total - N_noReaction)<<" % of [Total - NoReaction])\n";
  std::cout<<"Trig1 in no-reaction set (BH2[4-10] & m>=2)  : "<<N_trig1_in_noReaction
           <<"  ("<<pct(N_trig1_in_noReaction, N_noReaction)<<" % of NoReaction)\n";

  // ---- draw/save plots ----
  if(hasPDG && hNpiMinus){
    TCanvas* c1=new TCanvas("cNpi","N(pi-) & N(pi+) | no-reaction & no-beam",1000,500);
    c1->Divide(2,1);
    c1->cd(1); hNpiMinus->Draw("hist");
    c1->cd(2); hNpiPlus ->Draw("hist");

    TCanvas* c2=new TCanvas("cPairs","N(pi-) pairs | no-reaction & no-beam",1000,500);
    c2->Divide(2,1);
    c2->cd(1); h2_Npim_vs_Piplus->Draw("colz");
    c2->cd(2); h2_Npim_vs_P     ->Draw("colz");

    if(savePlots){
      c1->SaveAs("2pi_noReaction_noBeam_Npi.png");
      c2->SaveAs("2pi_noReaction_noBeam_pairs.png");
    }
  } else {
    std::cout<<"[INFO] PDG/origin branches missing -> pi histograms skipped.\n";
  }
}

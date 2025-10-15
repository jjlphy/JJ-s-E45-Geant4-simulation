// 2pi_study_v3.C
// 요구사항 정리:
// (1) 요약 카운트
//     - Total events
//     - BH2 [lo,hi]에 도달한 이벤트 수
//     - 그 중 "타겟에 닿지 않은(no-reaction)" 이벤트 수  (tgt_touch_flag==0 우선 사용)
// (2) (no-reaction 제외 후) BH2[lo,hi] && HTOF multiplicity>=2 이벤트 수
//     - HTOF 타일 ID는 copy_no 브랜치가 있으면 그것을, 없으면 StatusCode를 사용
// (3) (no-reaction 제외 후, 원빔 제외 가정) 반응입자 HTOF 타일 기반 2D 히트패턴
//     - pi- vs pi+ : 0..33 × 0..33
//     - pi- vs p   : 0..33 × 0..33
//     - 이벤트마다 유효 타일 집합(종별)에서 가능한 쌍(i,j)을 모두 1씩 누적
//
// 참고: 벡터<TParticle> 딕셔너리 생성 필요. (아래 DictGuard에서 자동 처리)
// HTOF 타일 판정: 우선순위 (HTOF_copyNo) → (TParticle::StatusCode()).
// BH2 세그 판정: 좌표→세그 매핑(MapBH2_WorldToSeg). edep-thr: *_edep 또는 Weight fallback.
// 타겟 접촉 여부: 우선 (tgt_touch_flag) 브랜치 사용, 없으면 TGT 히트 합 에너지 thr로 대체.

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
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>

namespace TPI3 {

// --- E72-like geometry constants ---
static const int    kNBH2Seg   = 15;     // BH2 seg indices: 0..14
static const double kBH2_x0    = -10.0;  // [mm]
static const double kBH2_pitch = 14.0;   // [mm]
static const double kBH2_total = kNBH2Seg * kBH2_pitch;
static const int    kNHTOF     = 34;     // HTOF tiles: 0..33

// ROOT needs a dictionary for vector<TParticle>
struct DictGuard {
  DictGuard(){ gSystem->Load("libPhysics");
               gInterpreter->GenerateDictionary("vector<TParticle>","TParticle.h;vector"); }
};
static DictGuard _dg;

// thresholds helper (MIP fraction × dE/dx × thickness[cm])
static inline double MIPThrMeV(double frac,double dEdx,double thk_mm){
  return frac * dEdx * (thk_mm*0.1); // mm → cm
}

// world (x,y,z) → BH2 seg id (0..14), -1 if out
static int MapBH2_WorldToSeg(double x, double, double){
  const double xloc = x - kBH2_x0;
  const double xmin = -0.5*kBH2_total, xmax = 0.5*kBH2_total;
  if (xloc < xmin || xloc >= xmax) return -1;
  int idx = int(std::floor((xloc - xmin)/kBH2_pitch));
  if(idx<0) idx=0; if(idx>=kNBH2Seg) idx=kNBH2Seg-1;
  return idx;
}

// 안전한 브랜치 존재 확인
static inline bool hasBranch(TTree* T, const char* name){
  return T && T->GetBranch(name)!=nullptr;
}

// 안전한 벡터 인덱싱
template<typename V> static inline bool hasIdx(const V* v, size_t i){
  return v && i < v->size();
}

} // namespace TPI3


// ================= PUBLIC ENTRY =================
void two_pi_study_v3(const char* filename,
                     const char* treename="g4hyptpc",
                     int bh2_lo=4, int bh2_hi=10,
                     bool savePlots=true,
                     double mipFrac=0.1,
                     double mipMeVperCm=2.0,
                     double BH2_thickness_mm=5.0,
                     double HTOF_thickness_mm=10.0,
                     double TGT_anyEdep_thr=0.0) // tgt_touch_flag 미존재 시 TGT 합edep≥thr로 대체(0→존재만)
{
  using namespace TPI3;

  // --- open & tree ---
  TFile* f=TFile::Open(filename,"READ");
  if(!f||f->IsZombie()){ std::cerr<<"[ERR] open "<<filename<<" failed\n"; return; }
  TTree* T=(TTree*)f->Get(treename);
  if(!T){ std::cerr<<"[ERR] tree "<<treename<<" not found\n"; return; }

  // 필수 브랜치(데이터 구성에 따라 없는 것도 허용하되 로직상 체크)
  const bool needBH2  = hasBranch(T,"BH2");
  const bool needHTOF = hasBranch(T,"HTOF");
  if(!needBH2 || !needHTOF){
    std::cerr<<"[ERR] need BH2 & HTOF (vector<TParticle>) branches.\n";
    return;
  }
  const bool hasPRM = hasBranch(T,"PRM");
  const bool hasSEC = hasBranch(T,"SEC");
  const bool hasTGT = hasBranch(T,"TGT");

  const bool hasBH2edep  = hasBranch(T,"BH2_edep");
  const bool hasHTOFedep = hasBranch(T,"HTOF_edep");
  const bool hasTGTedep  = hasBranch(T,"TGT_edep");

  const bool hasCopyNo   = hasBranch(T,"HTOF_copyNo");
  const bool hasTgtFlag  = hasBranch(T,"tgt_touch_flag");

  // --- set addresses ---
  std::vector<TParticle>* BH2=nullptr;   T->SetBranchAddress("BH2",&BH2);
  std::vector<TParticle>* HTOF=nullptr;  T->SetBranchAddress("HTOF",&HTOF);
  std::vector<TParticle>* PRM=nullptr;   if(hasPRM) T->SetBranchAddress("PRM",&PRM);
  std::vector<TParticle>* SEC=nullptr;   if(hasSEC) T->SetBranchAddress("SEC",&SEC);
  std::vector<TParticle>* TGT=nullptr;   if(hasTGT) T->SetBranchAddress("TGT",&TGT);

  std::vector<double>* BH2_edep=nullptr; if(hasBH2edep)  T->SetBranchAddress("BH2_edep",&BH2_edep);
  std::vector<double>* HTOF_edep=nullptr;if(hasHTOFedep) T->SetBranchAddress("HTOF_edep",&HTOF_edep);
  std::vector<double>* TGT_edep=nullptr; if(hasTGTedep)  T->SetBranchAddress("TGT_edep",&TGT_edep);

  std::vector<int>* HTOF_copyNo=nullptr; if(hasCopyNo) T->SetBranchAddress("HTOF_copyNo",&HTOF_copyNo);

  int tgt_touch_flag=0;
  if(hasTgtFlag) T->SetBranchAddress("tgt_touch_flag",&tgt_touch_flag);

  // --- thresholds ---
  const double thrBH2  = MIPThrMeV(mipFrac,mipMeVperCm,BH2_thickness_mm);   // ex) 0.10 MeV
  const double thrHTOF = MIPThrMeV(mipFrac,mipMeVperCm,HTOF_thickness_mm);  // ex) 0.20 MeV

  std::cout<<std::fixed<<std::setprecision(3);
  std::cout<<"[INFO] BH2 thr="<<thrBH2<<" MeV, HTOF thr="<<thrHTOF
           <<" MeV | BH2 range = "<<bh2_lo<<"-"<<bh2_hi<<"\n";
  if(hasTgtFlag)
    std::cout<<"[INFO] TGT-hit policy: use tgt_touch_flag(1=hit)\n";
  else
    std::cout<<"[INFO] TGT-hit policy: sum(TGT edep) >= "<<TGT_anyEdep_thr<<" MeV (0 → presence-only)\n";
  if(hasCopyNo)
    std::cout<<"[INFO] HTOF tile ID: using HTOF_copyNo\n";
  else
    std::cout<<"[WARN] HTOF_copyNo not found — falling back to StatusCode()\n";

  // --- counters ---
  Long64_t N_total=0;
  Long64_t N_BH2_in=0;                // BH2 [lo,hi] 도달
  Long64_t N_noTGT=0;                 // 그 중 타겟 미접촉
  Long64_t N_HTOFm2_after=0;          // (noTGT 제외 후) BH2[lo,hi] & HTOF mult>=2

  // --- histograms (0..33 x 0..33) ---
  TH2I* h2_pim_vs_piplus = new TH2I("h2_piMinus_vs_piPlus",
    "HTOF hitpattern: #pi^{-} tiles vs #pi^{+} tiles (reaction-only);HTOF tile (pi+);HTOF tile (pi-)",
    kNHTOF, -0.5, kNHTOF-0.5, kNHTOF, -0.5, kNHTOF-0.5);
  TH2I* h2_pim_vs_p = new TH2I("h2_piMinus_vs_p",
    "HTOF hitpattern: #pi^{-} tiles vs p tiles (reaction-only);HTOF tile (p);HTOF tile (pi-)",
    kNHTOF, -0.5, kNHTOF-0.5, kNHTOF, -0.5, kNHTOF-0.5);
  h2_pim_vs_piplus->SetDirectory(nullptr);
  h2_pim_vs_p     ->SetDirectory(nullptr);

  auto inBH2Range = [&](const std::set<int>& S){
    for(int s:S){ if(bh2_lo<=s && s<=bh2_hi) return true; }
    return false;
  };

  const Long64_t N=T->GetEntries();
  for(Long64_t ie=0; ie<N; ++ie){
    T->GetEntry(ie); N_total++;

    // ---------- BH2 유효 세그 집합 ----------
    std::map<int,double> bh2E; std::set<int> bh2Any;
    if(BH2){
      for(size_t i=0;i<BH2->size();++i){
        const TParticle& p=BH2->at(i);
        int sid = MapBH2_WorldToSeg(p.Vx(),p.Vy(),p.Vz());
        if(0<=sid && sid<kNBH2Seg){
          double ed = (hasBH2edep && hasIdx(BH2_edep,i)) ? BH2_edep->at(i) : p.GetWeight();
          bh2E[sid]+=ed; bh2Any.insert(sid);
        }
      }
    }
    std::set<int> bh2Valid;
    bool anyBH2=false; for(auto& kv:bh2E){ if(kv.second>0){ anyBH2=true; break; } }
    if(anyBH2){ for(auto& kv:bh2E){ if(kv.second>=thrBH2) bh2Valid.insert(kv.first); } }
    else       { bh2Valid = bh2Any; }

    const bool passBH2 = (!bh2Valid.empty()) && inBH2Range(bh2Valid);
    if(passBH2) N_BH2_in++; else continue;

    // ---------- 타겟 접촉 여부 ----------
    bool hasTGT = false;
    if(hasTgtFlag){
      hasTGT = (tgt_touch_flag!=0);
    }else{
      double E=0.0;
      if(TGT){
        for(size_t i=0;i<TGT->size();++i){
          double ed = (hasTGTedep && hasIdx(TGT_edep,i)) ? TGT_edep->at(i)
                                                         : TGT->at(i).GetWeight();
          E += ed;
        }
      }
      hasTGT = (TGT && !TGT->empty() && E>=TGT_anyEdep_thr);
    }
    if(!hasTGT){ N_noTGT++; continue; } // (1) no-reaction 집합 제외

    // ---------- HTOF 유효 타일 (종별 분리) ----------
    // ID: copy_no(선호) 또는 StatusCode
    std::map<int,double> ed_piMinus, ed_piPlus, ed_proton;  // 타일별 합edep
    std::set<int> any_piMinus, any_piPlus, any_proton;      // fallback(에너지 정보가 전혀 없을 때)
    int m_allValid=0;

    if(HTOF){
      for(size_t i=0;i<HTOF->size();++i){
        const TParticle& hp = HTOF->at(i);
        int tile = -1;
        if(hasCopyNo && hasIdx(HTOF_copyNo,i)) tile = HTOF_copyNo->at(i);
        else                                   tile = hp.GetStatusCode();
        if(tile<0 || tile>=kNHTOF) continue;

        const int pdg = hp.GetPdgCode();
        double ed = (hasHTOFedep && hasIdx(HTOF_edep,i)) ? HTOF_edep->at(i) : hp.GetWeight();

        // 종별 누적
        if(pdg==-211){ ed_piMinus[tile]+=ed; any_piMinus.insert(tile); }
        else if(pdg==+211){ ed_piPlus[tile]+=ed; any_piPlus.insert(tile); }
        else if(pdg==2212){ ed_proton[tile]+=ed; any_proton.insert(tile); }
      }
    }

    auto makeValidSet = [&](const std::map<int,double>& M, const std::set<int>& ANY)->std::set<int>{
      std::set<int> S; bool any=false; for(auto& kv:M){ if(kv.second>0){ any=true; break; } }
      if(any){ for(auto& kv:M){ if(kv.second>=thrHTOF) S.insert(kv.first); } }
      else    S = ANY; // 에너지 정보가 전혀 없으면 존재만으로
      return S;
    };

    const std::set<int> tiles_pim = makeValidSet(ed_piMinus, any_piMinus);
    const std::set<int> tiles_pip = makeValidSet(ed_piPlus , any_piPlus );
    const std::set<int> tiles_p   = makeValidSet(ed_proton , any_proton );

    // HTOF 전체 multiplicity(m>=2) 판단은 모든 유효 타일의 합집합으로
    std::set<int> tiles_all;
    tiles_all.insert(tiles_pim.begin(), tiles_pim.end());
    tiles_all.insert(tiles_pip.begin(), tiles_pip.end());
    tiles_all.insert(tiles_p.begin()  , tiles_p.end());
    m_allValid = (int)tiles_all.size();

    if(m_allValid>=2) N_HTOFm2_after++; // (2)

    // ---------- (3) 2D 히트패턴 누적 ----------
    // “원빔 제외”는 강제 2π/가드 로직에서 원빔이 HTOF까지 오지 않는다는 전제(=reaction hits만)로 처리.
    // 이벤트 내 종별 타일 집합의 모든 조합을 1씩 누적
    for(int j: tiles_pip){      // x축: pi+
      for(int i: tiles_pim){    // y축: pi-
        h2_pim_vs_piplus->Fill(j,i);
      }
    }
    for(int j: tiles_p){        // x축: p
      for(int i: tiles_pim){    // y축: pi-
        h2_pim_vs_p->Fill(j,i);
      }
    }
  } // loop

  // --- summary ---
  std::cout<<"\n===== 2pi_study_v3 summary (BH2 "<<bh2_lo<<"-"<<bh2_hi<<") =====\n";
  std::cout<<"Total events                                 : "<<N_total<<"\n";
  std::cout<<"BH2["<<bh2_lo<<"-"<<bh2_hi<<"] hit events (all)        : "<<N_BH2_in
           <<"  ("<<(N_total?100.0*double(N_BH2_in)/double(N_total):0.0)<<" % of Total)\n";
  std::cout<<"No-reaction (no TGT hit)                     : "<<N_noTGT
           <<"  ("<<(N_total?100.0*double(N_noTGT)/double(N_total):0.0)<<" % of Total)\n";
  const Long64_t denom_after = (N_BH2_in - N_noTGT);
  std::cout<<"(After removing no-reaction) BH2&HTOF(m>=2)  : "<<N_HTOFm2_after
           <<"  ("<<(denom_after>0?100.0*double(N_HTOFm2_after)/double(denom_after):0.0)
           <<" % of [BH2-in - NoReaction])\n";

  // --- draw & save ---
  TCanvas* cA=new TCanvas("cHTOF_pim_vs_piplus","HTOF pi- vs pi+ (reaction-only)",900,800);
  h2_pim_vs_piplus->Draw("colz");
  TCanvas* cB=new TCanvas("cHTOF_pim_vs_p","HTOF pi- vs p (reaction-only)",900,800);
  h2_pim_vs_p->Draw("colz");

  if(savePlots){
    cA->SaveAs("HTOF_hitpattern_piMinus_vs_piPlus.png");
    cB->SaveAs("HTOF_hitpattern_piMinus_vs_p.png");
  }
}

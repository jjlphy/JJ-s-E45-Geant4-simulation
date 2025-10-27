// -*- C++ -*-
// HTOF_pair_reactionChg_10_26.C  (2025-10-26, edited for jaejin)
// 목적: "리액션 하전입자(강제 2π 자식)" 히트만으로 HTOF 인접 페어 카운트.
// 선택 파이프라인:
//   1) tgt_touch_flag == 1
//   2) BH2 in [bh2_lo, bh2_hi] (분모)
//   3) HTOF: origin==kForced2PiChild 인 하전 히트만 사용 (유효=ΣEdep≥thr)
//      * UniqueID(origin)가 없으면 자동 탐지 → 경고 출력 후 폴백:
//        "하전 전부(원빔 포함 가능)" 로 집계(참고용)
//
// 사용법:
//   root -l
//   .L HTOF_pair_reactionChg_10_26.C+
//   HTOF_pair_reactionChg_10_26("E45_Nov_piplusn_098.root","g4hyptpc",4,10,0.10,2.0,5.0,10.0,true,"BH2_4_10");

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1I.h"
#include "TSystem.h"
#include "TInterpreter.h"
#include "TParticle.h"
#include "TString.h"
#include <vector>
#include <map>
#include <set>
#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>

namespace HPRC {

struct DictGuard {
  DictGuard(){
    gSystem->Load("libPhysics");
    gInterpreter->GenerateDictionary("vector<TParticle>","TParticle.h;vector");
  }
};
static DictGuard _dg;

// BH2 mapper (E72-like)
static const int    kNBH2Seg   = 15;     // 0..14
static const double kBH2_x0    = -10.0;  // [mm]
static const double kBH2_pitch = 14.0;   // [mm]
static const double kBH2_total = kNBH2Seg * kBH2_pitch;

static inline int MapBH2_WorldToSeg(double x, double, double){
  const double xloc = x - kBH2_x0;
  const double xmin = -0.5*kBH2_total, xmax = 0.5*kBH2_total;
  if (xloc < xmin || xloc >= xmax) return -1;
  int idx = int(std::floor((xloc - xmin)/kBH2_pitch));
  if (idx<0) idx=0;
  if (idx>=kNBH2Seg) idx=kNBH2Seg-1;
  return idx;
}

static inline double MIPThrMeV(double frac,double dEdx,double thk_mm){
  return frac * dEdx * (thk_mm*0.1); // mm → cm
}

// EOrigin codes (TrackTag.hh와 일치해야 함)
static const int ORI_UNKNOWN          = 0;
static const int ORI_PRIMARY_BEAM     = 2; // 프로젝트에서 실제 값 확인 필요
static const int ORI_FORCED2PI_CHILD  = 1; // kForced2PiChild (프로젝트 정의와 동기화)

// 하전 여부(보조): HTOFSD가 이미 neutral skip하므로 보통 필요없지만 방어용
static inline bool IsChargedByPDG(int pdg){
  // 빠른 커버: ±211, ±321, ±2212, ±11, ±13 ...
  if(pdg==0) return true; // PDG 미기록이면 SD에서 이미 charge 0 배제했다고 가정
  const int a = std::abs(pdg);
  if(a==211 || a==321 || a==2212 || a==11 || a==13) return true;
  // 나머지는 대충 홀수면 대체로 하전인 경우가 많지만, 확실치 않음 → true로 관대하게
  return true;
}

} // ns

void HTOF_pair_reactionChg_10_26(const char* filename="E45.root",
                                 const char* treename="g4hyptpc",
                                 int    bh2_lo=4, int bh2_hi=10,
                                 double mipFrac=0.10,
                                 double mipMeVperCm=2.0,
                                 double BH2_thickness_mm=5.0,
                                 double HTOF_thickness_mm=10.0,
                                 bool   save=true,
                                 const char* tag="BH2_4_10")
{
  using namespace HPRC;

  // ----- open -----
  TFile* f=TFile::Open(filename,"READ");
  if(!f||f->IsZombie()){ std::cerr<<"[ERR] open "<<filename<<" failed\n"; return; }
  TTree* T=(TTree*)f->Get(treename);
  if(!T){ std::cerr<<"[ERR] tree "<<treename<<" not found\n"; return; }

  if(!T->GetBranch("BH2") || !T->GetBranch("HTOF") || !T->GetBranch("tgt_touch_flag")){
    std::cerr<<"[ERR] need BH2, HTOF (vector<TParticle>) + tgt_touch_flag(int)\n"; return;
  }

  const bool hasBH2edep  = (T->GetBranch("BH2_edep")  != nullptr);
  const bool hasHTOFedep = (T->GetBranch("HTOF_edep") != nullptr);

  std::vector<TParticle>* BH2=nullptr;   T->SetBranchAddress("BH2",&BH2);
  std::vector<TParticle>* HTOF=nullptr;  T->SetBranchAddress("HTOF",&HTOF);
  std::vector<double>* BH2_edep=nullptr;   if(hasBH2edep)  T->SetBranchAddress("BH2_edep",&BH2_edep);
  std::vector<double>* HTOF_edep=nullptr;  if(hasHTOFedep) T->SetBranchAddress("HTOF_edep",&HTOF_edep);
  int tgt_touch_flag=0; T->SetBranchAddress("tgt_touch_flag",&tgt_touch_flag);

  const double thrBH2  = MIPThrMeV(mipFrac,mipMeVperCm,BH2_thickness_mm);
  const double thrHTOF = MIPThrMeV(mipFrac,mipMeVperCm,HTOF_thickness_mm);

  std::cout<<std::fixed<<std::setprecision(3)
           <<"[INFO] BH2 thr="<<thrBH2<<" MeV, HTOF thr="<<thrHTOF<<" MeV | "
           <<"BH2 range = "<<bh2_lo<<"-"<<bh2_hi<<"\n";
  std::cout<<"[INFO] Using ONLY reaction charged hits (origin==kForced2PiChild). "
              "If origin tag is missing, fallback=ALL charged (beam may leak).\n";

  const std::pair<int,int> pairs[] = {
    {13,14},{14,15},{15,16},{16,17},{17,18},{18,19},{19,20},
    {20,21},{21,22},{22,23},{23,24},{24,25},{25,26}
  };
  const int NP = sizeof(pairs)/sizeof(pairs[0]);

  Long64_t N_total=0, N_not_tgt=0, N_tgt=0, N_bh2In=0;
  Long64_t N_mp_ge2=0, N_mp_gt2=0;
  std::vector<Long64_t> pairCnt_ge2(NP,0), pairCnt_gt2(NP,0);

  TH1I* hPairs_ge2 = new TH1I("hPairs_ge2_react",
      Form("HTOF pairs (reaction charged) MP#geq2 | BH2 %d-%d;Pair;Events",bh2_lo,bh2_hi),
      NP, 0.5, NP+0.5);
  TH1I* hPairs_gt2 = new TH1I("hPairs_gt2_react",
      Form("HTOF pairs (reaction charged) MP>2   | BH2 %d-%d;Pair;Events",bh2_lo,bh2_hi),
      NP, 0.5, NP+0.5);
  hPairs_ge2->SetDirectory(nullptr);
  hPairs_gt2->SetDirectory(nullptr);
  for(int i=1;i<=NP;++i){
    hPairs_ge2->GetXaxis()->SetBinLabel(i, Form("(%d,%d)", pairs[i-1].first, pairs[i-1].second));
    hPairs_gt2->GetXaxis()->SetBinLabel(i, Form("(%d,%d)", pairs[i-1].first, pairs[i-1].second));
  }

  // ---- OriginTag 사용 가능 여부 빠른 감지 ----
  bool originSeen = false;  // 한 번이라도 ORI_FORCED2PI_CHILD을 보면 true
  {
    const Long64_t Npeek = std::min<Long64_t>(T->GetEntries(), 2000);
    for(Long64_t ie=0; ie<Npeek && !originSeen; ++ie){
      T->GetEntry(ie);
      if(!HTOF) continue;
      for(size_t i=0;i<HTOF->size();++i){
        const int ori = HTOF->at(i).GetUniqueID();
        if(ori == ORI_FORCED2PI_CHILD){ originSeen = true; break; }
      }
    }
    if(!originSeen){
      std::cerr<<"[WARN] No HTOF hit has UniqueID==kForced2PiChild in first "<<Npeek
               <<" events. Likely UniqueID not written in SD. "
               <<"FALLBACK to 'ALL charged hits' (beam may be included).\n";
    }
  }

  auto pct = [](long long a, long long b)->double{ return (b>0)? (100.0*(double)a/(double)b):0.0; };

  const Long64_t N = T->GetEntries();
  for(Long64_t ie=0; ie<N; ++ie){
    T->GetEntry(ie); N_total++;

    if(tgt_touch_flag!=1){ N_not_tgt++; continue; }
    N_tgt++;

    // ---- BH2 in-range (denominator) ----
    std::map<int,double> bh2E;
    if(BH2){
      for(size_t i=0;i<BH2->size();++i){
        const TParticle& p = BH2->at(i);
        const int sid = MapBH2_WorldToSeg(p.Vx(), p.Vy(), p.Vz());
        if(0<=sid && sid<kNBH2Seg){
          const double ed = (hasBH2edep && BH2_edep && i<BH2_edep->size())
                            ? BH2_edep->at(i) : p.GetWeight();
          bh2E[sid]+=ed;
        }
      }
    }
    bool bh2_in=false;
    for(const auto& kv:bh2E){
      if(kv.second>=thrBH2){
        const int s=kv.first;
        if(bh2_lo<=s && s<=bh2_hi){ bh2_in=true; break; }
      }
    }
    if(!bh2_in) continue;
    N_bh2In++;

    // ---- HTOF: reaction charged only (or fallback) ----
    std::map<int,double> tileE;
    if(HTOF){
      for(size_t i=0;i<HTOF->size();++i){
        const TParticle& p = HTOF->at(i);
        const int tile = p.GetStatusCode(); // copy-no
        if(tile<0 || tile>33) continue;

        const int ori = p.GetUniqueID();
        const int pdg = p.GetPdgCode(); // 방어용

        bool use=false;
        if(originSeen){
          // 정상 경로: 리액션 자식만
          use = (ori == ORI_FORCED2PI_CHILD) && IsChargedByPDG(pdg);
        }else{
          // 폴백: 하전 전부(원빔 포함 가능)
          use = IsChargedByPDG(pdg);
        }
        if(!use) continue;

        const double ed = (hasHTOFedep && HTOF_edep && i<HTOF_edep->size())
                          ? HTOF_edep->at(i) : p.GetWeight();
        tileE[tile]+=ed;
      }
    }

    std::set<int> tilesValid;
    for(const auto& kv:tileE) if(kv.second>=thrHTOF) tilesValid.insert(kv.first);
    const int mult = (int)tilesValid.size();

    if(mult>=2){
      N_mp_ge2++;
      for(int ip=0; ip<NP; ++ip){
        const auto [a,b]=pairs[ip];
        if(tilesValid.count(a) && tilesValid.count(b)){
          pairCnt_ge2[ip]++; hPairs_ge2->Fill(ip+1);
        }
      }
    }
    if(mult>2){
      N_mp_gt2++;
      for(int ip=0; ip<NP; ++ip){
        const auto [a,b]=pairs[ip];
        if(tilesValid.count(a) && tilesValid.count(b)){
          pairCnt_gt2[ip]++; hPairs_gt2->Fill(ip+1);
        }
      }
    }
  }

  // ---- summary ----
  std::cout<<"\n==== Summary (Reaction charged ONLY"
           <<(originSeen? "" : " | FALLBACK=ALL charged")<<") | BH2 "<<bh2_lo<<"-"<<bh2_hi<<" ====\n";
  std::cout<<"Total events                 : "<<N_total<<"\n";
  std::cout<<"Target-not-touched (excluded): "<<N_not_tgt
           <<"  ("<<std::setprecision(3)<<pct(N_not_tgt,N_total)<<" % of Total)\n";
  std::cout<<"Target-touched (used)        : "<<(N_total - N_not_tgt)<<"\n";
  std::cout<<"BH2 in range (denominator)   : "<<N_bh2In
           <<"  ("<<pct(N_bh2In,(N_total - N_not_tgt))<<" % of tgt-touched)\n";

  auto pr = [&](const char* cap, long long Nset, const std::vector<Long64_t>& cnt){
    std::cout<<cap<<Nset<<"  ("<<pct(Nset,N_bh2In)<<" % of BH2-range)\n";
    for(size_t ip=0; ip<cnt.size(); ++ip){
      const auto& pr = pairs[ip];
      std::cout<<"  Pair "<<std::setw(2)<<pr.first<<","<<std::setw(2)<<pr.second
               <<" : "<<cnt[ip]
               <<"  ("<<pct(cnt[ip],Nset)<<" % of "<<cap<<")\n";
    }
  };
  pr("HTOF multiplicity >= 2  : ", N_mp_ge2, pairCnt_ge2);
  pr("HTOF multiplicity >  2  : ", N_mp_gt2, pairCnt_gt2);

  // ---- draw ----
  TCanvas* c1 = new TCanvas("cPairs_ge2_react","Pairs (reaction charged) MP#geq2", 1000, 450);
  hPairs_ge2->SetLineWidth(2); hPairs_ge2->LabelsOption("v","X"); hPairs_ge2->Draw("hist");
  TCanvas* c2 = new TCanvas("cPairs_gt2_react","Pairs (reaction charged) MP>2", 1000, 450);
  hPairs_gt2->SetLineWidth(2); hPairs_gt2->LabelsOption("v","X"); hPairs_gt2->Draw("hist");

  if(save){
    TString t(tag); if(t.IsNull()) t = Form("BH2_%d_%d", bh2_lo, bh2_hi);
    c1->SaveAs("HTOF_pairs_react_MPge2_"+t+".png");
    c2->SaveAs("HTOF_pairs_react_MPgt2_"+t+".png");
  }
}

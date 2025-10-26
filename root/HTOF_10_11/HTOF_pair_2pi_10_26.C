// HTOF_pair_2pi_10_26.C
// 목적: "리액션 π-"(forced-2π child) 히트만으로 HTOF 인접 페어 카운트.
// 선택 파이프라인:
//   0) 파일/트리 로드
//   1) tgt_touch_flag == 1 (타겟에 닿은 이벤트만 채택; 미도달은 제외/집계)
//   2) BH2 in [bh2_lo, bh2_hi] (분모)
//   3) 리액션 π-(origin==kForced2PiChild && pdg==-211) 히트만으로 HTOF 유효 타일 집합 구성
//   4) MP>=2 / MP>2 별로 인접 페어 카운트(이벤트당 유일)
//
// 카운트 보고:
//   - Total events
//   - Target-not-touched events (제외된 수)
//   - BH2 in-range (분모)
//   - MP>=2 이벤트 수 및 각 페어별 개수
//   - MP>2 이벤트 수 및 각 페어별 개수
//
// 정책/브랜치:
//   - HTOF ID = TParticle::StatusCode()  (copy-no, 0..33)
//   - HTOF Origin = TParticle::UniqueID() (EOrigin 코드; 기본 가정: forced-2π child == 1)
//   - HTOF pdg = TParticle::GetPdgCode()  (반응 π- 필터: -211)
//   - BH2 seg id = 좌표→세그먼트 맵(E72-like, 기존과 동일)
//   - *_edep 브랜치가 있으면 그 합, 없으면 TParticle::Weight() 합 사용
//   - "유효(hit)" = 같은 ID 내 에너지합 ≥ threshold (no fallback)
//   - tgt_touch_flag(int) 브랜치 필요(강제2π 모드에서 세팅)
//
// 사용법:
//   root -l
//   .L HTOF_pair_2pi_10_26.C+
//   HTOF_pair_2pi_10_26("E45_2pi.root","g4hyptpc",
//                       /*bh2_lo=*/4, /*bh2_hi=*/10,
//                       /*mipFrac=*/0.10, /*mipMeVperCm=*/2.0,
//                       /*BH2_thk_mm=*/5.0, /*HTOF_thk_mm=*/10.0,
//                       /*save=*/true, /*tag=*/"BH2_4_10");
//
// 주의: 만약 UniqueID의 Origin 코드가 다르면 ORIGIN_FORCED2PI 상수를 수정하세요.

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1I.h"
#include "TLegend.h"
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

namespace P2PI_10_26 {

// ----- dictionary for vector<TParticle> -----
struct DictGuard {
  DictGuard(){
    gSystem->Load("libPhysics");
    gInterpreter->GenerateDictionary("vector<TParticle>","TParticle.h;vector");
  }
};
static DictGuard _dg;

// ----- BH2 geometry mapping (E72-like) -----
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

// ----- Origin code 가정 (UniqueID) -----
// TrackTag.hh의 EOrigin 할당에 맞춰 수정 가능
static const int ORIGIN_FORCED2PI_CHILD = 1; // kForced2PiChild
static const int PDG_PIM = -211;

} // namespace P2PI_10_26


void HTOF_pair_2pi_10_26(const char* filename="E45_2pi.root",
                          const char* treename="g4hyptpc",
                          int    bh2_lo=4, int bh2_hi=10,
                          double mipFrac=0.10,
                          double mipMeVperCm=2.0,
                          double BH2_thickness_mm=5.0,
                          double HTOF_thickness_mm=10.0,
                          bool   save=true,
                          const char* tag="BH2_4_10")
{
  using namespace P2PI_10_26;

  // ---------------- Open file & tree ----------------
  TFile* f=TFile::Open(filename,"READ");
  if(!f||f->IsZombie()){ std::cerr<<"[ERR] open "<<filename<<" failed\n"; return; }
  TTree* T=(TTree*)f->Get(treename);
  if(!T){ std::cerr<<"[ERR] tree "<<treename<<" not found\n"; return; }

  // Required branches
  if(!T->GetBranch("BH2") || !T->GetBranch("HTOF") || !T->GetBranch("tgt_touch_flag")){
    std::cerr<<"[ERR] need BH2, HTOF (vector<TParticle>) and tgt_touch_flag (int)\n"; return;
  }

  // Optional edep branches
  const bool hasBH2edep  = (T->GetBranch("BH2_edep")  != nullptr);
  const bool hasHTOFedep = (T->GetBranch("HTOF_edep") != nullptr);

  // Set addresses
  std::vector<TParticle>* BH2=nullptr;   T->SetBranchAddress("BH2",&BH2);
  std::vector<TParticle>* HTOF=nullptr;  T->SetBranchAddress("HTOF",&HTOF);
  std::vector<double>* BH2_edep=nullptr;   if(hasBH2edep)  T->SetBranchAddress("BH2_edep",&BH2_edep);
  std::vector<double>* HTOF_edep=nullptr;  if(hasHTOFedep) T->SetBranchAddress("HTOF_edep",&HTOF_edep);
  int tgt_touch_flag=0; T->SetBranchAddress("tgt_touch_flag",&tgt_touch_flag);

  // ---------------- Thresholds & pairs ----------------
  const double thrBH2  = MIPThrMeV(mipFrac,mipMeVperCm,BH2_thickness_mm);   // ~0.10 MeV
  const double thrHTOF = MIPThrMeV(mipFrac,mipMeVperCm,HTOF_thickness_mm);  // ~0.20 MeV

  std::cout<<std::fixed<<std::setprecision(3)
           <<"[INFO] BH2 thr="<<thrBH2<<" MeV, HTOF thr="<<thrHTOF<<" MeV | "
           <<"BH2 range = "<<bh2_lo<<"-"<<bh2_hi<<"\n";
  std::cout<<"[INFO] HTOF tile = StatusCode(copy-no), Origin = UniqueID (kForced2PiChild="<<ORIGIN_FORCED2PI_CHILD<<")\n";
  std::cout<<"[INFO] Counting pairs using ONLY reaction pi- hits (origin==forced2pi & pdg==-211).\n";

  // Pairs to check (13..26 adjacent)
  const std::pair<int,int> pairs[] = {
    {13,14},{14,15},{15,16},{16,17},{17,18},{18,19},{19,20},
    {20,21},{21,22},{22,23},{23,24},{24,25},{25,26}
  };
  const int NP = sizeof(pairs)/sizeof(pairs[0]);

  // ---------------- Accumulators ----------------
  Long64_t N_total=0, N_not_tgt=0, N_tgt=0;
  Long64_t N_bh2In=0;
  Long64_t N_mp_ge2=0, N_mp_gt2=0;

  // Pair counters (unique-per-event per set)
  std::vector<Long64_t> pairCnt_ge2(NP,0), pairCnt_gt2(NP,0);

  // Histograms (optional)
  TH1I* hPairs_ge2 = new TH1I("hPairs_ge2_2pi",
      Form("HTOF pairs (reaction #pi^{-}) MP#geq2 | BH2 %d-%d;Pair;Events",bh2_lo,bh2_hi),
      NP, 0.5, NP+0.5);
  TH1I* hPairs_gt2 = new TH1I("hPairs_gt2_2pi",
      Form("HTOF pairs (reaction #pi^{-}) MP>2 | BH2 %d-%d;Pair;Events",bh2_lo,bh2_hi),
      NP, 0.5, NP+0.5);
  hPairs_ge2->SetDirectory(nullptr);
  hPairs_gt2->SetDirectory(nullptr);
  for(int i=1;i<=NP;++i){
    hPairs_ge2->GetXaxis()->SetBinLabel(i, Form("(%d,%d)", pairs[i-1].first, pairs[i-1].second));
    hPairs_gt2->GetXaxis()->SetBinLabel(i, Form("(%d,%d)", pairs[i-1].first, pairs[i-1].second));
  }

  // ---------------- Event loop ----------------
  const Long64_t N = T->GetEntries();
  for(Long64_t ie=0; ie<N; ++ie){
    T->GetEntry(ie); N_total++;

    // 1) target touch requirement
    if(tgt_touch_flag!=1){ N_not_tgt++; continue; }
    N_tgt++;

    // ---- BH2 valid segments (for denominator) ----
    std::map<int,double> bh2E;
    if(BH2){
      const size_t n = BH2->size();
      for(size_t i=0;i<n;++i){
        const TParticle& p = BH2->at(i);
        const int sid = MapBH2_WorldToSeg(p.Vx(), p.Vy(), p.Vz());
        if(0<=sid && sid<kNBH2Seg){
          const double ed = (hasBH2edep && BH2_edep && i<BH2_edep->size())
                            ? BH2_edep->at(i) : p.GetWeight();
          bh2E[sid] += ed;
        }
      }
    }
    bool bh2_in_range=false;
    for(const auto& kv : bh2E){
      if(kv.second >= thrBH2){
        const int s = kv.first;
        if(bh2_lo<=s && s<=bh2_hi){ bh2_in_range=true; break; }
      }
    }
    if(!bh2_in_range) continue; // event rejected
    N_bh2In++;

    // ---- HTOF valid tiles using ONLY reaction pi- hits ----
    std::map<int,double> tileE; // tile -> sumE
    if(HTOF){
      const size_t n = HTOF->size();
      for(size_t i=0;i<n;++i){
        const TParticle& p = HTOF->at(i);
        const int tile   = p.GetStatusCode();  // copy-no → tile
        if(tile < 0 || tile > 33) continue;

        // 원천 필터: reaction child AND pi-
        const int origin = p.GetUniqueID();   // OriginTag code
        const int pdg    = p.GetPdgCode();    // PDG code at hit
        if(origin != ORIGIN_FORCED2PI_CHILD) continue;
        if(pdg    != PDG_PIM) continue;

        const double ed = (hasHTOFedep && HTOF_edep && i<HTOF_edep->size())
                          ? HTOF_edep->at(i) : p.GetWeight();
        tileE[tile] += ed;
      }
    }

    // Build valid-tile set & multiplicity (reaction π- only)
    std::set<int> tilesValid;
    for(const auto& kv : tileE){
      if(kv.second >= thrHTOF) tilesValid.insert(kv.first);
    }
    const int mult = (int)tilesValid.size();

    // ---- MP ≥ 2 set ----
    if(mult >= 2){
      N_mp_ge2++;
      for(int ip=0; ip<NP; ++ip){
        const auto [a,b] = pairs[ip];
        if(tilesValid.count(a) && tilesValid.count(b)){
          pairCnt_ge2[ip]++; hPairs_ge2->Fill(ip+1);
        }
      }
    }

    // ---- MP > 2 set ----
    if(mult > 2){
      N_mp_gt2++;
      for(int ip=0; ip<NP; ++ip){
        const auto [a,b] = pairs[ip];
        if(tilesValid.count(a) && tilesValid.count(b)){
          pairCnt_gt2[ip]++; hPairs_gt2->Fill(ip+1);
        }
      }
    }
  }

  // ---------------- Print summary ----------------
  auto pct = [](long long a, long long b)->double{ return (b>0)? (100.0*(double)a/(double)b):0.0; };

  std::cout<<"\n==== Summary (Reaction pi- only) | BH2 "<<bh2_lo<<"-"<<bh2_hi<<" ====\n";
  std::cout<<"Total events                 : "<<N_total<<"\n";
  std::cout<<"Target-not-touched (excluded): "<<N_not_tgt
           <<"  ("<<std::setprecision(3)<<pct(N_not_tgt,N_total)<<" % of Total)\n";
  std::cout<<"Target-touched (used)        : "<<(N_total - N_not_tgt)<<"\n";
  std::cout<<"BH2 in range (denominator)   : "<<N_bh2In
           <<"  ("<<pct(N_bh2In, (N_total - N_not_tgt))<<" % of tgt-touched)\n";

  std::cout<<"HTOF multiplicity >= 2       : "<<N_mp_ge2
           <<"  ("<<pct(N_mp_ge2,N_bh2In)<<" % of BH2-range)\n";
  for(int ip=0; ip<NP; ++ip){
    std::cout<<"  Pair "<<std::setw(2)<<pairs[ip].first<<","<<std::setw(2)<<pairs[ip].second
             <<" : "<<pairCnt_ge2[ip]<<"\n";
  }

  std::cout<<"HTOF multiplicity >  2       : "<<N_mp_gt2
           <<"  ("<<pct(N_mp_gt2,N_bh2In)<<" % of BH2-range)\n";
  for(int ip=0; ip<NP; ++ip){
    std::cout<<"  Pair "<<std::setw(2)<<pairs[ip].first<<","<<std::setw(2)<<pairs[ip].second
             <<" : "<<pairCnt_gt2[ip]<<"\n";
  }

  // ---------------- Draw (optional) ----------------
  TCanvas* c1 = new TCanvas("cPairs_ge2_2pi","Pairs (reaction pi-) MP#geq2", 1000, 450);
  hPairs_ge2->SetLineWidth(2);
  hPairs_ge2->LabelsOption("v","X");
  hPairs_ge2->Draw("hist");

  TCanvas* c2 = new TCanvas("cPairs_gt2_2pi","Pairs (reaction pi-) MP>2", 1000, 450);
  hPairs_gt2->SetLineWidth(2);
  hPairs_gt2->LabelsOption("v","X");
  hPairs_gt2->Draw("hist");

  if(save){
    TString t(tag);
    if(t.IsNull()) t = Form("BH2_%d_%d", bh2_lo, bh2_hi);
    c1->SaveAs("HTOF_pairs2pi_MPge2_"+t+".png");
    c2->SaveAs("HTOF_pairs2pi_MPgt2_"+t+".png");
  }
}

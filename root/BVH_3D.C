// BVH_3D_english_dualden.C
// - 분모 A: (BH2 && BVH_U)
// - 분모 B: (BH2 only)
// - 분자:  (U&D가 존재) && (mask 매칭)  // mask는 1st pass에서 event_threshold 이상 셀
// - ecutBH2 / ecutUMeV / ecutDMeV 는 히트 인정용 에너지 컷 (그대로 유지)

//-----작동법__BVH_3D_dualden("E45_BVH1_60mm_54.root", 0.10, 0.04, 0.04, /*event_threshold=*/1);

#include "Rtypes.h"
#include "TFile.h"
#include "TTree.h"
#include "TParticle.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TBox.h"
#include "TLine.h"
#include "TString.h"
#include "TPad.h"
#include "TROOT.h"
#include "TInterpreter.h"
#include "TSystem.h"
#include <vector>
#include <unordered_set>
#include <algorithm>
#include <iostream>
#include <iomanip>

static const int N_BH2  = 15; // 0..14
static const int N_BVHU = 22; // 0..14
static const int N_BVHD = 32; // 0..53

// ===== 유틸: 컷을 넘는 고유 세그먼트만 추출 =====
static inline void get_unique_hits(const std::vector<TParticle>* v,
                                   int nmax, double cutMeV,
                                   std::vector<int>& out)
{
  out.clear();
  if(!v) return;
  std::unordered_set<int> s; s.reserve(8);
  for(const auto& p : *v){
    if(p.GetWeight() <= cutMeV) continue; // weight = Edep[MeV]
    int id = p.GetMother(1);
    if(0<=id && id<nmax) s.insert(id);
  }
  out.assign(s.begin(), s.end());
  std::sort(out.begin(), out.end());
}

void BVH_3D_dualden(const char* fname,
                    double ecutBH2MeV = 0.10,
                    double ecutUMeV   = 0.04,
                    double ecutDMeV   = 0.04,
                    int    event_threshold = 1)
{
  gSystem->Load("libPhysics");
  gInterpreter->GenerateDictionary("vector<TParticle>", "TParticle.h;vector");
  gStyle->SetOptStat(0);

  // 1) 파일/트리
  TFile* f = TFile::Open(fname, "READ");
  if(!f || f->IsZombie()){ std::cerr << "[Error] Cannot open file: " << fname << std::endl; return; }
  TTree* tr = (TTree*)f->Get("g4hyptpc");
  if(!tr){ std::cerr << "[Error] Cannot find TTree 'g4hyptpc'." << std::endl; f->Close(); return; }

  std::vector<TParticle> *BH2=nullptr, *BVHU=nullptr, *BVHD=nullptr;
  tr->SetBranchAddress("BH2",   &BH2);
  tr->SetBranchAddress("BVH_U", &BVHU);
  tr->SetBranchAddress("BVH_D", &BVHD);

  // 2) BH2별 U×D 히스토 (마스크 만들기용)
  std::vector<TH2F*> h_bvh_ud(N_BH2);
  for (int h = 0; h < N_BH2; ++h) {
    TString h_name  = TString::Format("h_bvh_ud_bh2_%d", h);
    TString h_title = TString::Format("BH2=%d;BVH_U Seg;BVH_D Seg", h);
    h_bvh_ud[h] = new TH2F(h_name, h_title,
                           N_BVHU, -0.5, N_BVHU-0.5,
                           N_BVHD, -0.5, N_BVHD-0.5);
    h_bvh_ud[h]->SetDirectory(nullptr);
  }

  // 3) 1st pass: U×D 채우기 (H,U,D 모두 있는 이벤트만)
  const Long64_t N = tr->GetEntries();
  std::vector<int> hitsH, hitsU, hitsD;
  hitsH.reserve(8); hitsU.reserve(8); hitsD.reserve(8);

  Long64_t cnt_BH2_only = 0;
  Long64_t cnt_BH2_and_U = 0;
  Long64_t cnt_BH2_U_and_D = 0;

  for(Long64_t i=0;i<N;++i){
    if(i && (i%100000==0)) std::cout<<"  [1st] "<<i<<"/"<<N<<"\r"<<std::flush;
    tr->GetEntry(i);

    get_unique_hits(BH2 ,N_BH2 , ecutBH2MeV, hitsH);
    get_unique_hits(BVHU,N_BVHU, ecutUMeV  , hitsU);
    get_unique_hits(BVHD,N_BVHD, ecutDMeV  , hitsD);

    const bool hasH=!hitsH.empty(), hasU=!hitsU.empty(), hasD=!hitsD.empty();

    if(hasH && !hasU) cnt_BH2_only++;
    if(hasH &&  hasU) cnt_BH2_and_U++;
    if(hasH &&  hasU && hasD) cnt_BH2_U_and_D++;

    if(!(hasH && hasU && hasD)) continue;
    for(int h: hitsH){
      TH2F* H = h_bvh_ud[h];
      for(int u: hitsU) for(int d: hitsD) H->Fill(u,d);
    }
  }
  std::cout << "\n[Info] 1st pass finished.\n";

  // 4) 마스크 구축 (빈도 >= event_threshold)
  std::vector<std::vector<char>> mask(N_BH2, std::vector<char>(N_BVHU*N_BVHD, 0));
  for(int h=0; h<N_BH2; ++h){
    TH2F* H = h_bvh_ud[h];
    for(int bx=1; bx<=H->GetNbinsX(); ++bx){
      for(int by=1; by<=H->GetNbinsY(); ++by){
        if(H->GetBinContent(bx,by) >= event_threshold){
          int u = bx-1, d = by-1;
          mask[h][u*N_BVHD + d] = 1;
        }
      }
    }
  }

  // 5) 2nd pass: 두 가지 분모를 동시에 집계
  Long64_t denom_UH_all = 0, num_UH_all = 0;                 // 분모 A: BH2&&U
  Long64_t denom_H_all  = 0, num_H_all  = 0;                 // 분모 B: BH2 only
  std::vector<Long64_t> denom_UH_h(N_BH2,0), num_UH_h(N_BH2,0);
  std::vector<Long64_t> denom_H_h (N_BH2,0), num_H_h (N_BH2,0);

  for(Long64_t i=0;i<N;++i){
    if(i && (i%100000==0)) std::cout<<"  [2nd] "<<i<<"/"<<N<<"\r"<<std::flush;
    tr->GetEntry(i);

    get_unique_hits(BH2 ,N_BH2 , ecutBH2MeV, hitsH);
    get_unique_hits(BVHU,N_BVHU, ecutUMeV  , hitsU);
    get_unique_hits(BVHD,N_BVHD, ecutDMeV  , hitsD);

    const bool hasH=!hitsH.empty(), hasU=!hitsU.empty(), hasD=!hitsD.empty();

    // --- 분모 B: BH2 only ---
    if(hasH){
      denom_H_all++;
      for(int h: hitsH) denom_H_h[h]++;
    }

    // --- 분모 A: BH2 && U ---
    if(hasH && hasU){
      denom_UH_all++;
      bool global_match = false;

      for(int h: hitsH){
        bool matched_h = false;
        denom_UH_h[h]++;

        if(hasD){
          for(int u: hitsU){
            const int base = u * N_BVHD;
            for(int d: hitsD){
              if(mask[h][base + d]){ matched_h = true; global_match = true; break; }
            }
            if(matched_h) break;
          }
        }

        if(matched_h){
          num_UH_h[h]++;     // 분모 A의 분자
          num_H_h[h]++;      // 분모 B의 분자도 동일히 증가 (U&D&mask 매칭이면 BH2-only에서도 잘린 것)
        }
      }

      if(global_match){
        num_UH_all++;
        num_H_all++;
      }
    }
    // hasH && !hasU 인 이벤트는 분모B에만 포함(분자 불가능→U가 없으므로 mask 매칭 불가)
  }
  std::cout << "\n[Info] 2nd pass finished.\n";

  // 6) 출력
  auto pct = [](Long64_t a, Long64_t b)->double{ return (b>0)? 100.0*double(a)/double(b) : 0.0; };

  std::cout << "\n========== Summary (energy cuts kept, mask threshold="<<event_threshold<<") ==========\n";
  std::cout << std::fixed << std::setprecision(3);
  std::cout << "Input file                         : " << fname << "\n";
  std::cout << "Total entries                      : " << N << "\n";
  std::cout << "Cuts (MeV)    BH2="<<ecutBH2MeV<<"  U="<<ecutUMeV<<"  D="<<ecutDMeV << "\n";
  std::cout << "First-pass counts (for info):\n";
  std::cout << "  BH2 total (=BH2_only+BH2&&U)     : " << (cnt_BH2_only+cnt_BH2_and_U) << "\n";
  std::cout << "  BH2 && U                         : " << cnt_BH2_and_U << "\n";
  std::cout << "  BH2 && U && D                    : " << cnt_BH2_U_and_D << "\n";

  std::cout << "\n-- Denominator A:  (BH2 && BVH_U) --\n";
  std::cout << "  denom_UH_all = " << denom_UH_all << "\n";
  std::cout << "  num_UH_all   = " << num_UH_all   << "\n";
  std::cout << "  Veto rate A  = " << pct(num_UH_all, denom_UH_all) << " %\n";

  std::cout << "\n-- Denominator B:  (BH2 only) --\n";
  std::cout << "  denom_H_all  = " << denom_H_all << "\n";
  std::cout << "  num_H_all    = " << num_H_all   << "\n";
  std::cout << "  Veto rate B  = " << pct(num_H_all, denom_H_all) << " %\n";

  // 필요하면 BH2별 출력도 추가 가능 (주석 해제)
/*
  std::cout << "\nPer-BH2 (for Denominator A: BH2&&U)\n";
  std::cout << "  h :  num/den = rate(%)\n";
  for(int h=0; h<N_BH2; ++h){
    std::cout << " " << std::setw(2) << h << " : "
              << std::setw(9) << num_UH_h[h] << " / " << std::setw(9) << denom_UH_h[h]
              << " = " << std::setw(7) << pct(num_UH_h[h],denom_UH_h[h]) << "\n";
  }
  std::cout << "\nPer-BH2 (for Denominator B: BH2 only)\n";
  for(int h=0; h<N_BH2; ++h){
    std::cout << " " << std::setw(2) << h << " : "
              << std::setw(9) << num_H_h[h] << " / " << std::setw(9) << denom_H_h[h]
              << " = " << std::setw(7) << pct(num_H_h[h],denom_H_h[h]) << "\n";
  }
*/
}

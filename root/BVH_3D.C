// BVH_3D_english.C  (denominator = BH2 && BVH_U)
//
// - Reads branches:  BH2 (15 seg), BVH_U (22 seg), BVH_D (32 seg)
// - Energy cuts (MeV): BH2>=ecutBH2MeV, BVH_U>=ecutUMeV, BVH_D>=ecutDMeV
// - 1st pass: per-BH2 2D hist (U x D) 채우기  →  D가 실제로 히트한 경우만 U×D 셀 카운트
// - 마스크(mask[h][u,d])는 event_threshold 초과인 셀만 1로 표기
// - 2nd pass: **분모는 (BH2 && BVH_U) 이벤트**로 집계, 분자는 (hitsD 및 마스크 매칭)인 이벤트
//
// Summary 출력:
//   - Total events (TTree entries)
//   - BH2-only events (BH2 && !BVH_U)
//   - BH2 && BVH_U events  [= veto rate의 전역 분모]
//   - (참고) BH2 && BVH_U && BVH_D events
//
// 사용법:
//   .L BVH_3D_english.C+
//   run_analysis();  // 또는 BVH_3D("파일.root", 0.10, 0.04, 0.04, 1);

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
static const int N_BVHU = 15; // 0..21
static const int N_BVHD = 32; // 0..31

// 고유 세그먼트 집합 뽑기 (가중치=Edep[MeV]가 cut 이상인 것만)
static inline void get_unique_hits(const std::vector<TParticle>* v,
                                   int nmax, double cutMeV,
                                   std::vector<int>& out)
{
  out.clear();
  if(!v) return;
  std::unordered_set<int> s; s.reserve(8);
  for(const auto& p : *v){
    if(p.GetWeight() <= cutMeV) continue;
    int id = p.GetMother(1);
    if(0<=id && id<nmax) s.insert(id);
  }
  out.assign(s.begin(), s.end());
  std::sort(out.begin(), out.end());
}

void BVH_3D(const char* fname = "E45_BVH1_60mm.root",
            double ecutBH2MeV = 0.10,
            double ecutUMeV   = 0.04,
            double ecutDMeV   = 0.04,
            int    event_threshold = 1)
{
  gSystem->Load("libPhysics");
  gInterpreter->GenerateDictionary("vector<TParticle>", "TParticle.h;vector");

  // 1) 파일/트리
  TFile* f = TFile::Open(fname, "READ");
  if(!f || f->IsZombie()){ std::cerr << "[Error] Cannot open file: " << fname << std::endl; return; }
  TTree* tr = (TTree*)f->Get("g4hyptpc");
  if(!tr){ std::cerr << "[Error] Cannot find TTree 'g4hyptpc'." << std::endl; f->Close(); return; }

  std::vector<TParticle> *BH2=nullptr, *BVHU=nullptr, *BVHD=nullptr;
  tr->SetBranchAddress("BH2",   &BH2);
  tr->SetBranchAddress("BVH_U", &BVHU);
  tr->SetBranchAddress("BVH_D", &BVHD);

  // 2) BH2별 U×D 히스토그램
  std::vector<TH2F*> h_bvh_ud(N_BH2);
  for (int h = 0; h < N_BH2; ++h) {
    TString h_name  = TString::Format("h_bvh_ud_bh2_%d", h);
    TString h_title = TString::Format("BH2=%d;BVH_U Seg;BVH_D Seg", h);
    h_bvh_ud[h] = new TH2F(h_name, h_title,
                           N_BVHU, -0.5, N_BVHU-0.5,
                           N_BVHD, -0.5, N_BVHD-0.5);
  }

  // 3) 1st pass: U×D 빈도 채우기 (D 히트 존재할 때만)
  const Long64_t N = tr->GetEntries();
  std::vector<int> hitsH, hitsU, hitsD;
  hitsH.reserve(8); hitsU.reserve(8); hitsD.reserve(8);

  // 통계 집계용 카운터
  Long64_t cnt_BH2_only        = 0; // BH2 && !U
  Long64_t cnt_BH2_and_U       = 0; // BH2 && U   (veto rate 전역 분모)
  Long64_t cnt_BH2_U_and_D     = 0; // BH2 && U && D (참고)

  std::cout << "[Info] 1st pass: filling U×D per BH2, also counting categories..." << std::endl;
  for(Long64_t i = 0; i < N; ++i){
    if(i && (i % 100000 == 0)) std::cout << "  " << i << " / " << N << "\r" << std::flush;
    tr->GetEntry(i);

    get_unique_hits(BH2,  N_BH2,  ecutBH2MeV, hitsH);
    get_unique_hits(BVHU, N_BVHU, ecutUMeV,   hitsU);
    get_unique_hits(BVHD, N_BVHD, ecutDMeV,   hitsD);

    const bool hasH = !hitsH.empty();
    const bool hasU = !hitsU.empty();
    const bool hasD = !hitsD.empty();

    if(hasH && !hasU) cnt_BH2_only++;
    if(hasH &&  hasU) cnt_BH2_and_U++;
    if(hasH &&  hasU && hasD) cnt_BH2_U_and_D++;

    // U×D 2D 빈도는 실제 D가 있을 때만 채운다 (매트릭스는 U와 D의 동시 상관)
    if(!(hasH && hasU && hasD)) continue;

    for (int h : hitsH) {
      TH2F* H = h_bvh_ud[h];
      for (int u : hitsU) {
        for (int d : hitsD) {
          H->Fill(u, d);
        }
      }
    }
  }
  std::cout << "\n[Info] 1st pass finished." << std::endl;

  // 4) 마스크 구축(빈도 > event_threshold)
  std::vector< std::vector<char> > mask(N_BH2, std::vector<char>(N_BVHU * N_BVHD, 0));
  for(int h=0; h<N_BH2; ++h){
    TH2F* H = h_bvh_ud[h];
    for(int bx=1; bx<=H->GetNbinsX(); ++bx){
      for(int by=1; by<=H->GetNbinsY(); ++by){
        double c = H->GetBinContent(bx,by);
        if(c >= event_threshold){
          int u = bx-1;
          int d = by-1;
          mask[h][u*N_BVHD + d] = 1;
        }
      }
    }
  }

  // 5) 2nd pass: **분모 = BH2 && BVH_U**, 분자 = (마스크 매칭 & D히트 존재)
  Long64_t denom_all = 0, num_all = 0;
  std::vector<Long64_t> denom_h(N_BH2, 0), num_h(N_BH2, 0);

  tr->ResetBranchAddresses();
  tr->SetBranchAddress("BH2",   &BH2);
  tr->SetBranchAddress("BVH_U", &BVHU);
  tr->SetBranchAddress("BVH_D", &BVHD);

  std::cout << "[Info] 2nd pass: event-level veto counting (denom = BH2 && U)..." << std::endl;
  for(Long64_t i = 0; i < N; ++i){
    if(i && (i % 100000 == 0)) std::cout << "  " << i << " / " << N << "\r" << std::flush;
    tr->GetEntry(i);

    get_unique_hits(BH2,  N_BH2,  ecutBH2MeV, hitsH);
    get_unique_hits(BVHU, N_BVHU, ecutUMeV,   hitsU);
    get_unique_hits(BVHD, N_BVHD, ecutDMeV,   hitsD);

    const bool hasH = !hitsH.empty();
    const bool hasU = !hitsU.empty();
    const bool hasD = !hitsD.empty();

    // 분모: BH2 && U 가 히트한 이벤트만
    if(!(hasH && hasU)) continue;
    denom_all++;

    bool global_match = false;

    for(int h : hitsH){
      bool matched_h = false;
      denom_h[h]++;

      if(hasD){
        for(int u : hitsU){
          const int base = u * N_BVHD;
          for(int d : hitsD){
            if(mask[h][base + d]){
              matched_h = true;
              global_match = true;
              break;
            }
          }
          if(matched_h) break;
        }
      }

      if(matched_h) num_h[h]++;
    }

    if(global_match) num_all++;
  }
  std::cout << "\n[Info] 2nd pass finished." << std::endl;

  // 6) 그리기
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kViridis);
  TLine* grid_line = new TLine();
  grid_line->SetLineStyle(kDotted);
  grid_line->SetLineColor(kGray+1);

  std::vector<TCanvas*> canvases;
  for (int h = 0; h < N_BH2; ++h) {
    if (h % 3 == 0) {
      TString c_name  = TString::Format("c_bh2_%d_to_%d", h, std::min(h+2,N_BH2-1));
      TString c_title = TString::Format("BH2 Segments %d to %d", h, std::min(h+2,N_BH2-1));
      TCanvas* c = new TCanvas(c_name, c_title, 1500, 500);
      c->Divide(3, 1);
      canvases.push_back(c);
    }

    canvases.back()->cd((h % 3) + 1);
    TH2F* H = h_bvh_ud[h];
    H->GetXaxis()->SetTitle("BVH_U Seg");
    H->GetYaxis()->SetTitle("BVH_D Seg");
    H->Draw("COLZ");
    gPad->Update();

    const double x_min = H->GetXaxis()->GetXmin();
    const double x_max = H->GetXaxis()->GetXmax();
    const double y_min = H->GetYaxis()->GetXmin();
    const double y_max = H->GetYaxis()->GetXmax();

    for (int i = 0; i <= N_BVHU; ++i) {
      const double x_pos = H->GetXaxis()->GetBinLowEdge(i+1);
      grid_line->DrawLine(x_pos, y_min, x_pos, y_max);
    }
    for (int j = 0; j <= N_BVHD; ++j) {
      const double y_pos = H->GetYaxis()->GetBinLowEdge(j+1);
      grid_line->DrawLine(x_min, y_pos, x_max, y_pos);
    }

    for (int bx = 1; bx <= H->GetNbinsX(); ++bx) {
      for (int by = 1; by <= H->GetNbinsY(); ++by) {
        const double count = H->GetBinContent(bx, by);
        if (count >= event_threshold) {
          const double x1 = H->GetXaxis()->GetBinLowEdge(bx);
          const double x2 = H->GetXaxis()->GetBinUpEdge(bx);
          const double y1 = H->GetYaxis()->GetBinLowEdge(by);
          const double y2 = H->GetYaxis()->GetBinUpEdge(by);
          auto *box = new TBox(x1, y1, x2, y2);
          box->SetFillStyle(0);
          box->SetLineColor(kRed);
          box->SetLineWidth(2);
          box->Draw("SAME");
        }
      }
    }
    gPad->RedrawAxis();
  }

  std::cout << "\n========== Summary (denominator = BH2 && BVH_U) ==========\n";
  std::cout << std::fixed << std::setprecision(3);
  std::cout << "Input file                        : " << fname << "\n";
  std::cout << "Total events in TTree             : " << N << "\n";

  auto pct = [&](Long64_t x){ return (N>0) ? (100.0 * (double)x / (double)N) : 0.0; };

  Long64_t cnt_BH2_total = cnt_BH2_only + cnt_BH2_and_U;

  std::cout << "BH2 total events                  : " << cnt_BH2_total
            << " (" << pct(cnt_BH2_total) << "%)\n";
  std::cout << "BH2 && BVH_U events   [denom_all] : " << cnt_BH2_and_U
            << " (" << pct(cnt_BH2_and_U) << "%)\n";
  std::cout << "BH2 && U && D (for reference)     : " << cnt_BH2_U_and_D
            << " (" << pct(cnt_BH2_U_and_D) << "%)\n";

  std::cout << "Cuts (MeV):  BH2=" << ecutBH2MeV
            << "  U=" << ecutUMeV
            << "  D=" << ecutDMeV
            << "  | Highlight threshold(counts >=) " << event_threshold << "\n";

  std::cout << "Event-level counts (using denom = BH2&&U):\n";
  std::cout << "  denom_all = " << denom_all << "\n";
  std::cout << "  num_all   = " << num_all   << "\n";
  if(denom_all>0){
    double rate_all = (double)num_all / (double)denom_all * 100.0;
    std::cout << "Global veto rate (ANY BH2): " << num_all << " / " << denom_all
              << "  = " << rate_all << " %\n";
  }

  std::cout << "\nPer-BH2 segment veto rates (denom_h = events with BH2=h && U):\n";
  std::cout << "  h :  numerator / denominator  =  rate(%)\n";
  for(int h=0; h<N_BH2; ++h){
    if(denom_h[h]>0){
      double r = (double)num_h[h] / (double)denom_h[h] * 100.0;
      std::cout << " " << std::setw(2) << h << " :  "
                << std::setw(9) << num_h[h] << " / " << std::setw(9) << denom_h[h]
                << "  =  " << std::setw(7) << std::setprecision(3) << r << "\n";
    }else{
      std::cout << " " << std::setw(2) << h << " :  "
                << std::setw(9) << 0 << " / " << std::setw(9) << 0
                << "  =      n/a\n";
    }
  }
  std::cout << "===========================================================\n";
  std::cout << "[Info] Analysis & visualization done. "
            << canvases.size() << " canvas windows were created.\n";

  tr->ResetBranchAddresses();
}

void run_analysis() {
  int event_threshold = 1; // 빨간 박스 기준 (>= counts)
  BVH_3D("E45_BVH1_60mm.root", 0.10, 0.04, 0.04, event_threshold);
  std::cout << "[Config] Red boxes mark bins with counts >= " << event_threshold << ".\n";
}

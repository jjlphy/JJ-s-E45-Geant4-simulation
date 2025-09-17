// BVH_3D_english.C  (final, with English output & veto-rate calculation)
//
// - Reads branches:  BH2 (15 seg), BVH_U (22 seg), BVH_D (32 seg)
// - Energy cuts (MeV): BH2>=ecutBH2MeV, BVH_U>=ecutUMeV, BVH_D>=ecutDMeV
// - Fills 22x32 2D histograms (U x D) per BH2 segment
// - Draws dotted 22x32 grid and red boxes where bin content > event_threshold
// - Computes Veto rates with event-level numerators/denominators:
//      denom_all: #events with ≥1 BH2 hit (above BH2 threshold)
//      num_all  : #events that fall in ANY highlighted (U,D) bin for ANY hit BH2 seg
//      denom_h  : #events with a hit in BH2=h
//      num_h    : #events that fall in a highlighted (U,D) bin for that specific BH2=h
//
// Usage in ROOT:
//     .L BVH_3D_english.C+
//     run_analysis();                   // default
//     // or tweak thresholds:
//     BVH_3D("E45_BVH4.root", 0.10, 0.04, 0.04, 300);
//

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
#include <vector>
#include <unordered_set>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include "TInterpreter.h"
#include "TSystem.h"


// -------- Detector Segmentation --------
static const int N_BH2  = 15; // 0..14
static const int N_BVHU = 22; // 0..21
static const int N_BVHD = 40; // 0..31

// -------- Helper: unique segment indices passing e-dep cut --------
static inline void get_unique_hits(const std::vector<TParticle>* v,
                                   int nmax, double cutMeV,
                                   std::vector<int>& out)
{
  out.clear();
  if(!v) return;
  std::unordered_set<int> s;
  s.reserve(8);
  for(const auto& p : *v){
    if(p.GetWeight() <= cutMeV) continue;
    int id = p.GetMother(1);
    if(0<=id && id<nmax) s.insert(id);
  }
  out.assign(s.begin(), s.end());
  std::sort(out.begin(), out.end());
}

// ================== Main Analysis Function ==================
void BVH_3D(const char* fname = "E45_segment_40_420.root",
            double ecutBH2MeV = 0.10,
            double ecutUMeV   = 0.04,
            double ecutDMeV   = 0.04,
            int    event_threshold = 1)
{
    gSystem->Load("libPhysics");  // TParticle이 있는 라이브러리
gInterpreter->GenerateDictionary("vector<TParticle>", "TParticle.h;vector");

  // 1) Open file and TTree
  TFile* f = TFile::Open(fname, "READ");
  if(!f || f->IsZombie()){ std::cerr << "[Error] Cannot open file: " << fname << std::endl; return; }
  TTree* tr = (TTree*)f->Get("g4hyptpc");
  if(!tr){ std::cerr << "[Error] Cannot find TTree 'g4hyptpc'." << std::endl; f->Close(); return; }

  std::vector<TParticle> *BH2=nullptr, *BVHU=nullptr, *BVHD=nullptr;
  tr->SetBranchAddress("BH2",   &BH2);
  tr->SetBranchAddress("BVH_U", &BVHU);
  tr->SetBranchAddress("BVH_D", &BVHD); // fixed typo: "BVH_D"

  // 2) Create histograms (one per BH2 segment)
  std::vector<TH2F*> h_bvh_ud(N_BH2);
  for (int h = 0; h < N_BH2; ++h) {
    TString h_name  = TString::Format("h_bvh_ud_bh2_%d", h);
    TString h_title = TString::Format("BH2=%d;BVH_U Seg;BVH_D Seg", h);
    h_bvh_ud[h] = new TH2F(h_name, h_title,
                           N_BVHU, -0.5, N_BVHU-0.5,
                           N_BVHD, -0.5, N_BVHD-0.5);
  }

  // 3) First pass: fill histograms
  const Long64_t N = tr->GetEntries();
  std::vector<int> hitsH, hitsU, hitsD;
  hitsH.reserve(8); hitsU.reserve(8); hitsD.reserve(8);

  std::cout << "[Info] Processing " << N << " events..." << std::endl;
  for(Long64_t i = 0; i < N; ++i){
    if(i && (i % 100000 == 0)) std::cout << "  " << i << " / " << N << "\r" << std::flush;
    tr->GetEntry(i);

    get_unique_hits(BH2,  N_BH2,  ecutBH2MeV, hitsH);
    get_unique_hits(BVHU, N_BVHU, ecutUMeV,   hitsU);
    get_unique_hits(BVHD, N_BVHD, ecutDMeV,   hitsD);

    if (hitsH.empty() || hitsU.empty() || hitsD.empty()) continue;

    for (int h : hitsH) {
      TH2F* H = h_bvh_ud[h];
      for (int u : hitsU) {
        for (int d : hitsD) {
          H->Fill(u, d);
        }
      }
    }
  }
  std::cout << "\n[Info] Filling finished." << std::endl;

  // 4) Build highlighted masks for bins above threshold (per BH2)
  //    mask[h][u*N_BVHD + d] == true if that bin is highlighted
  std::vector< std::vector<char> > mask(N_BH2, std::vector<char>(N_BVHU * N_BVHD, 0));
  for(int h=0; h<N_BH2; ++h){
    TH2F* H = h_bvh_ud[h];
    for(int bx=1; bx<=H->GetNbinsX(); ++bx){
      for(int by=1; by<=H->GetNbinsY(); ++by){
        double c = H->GetBinContent(bx,by);
        if(c > event_threshold){
          int u = bx-1; // 0..N_BVHU-1
          int d = by-1; // 0..N_BVHD-1
          mask[h][u*N_BVHD + d] = 1;
        }
      }
    }
  }

  // 5) Second pass: event-level veto counting
  Long64_t denom_all = 0, num_all = 0;
  std::vector<Long64_t> denom_h(N_BH2, 0), num_h(N_BH2, 0);

  tr->ResetBranchAddresses();
  tr->SetBranchAddress("BH2",   &BH2);
  tr->SetBranchAddress("BVH_U", &BVHU);
  tr->SetBranchAddress("BVH_D", &BVHD);

  std::cout << "[Info] Second pass for event-level veto rates..." << std::endl;
  for(Long64_t i = 0; i < N; ++i){
    if(i && (i % 100000 == 0)) std::cout << "  " << i << " / " << N << "\r" << std::flush;
    tr->GetEntry(i);

    get_unique_hits(BH2,  N_BH2,  ecutBH2MeV, hitsH);
    get_unique_hits(BVHU, N_BVHU, ecutUMeV,   hitsU);
    get_unique_hits(BVHD, N_BVHD, ecutDMeV,   hitsD);

    if(hitsH.empty()) continue;
    denom_all++;

    // Global flag for "this event falls into any highlighted bin for any BH2 seg it hit"
    bool global_match = false;

    for(int h : hitsH){
      bool matched_h = false;
      denom_h[h]++;

      if(!hitsU.empty() && !hitsD.empty()){
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
  std::cout << "\n[Info] Event-level counting finished." << std::endl;

  // 6) Draw histograms with 22x32 grid and highlight boxes
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

    const double x_min = H->GetXaxis()->GetXmin(); // -0.5
    const double x_max = H->GetXaxis()->GetXmax(); // 21.5
    const double y_min = H->GetYaxis()->GetXmin(); // -0.5
    const double y_max = H->GetYaxis()->GetXmax(); // 31.5

    // vertical grid lines (U bins)
    for (int i = 0; i <= N_BVHU; ++i) {
      const double x_pos = H->GetXaxis()->GetBinLowEdge(i+1); // from -0.5 in steps of 1
      grid_line->DrawLine(x_pos, y_min, x_pos, y_max);
    }
    // horizontal grid lines (D bins)
    for (int j = 0; j <= N_BVHD; ++j) {
      const double y_pos = H->GetYaxis()->GetBinLowEdge(j+1);
      grid_line->DrawLine(x_min, y_pos, x_max, y_pos);
    }

    // red boxes for bins > threshold
    for (int bx = 1; bx <= H->GetNbinsX(); ++bx) {
      for (int by = 1; by <= H->GetNbinsY(); ++by) {
        const double count = H->GetBinContent(bx, by);
        if (count > event_threshold) {
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

  // 7) Print summary (English)
  std::cout << "\n========== Summary ==========\n";
  std::cout << std::fixed << std::setprecision(3);
  std::cout << "Input file                  : " << fname << "\n";
  std::cout << "BH2 energy cut (MeV)       : " << ecutBH2MeV << "\n";
  std::cout << "BVH_U energy cut (MeV)     : " << ecutUMeV   << "\n";
  std::cout << "BVH_D energy cut (MeV)     : " << ecutDMeV   << "\n";
  std::cout << "Highlight threshold (counts): " << event_threshold << "\n";
  std::cout << "Events with BH2 hit (denom): " << denom_all << "\n";
  if(denom_all>0){
    double rate_all = (double)num_all / (double)denom_all;
    std::cout << "Global veto rate (ANY BH2) : " << num_all << " / " << denom_all
              << "  = " << rate_all*100.0 << " %\n";
  }

  std::cout << "\nPer-BH2 segment veto rates:\n";
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
  std::cout << "=============================\n";
  std::cout << "[Info] Analysis & visualization done. " << canvases.size()
            << " canvas windows were created." << std::endl;

  // Keep tree/branches alive while canvases are open
  tr->ResetBranchAddresses();
}

// ================== Convenience wrapper ==================
void run_analysis() {
  int event_threshold = 1; // (1) change here to highlight bins >=300 counts
  BVH_3D("E45_segment_40_420.root", 0.10, 0.04, 0.04, event_threshold);
  std::cout << "[Config] Red boxes mark bins with counts >= " << event_threshold << "." << std::endl;
}

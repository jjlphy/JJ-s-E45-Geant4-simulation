// BVH_3D_ver3.C
// (exclusion/inclusion + global scale + denom=BH2&&BVH_U + palette tuning
//  + Dilation in U/SCH for veto matching
//  + Red boxes reflect dilation
//  + Save per-pad histogram PNG
//  + Print & SAVE all overlayed {BH2, BVH_U, BVH_D} triplets to overlay_triplets.txt)

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
#include "TGaxis.h"
#include "TPaletteAxis.h"
#include <vector>
#include <tuple>
#include <unordered_set>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>  // [NEW] for saving overlay triplets

// -------- Detector Segmentation --------
static const int N_BH2  = 15; // 0..14
static const int N_BVHU = 22; // 0..14
static const int N_SCH  = 64; // 0..63

// ====== Matching options ======
static const bool use_sch_dilation     = true; // SCH 겹침 보정 사용
static const int  sch_dilate_halfwidth = 1;    // ±1 strip 확장

// ====== U-side (BVH_U) dilation ======
static const bool use_u_dilation       = true; // BVH_U도 넉넉히 ±1 포함
static const int  u_dilate_halfwidth   = 1;    // ±1 strip 확장

// ====== Red-box drawing options ======
static const bool draw_dilated_boxes   = true; // 팽창된 이웃도 빨간 박스로 표시
static const bool distinguish_dilated  = true; // 원래셀=실선/굵게, 팽창셀=점선/얇게

// ====== Exclusion / Inclusion lists (BH2, BVH_U, BVH_D=SCH) ======
static std::vector<std::tuple<int,int,int>> kExcludedTriplets = { {0,1,31},{1,2,39},{2,8,1},{2,8,2}, {3,12,0}, {3,5,8}, {4,7,55}, {5,7,3}, {5,7, 7}, {5,8,8}, {6,9,60}, {6,9,63}, {6,9,4}, {6,9,60}, {6,9,63}, {7,11,62}, {7, 19, 31}, {8,11,58}, {8,16,10}, {8,16,11}, {8,10,12},{8,11,13},{8,12,10},{8,12,11}
,{9,8,54} ,{11,17,44}, {6,9,51}, {6, 9, 48} ,{7,11,49}, {7,10,50}, {4, 7, 4}, {4, 7, 7}, {5,5,42}
};
static std::vector<std::tuple<int,int,int>> kIncludedTriplets = {
  // {3, 4, 17},
};
static inline int key3(int h,int u,int s){ return (h<<16) | (u<<8) | s; }

// ===== Helper: branch auto-bind =====
static bool bind_branch(TTree* tr, const char* preferred,
                        const std::vector<const char*>& fallbacks,
                        std::vector<TParticle>*& ptr) {
  if (tr->GetBranch(preferred)) { tr->SetBranchAddress(preferred, &ptr); return true; }
  for (auto nm : fallbacks) {
    if (tr->GetBranch(nm)) { tr->SetBranchAddress(nm, &ptr);
      Warning("BVH_3D","using fallback branch '%s' for '%s'", nm, preferred);
      return true;
    }
  }
  Error("BVH_3D","branch '%s' not found (no fallback matched)", preferred);
  return false;
}

// ===== Helper: unique segment indices passing e-dep cut =====
static inline void get_unique_hits(const std::vector<TParticle>* v,
                                   int nmax, double cutMeV,
                                   std::vector<int>& out)
{
  out.clear();
  if(!v) return;
  std::unordered_set<int> s;
  s.reserve(8);
  for(const auto& p : *v){
    const double edepMeV = p.GetWeight();
    if(edepMeV < cutMeV) continue;
    int seg = p.GetMother(1);
    if(0<=seg && seg<nmax) s.insert(seg);
  }
  out.assign(s.begin(), s.end());
  std::sort(out.begin(), out.end());
}

// ================== Main Analysis Function ==================
void BVH_3D(const char* fname = "../E45_Beam_Ver3.root",
            double ecutBH2MeV   = 0.10,
            double ecutUMeV     = 0.04,
            double ecutSCHMeV   = 0.04,
            int    event_threshold = 1)
{
  gSystem->Load("libPhysics");
  gInterpreter->GenerateDictionary("vector<TParticle>", "TParticle.h;vector");

  // 0) Build exclusion/inclusion sets (fast lookup)
  std::unordered_set<int> EX, INC;
  EX.reserve(kExcludedTriplets.size()*2+8);
  INC.reserve(kIncludedTriplets.size()*2+8);
  for(const auto& t : kExcludedTriplets){
    int h,u,s; std::tie(h,u,s)=t;
    if(0<=h && h<N_BH2 && 0<=u && u<N_BVHU && 0<=s && s<N_SCH)
      EX.insert(key3(h,u,s));
  }
  for(const auto& t : kIncludedTriplets){
    int h,u,s; std::tie(h,u,s)=t;
    if(0<=h && h<N_BH2 && 0<=u && u<N_BVHU && 0<=s && s<N_SCH)
      INC.insert(key3(h,u,s));
  }

  // 1) Open file and TTree
  TFile* f = TFile::Open(fname, "READ");
  if(!f || f->IsZombie()){ std::cerr << "[Error] Cannot open file: " << fname << std::endl; return; }
  TTree* tr = (TTree*)f->Get("g4hyptpc");
  if(!tr){ std::cerr << "[Error] Cannot find TTree 'g4hyptpc'." << std::endl; f->Close(); return; }

  std::vector<TParticle> *BH2=nullptr, *BVHU=nullptr, *SCHv=nullptr;

  // 2) Branch binding with fallbacks
  bool ok = true;
  ok &= bind_branch(tr, "BH2",   {},                           BH2);
  ok &= bind_branch(tr, "BVH_U", {"BVH"},                      BVHU);
  ok &= bind_branch(tr, "SCH",   {"/SCH","SCH1","SCHHits"},    SCHv);
  if(!ok){ f->Close(); return; }

  // 3) Create histograms (one per BH2 segment)
  std::vector<TH2F*> h_ud(N_BH2);
  for (int h = 0; h < N_BH2; ++h) {
    TString h_name  = TString::Format("h_bvhu_sch_bh2_%d", h);
    TString h_title = TString::Format("BH2=%d;#BVH1 Seg;#BVH2 Seg");
    h_ud[h] = new TH2F(h_name, h_title,
                       N_BVHU, -0.5, N_BVHU-0.5,
                       N_SCH,  -0.5, N_SCH -0.5);
  }

  // 4) First pass: fill histograms (exclude during fill)
  const Long64_t N = tr->GetEntries();
  std::vector<int> hitsH, hitsU, hitsS;
  hitsH.reserve(8); hitsU.reserve(8); hitsS.reserve(8);

  std::cout << "[Info] Processing " << N << " events..." << std::endl;
  for(Long64_t i = 0; i < N; ++i){
    if(i && (i % 100000 == 0)) std::cout << "  " << i << " / " << N << "\r" << std::flush;
    tr->GetEntry(i);

    get_unique_hits(BH2,  N_BH2,  ecutBH2MeV, hitsH);
    get_unique_hits(BVHU, N_BVHU, ecutUMeV,   hitsU);
    get_unique_hits(SCHv, N_SCH,  ecutSCHMeV, hitsS);

    if (hitsH.empty() || hitsU.empty() || hitsS.empty()) continue;

    for (int h : hitsH) {
      TH2F* H = h_ud[h];
      for (int u : hitsU){
        for (int s : hitsS){
          if(EX.count(key3(h,u,s))) continue; // exclude는 fill 단계부터 배제
          H->Fill(u, s);
        }
      }
    }
  }
  std::cout << "\n[Info] Filling finished." << std::endl;

  // === Global color scale ===
  double globalMax = 0.0;
  for (int h = 0; h < N_BH2; ++h) {
    if (h_ud[h]) {
      const double m = h_ud[h]->GetMaximum();
      if (m > globalMax) globalMax = m;
    }
  }
  if (globalMax <= 0) globalMax = 1.0;
  std::cout << "[Info] Global Z max = " << globalMax << std::endl;

  // 5) Build masks with dilation (U±, SCH±). Exclusion wins over inclusion.
  std::vector< std::vector<char> > base_mask(N_BH2, std::vector<char>(N_BVHU * N_SCH, 0));
  std::vector< std::vector<char> > mask      (N_BH2, std::vector<char>(N_BVHU * N_SCH, 0));

  for(int h=0; h<N_BH2; ++h){
    TH2F* H = h_ud[h];
    for(int bx=1; bx<=H->GetNbinsX(); ++bx){
      for(int by=1; by<=H->GetNbinsY(); ++by){
        const int u = bx-1;
        const int s = by-1;
        const int k = key3(h,u,s);
        if(EX.count(k)) continue;  // 제외는 전부 패스

        const bool is_included = INC.count(k); // 강제 포함 여부
        const bool pass_thresh = (H->GetBinContent(bx,by) >= event_threshold);

        if(is_included || pass_thresh){
          base_mask[h][u*N_SCH + s] = 1;

          const int u_lo = use_u_dilation ? std::max(0, u - u_dilate_halfwidth) : u;
          const int u_hi = use_u_dilation ? std::min(N_BVHU-1, u + u_dilate_halfwidth) : u;
          const int s_lo = use_sch_dilation ? std::max(0, s - sch_dilate_halfwidth) : s;
          const int s_hi = use_sch_dilation ? std::min(N_SCH-1, s + sch_dilate_halfwidth) : s;

          for(int uu = u_lo; uu <= u_hi; ++uu){
            for(int ss = s_lo; ss <= s_hi; ++ss){
              const int kd = key3(h,uu,ss);
              if(!EX.count(kd))
                mask[h][uu*N_SCH + ss] = 1;
            }
          }
        }
      }
    }
  }

  // 6) Second pass: event-level veto counting (denominator = BH2 && BVH_U)
  Long64_t denom_all = 0, num_all = 0;
  std::vector<Long64_t> denom_h(N_BH2, 0), num_h(N_BH2, 0);

  tr->ResetBranchAddresses();
  bind_branch(tr, "BH2",   {},                        BH2);
  bind_branch(tr, "BVH_U", {"BVH"},                   BVHU);
  bind_branch(tr, "SCH",   {"/SCH","SCH1","SCHHits"}, SCHv);

  std::cout << "[Info] Second pass for event-level veto rates..." << std::endl;
  for(Long64_t i = 0; i < N; ++i){
    if(i && (i % 100000 == 0)) std::cout << "  " << i << " / " << N << "\r" << std::flush;
    tr->GetEntry(i);

    get_unique_hits(BH2,  N_BH2,  ecutBH2MeV, hitsH);
    get_unique_hits(BVHU, N_BVHU, ecutUMeV,   hitsU);
    get_unique_hits(SCHv, N_SCH,  ecutSCHMeV, hitsS);

    if(hitsH.empty() || hitsU.empty()) continue; // denom = BH2 && BVH_U
    denom_all++;

    bool global_match = false;
    for(int h : hitsH){
      bool matched_h = false;
      denom_h[h]++;

      if(!hitsS.empty()){
        for(int u : hitsU){
          const int base = u * N_SCH;
          for(int s : hitsS){
            const int k = key3(h,u,s);
            if(EX.count(k)) continue;           // 제외는 무시
            if(mask[h][base + s]) { matched_h = true; global_match = true; break; }
          }
          if(matched_h) break;
        }
      }
      if(matched_h) num_h[h]++;
    }
    if(global_match) num_all++;
  }
  std::cout << "\n[Info] Event-level counting finished." << std::endl;

  // 7) Global style for palette & axes
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kViridis);
  gStyle->SetNumberContours(50);
  TGaxis::SetMaxDigits(3);

  // === 오버레이 좌표 로깅 & 파일 저장 준비 ===
  static const bool log_overlay_coords   = true;                    // 터미널 출력
  static const bool save_overlay_to_file = true;                    // [NEW] 파일로 저장
  static const char* overlay_filename    = "overlay_triplets_ver3.txt";  // [NEW]

  std::unordered_set<int> overlay_set;
  overlay_set.reserve(4096);
  std::vector<std::tuple<int,int,int>> overlay_list;

  // 8) Draw histograms + grid + red boxes (mask 반영)
  TLine* grid_line = new TLine();
  grid_line->SetLineStyle(kDotted);
  grid_line->SetLineColor(kGray+1);

  std::vector<TCanvas*> canvases;
  for (int h = 0; h < N_BH2; ++h) {
    if (h % 3 == 0) {
      TString c_name  = TString::Format("c_bh2_%02d_to_%02d", h, std::min(h+2,N_BH2-1));
      TString c_title = TString::Format("BH2 Segments %d to %d", h, std::min(h+2,N_BH2-1));
      TCanvas* c = new TCanvas(c_name, c_title, 1500, 500);
      c->Divide(3, 1);
      canvases.push_back(c);
    }

    TCanvas* c = canvases.back(); c->cd((h % 3) + 1);
    gPad->SetRightMargin(0.16);

    TH2F* H = h_ud[h];
    H->SetMinimum(0);
    H->SetMaximum(globalMax);
    H->GetXaxis()->SetTitle("#BVH1 Seg");
    H->GetYaxis()->SetTitle("#BVH2 Seg");
    H->Draw("COLZ");
    gPad->Update();

    if (TPaletteAxis* pal = (TPaletteAxis*)H->GetListOfFunctions()->FindObject("palette")) {
      pal->SetX1NDC(pal->GetX1NDC() - 0.005);
      pal->SetX2NDC(pal->GetX2NDC() - 0.005);
      pal->SetLabelFont(42);
      pal->SetLabelSize(0.020);
      pal->SetTitleSize(0.020);
      pal->SetTitleOffset(0.9);
      H->GetZaxis()->SetLabelSize(0.020);
      H->GetZaxis()->SetTitleSize(0.020);
      H->GetZaxis()->SetTitleOffset(1.0);
    }

    const double x_min = H->GetXaxis()->GetXmin(), x_max = H->GetXaxis()->GetXmax();
    const double y_min = H->GetYaxis()->GetXmin(), y_max = H->GetYaxis()->GetXmax();

    for (int i = 0; i <= N_BVHU; ++i) {
      const double x_pos = H->GetXaxis()->GetBinLowEdge(i+1);
      grid_line->DrawLine(x_pos, y_min, x_pos, y_max);
    }
    for (int j = 0; j <= N_SCH; ++j) {
      const double y_pos = H->GetYaxis()->GetBinLowEdge(j+1);
      grid_line->DrawLine(x_min, y_pos, x_max, y_pos);
    }

    // 빨간 박스: mask(팽창 포함) 기준으로 표시 + 좌표 수집
    for (int bx = 1; bx <= H->GetNbinsX(); ++bx) {
      for (int by = 1; by <= H->GetNbinsY(); ++by) {
        const int u = bx-1, s = by-1;
        const int k = key3(h,u,s);
        if(EX.count(k)) continue; // 제외 조합은 표시/로그하지 않음

        const bool in_base = (base_mask[h][u*N_SCH + s] != 0);
        const bool in_mask = (mask[h][u*N_SCH + s]     != 0);

        if (!draw_dilated_boxes) {
          if (!in_base) continue;
        } else {
          if (!in_mask) continue;
        }

        // draw red box
        const double x1 = H->GetXaxis()->GetBinLowEdge(bx);
        const double x2 = H->GetXaxis()->GetBinUpEdge(bx);
        const double y1 = H->GetYaxis()->GetBinLowEdge(by);
        const double y2 = H->GetYaxis()->GetBinUpEdge(by);
        auto *box = new TBox(x1, y1, x2, y2);
        box->SetFillStyle(0);
        box->SetLineColor(kRed);
        if (draw_dilated_boxes && distinguish_dilated) {
          if (in_mask && !in_base) { box->SetLineStyle(kDashed); box->SetLineWidth(2); }
          else                     { box->SetLineStyle(kSolid);  box->SetLineWidth(3); }
        } else {
          box->SetLineStyle(kSolid);
          box->SetLineWidth(2);
        }
        box->Draw("SAME");

        // 수집(중복 제거)
        if (overlay_set.insert(k).second) {
          overlay_list.emplace_back(h,u,s);
        }
      }
    }
    gPad->RedrawAxis();

    // === 각 BH2 세그먼트별 pad만 PNG로 저장 ===
    {
      TString single = TString::Format("bvhu_sch_bh2_%02d.png", h);
      gPad->Update();
      gPad->SaveAs(single);
    }

    if (h % 3 == 2 || h == N_BH2-1) {
      TString out = TString::Format("bvhu_sch_bh2_%02d_to_%02d.png",
                                    h - (h%3), std::min(h - (h%3) + 2,N_BH2-1));
      canvases.back()->SaveAs(out);
    }
  }

  // === 정렬 후 터미널 출력 & 파일 저장 ===
  std::sort(overlay_list.begin(), overlay_list.end(),
            [](const std::tuple<int,int,int>& a, const std::tuple<int,int,int>& b){
              if (std::get<0>(a) != std::get<0>(b)) return std::get<0>(a) < std::get<0>(b);
              if (std::get<1>(a) != std::get<1>(b)) return std::get<1>(a) < std::get<1>(b);
              return std::get<2>(a) < std::get<2>(b);
            });

  if (log_overlay_coords) {
    std::cout << "\n[Overlay Triplets] (" << overlay_list.size() << " cells)\n";
    std::cout << "{BH2, BVH_U, BVH_D} list matching current overlay rule ("
              << (draw_dilated_boxes ? "mask+dilation" : "base only")
              << ", EX excluded):\n";
    for (const auto& t : overlay_list) {
      int h,u,s; std::tie(h,u,s)=t;
      std::cout << "{" << h << ", " << u << ", " << s << "}\n";
    }
    std::cout << "[End of Overlay Triplets]\n\n";
  }

  if (save_overlay_to_file) {
    std::ofstream ofs(overlay_filename);
    if (!ofs.is_open()) {
      Error("BVH_3D", "Failed to open '%s' for writing.", overlay_filename);
    } else {
      ofs << "# Overlay Triplets exported by BVH_3D_ver4.C\n";
      ofs << "# Format: {BH2, BVH_U, BVH_D}\n";
      ofs << "# Rule: " << (draw_dilated_boxes ? "mask+dilation" : "base only") << ", EX excluded\n";
      for (const auto& t : overlay_list) {
        int h,u,s; std::tie(h,u,s)=t;
        ofs << "{" << h << ", " << u << ", " << s << "}\n";
      }
      ofs.close();
      std::cout << "[Info] Saved overlay triplets to '" << overlay_filename
                << "' (" << overlay_list.size() << " lines)\n";
    }
  }

  // 9) Summary
  std::cout << "\n========== Summary ==========\n";
  std::cout << std::fixed << std::setprecision(3);
  std::cout << "Input file                   : " << fname << "\n";
  std::cout << "BH2 energy cut (MeV)        : " << ecutBH2MeV << "\n";
  std::cout << "BVH_U energy cut (MeV)      : " << ecutUMeV     << "\n";
  std::cout << "SCH energy cut (MeV)        : " << ecutSCHMeV   << "\n";
  std::cout << "Highlight threshold (counts): " << event_threshold << " (>=)\n";
  std::cout << "U dilation                  : " << (use_u_dilation ? "ON" : "OFF")
            << "  (±" << u_dilate_halfwidth << ")\n";
  std::cout << "SCH dilation                : " << (use_sch_dilation ? "ON" : "OFF")
            << "  (±" << sch_dilate_halfwidth << ")\n";
  std::cout << "Events with BH2 && BVH_U    : " << denom_all << "  (denominator)\n";
  if(denom_all>0){
    double rate_all = (double)num_all / (double)denom_all;
    std::cout << "Global veto rate (ANY BH2)  : " << num_all << " / " << denom_all
              << "  = " << rate_all*100.0 << " %\n";
  }
  std::cout << "\nPer-BH2 segment veto rates:\n";
  std::cout << "  h :  numerator / denominator  =  rate(%)\n";
  for(int h=0; h<N_BH2; ++h){
    if(denom_h[h]>0){
      double r = (double)num_h[h] / (double)denom_h[h] * 100.0;
      std::cout << " " << std::setw(2) << h << " :  "
                << std::setw(9) << num_h[h] << " / " << std::setw(9) << denom_h[h]
                << "  =  " << std::setw(7) << r << "\n";
    }else{
      std::cout << " " << std::setw(2) << h << " :  "
                << std::setw(9) << 0 << " / " << std::setw(9) << 0
                << "  =      n/a\n";
    }
  }
  std::cout << "=============================\n";
  std::cout << "[Info] Analysis & visualization done. "
            << canvases.size() << " canvas windows were created." << std::endl;

  tr->ResetBranchAddresses();
}

// ================== Convenience wrapper ==================
void run_analysis() {
  int event_threshold = 1; // 필요 시 100 등으로 상향
  BVH_3D("../E45_Beam_Ver3.root", 0.10, 0.04, 0.04, event_threshold);
  std::cout << "[Config] Red boxes mark bins with "
            << (draw_dilated_boxes ? "mask (dilated)" : "base (non-dilated)")
            << " cells >= threshold.  (U dilation="
            << (use_u_dilation ? "ON" : "OFF")
            << ", ±" << u_dilate_halfwidth
            << "; SCH dilation="
            << (use_sch_dilation ? "ON" : "OFF")
            << ", ±" << sch_dilate_halfwidth
            << "; distinguish_dilated="
            << (distinguish_dilated ? "ON" : "OFF") << ")\n";
}

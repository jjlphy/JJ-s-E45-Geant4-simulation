// VP_view.C  (Reader-based; overlay∩overlay intersections; extra #2/#3 baselines)
// usage:
//   root -l
//   .L VP_view.C+
//   VP_view("../E45_VP5.root");

#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TBox.h"
#include "TLatex.h"
#include "TString.h"
#include "TParticle.h"
#include <map>
#include <vector>
#include <iostream>
#include <cmath>

// =============================
// ========= CONFIG ============
// =============================

// (0) 출력용 정보
static const Long64_t kTotalEventsInfo = 1000000; // 표기용

// (1) VP z 위치 & 두께
struct VPZ { double z, dz; };
static std::map<int, VPZ> kVPZ = {
  {1, {-1000.0, 1.0}},
  {2, { -900.0, 1.0}},
  {6, {  900.0, 1.0}},
  {7, { 1200.0, 1.0}},
  {8, { 1500.0, 1.0}}
};

// (2) ★ VP별 히스토그램 범위 & bin 수 ★  { xmin, xmax, ymin, ymax, nxbins, nybins }
struct VPRange { double xmin, xmax, ymin, ymax; int nxbins, nybins; };
static std::map<int, VPRange> kVPRange = {
  {1, { -400,  200,  -200,  120,  440, 240 }},
  {2, { -400,  200,  -200,  120,  440, 240 }},
  {6, {  100,  800,  -200,  200,  700, 400 }},
  {7, {  150,  850,  -220,  220,  700, 440 }},
  {8, {  200,  900,  -240,  240,  700, 480 }},
};

// (3) ★ VP별 오버레이 (세그 개수/폭/높이/라벨 간격) ★
struct OverlaySpec { int nSeg; double segW, segH; int labelEvery; };
static std::map<int, OverlaySpec> kOverlay = {
  {1, {22, 10.0,  140.0, 5}},   // Upstream(BVH_U)
  {2, {22, 10.0,  140.0, 5}},
  {6, {32, 10.0, 140.0,10}},   // Downstream(BVH_D)
  {7, {32, 10.0, 140.0,10}},
  {8, {32, 10.0, 140.0,10}},
};

// (4) ★ VP별 오버레이 중심 ★
struct Center { double x, y; };
static std::map<int, Center> kCenter = {
  {1, {-20.0, 0.0}},
  {2, {-10.0, 0.0}},
  {6, {340.0, 0.0}},
  {7, {450.0, 0.0}},
  {8, {600.0, 0.0}},
};

// (5) BH2 히스토그램 범위
struct Range2D { double xmin, xmax, ymin, ymax; int nxbins, nybins; };
static Range2D kBH2Range = { -200, 150, -100, 100, 350, 200 };

// =============================
// ========= INTERNALS =========
// =============================

struct Box { double xL, xR, yL, yR; };
static Box MakeBox(const OverlaySpec& P, const Center& C) {
  const double totalW = P.nSeg * P.segW;
  return { C.x - 0.5*totalW, C.x + 0.5*totalW, C.y - 0.5*P.segH, C.y + 0.5*P.segH };
}

static TH2D* BookVPHist(int vp, const char* titlePrefix)
{
  const auto& R = kVPRange.at(vp);
  TString name  = Form("h_vp%d", vp);
  TString title = Form("%s (VP%d, z=%.0f mm);X [mm];Y [mm]", titlePrefix, vp, kVPZ.at(vp).z);
  return new TH2D(name, title, R.nxbins, R.xmin, R.xmax, R.nybins, R.ymin, R.ymax);
}

static void DrawOverlay(int vp)
{
  const auto& O = kOverlay.at(vp);
  const auto& C = kCenter.at(vp);

  const double totalW  = O.nSeg * O.segW;
  const double start_x = C.x - 0.5 * totalW;
  const int    center_idx = O.nSeg / 2;

  for (int i=0; i<O.nSeg; ++i) {
    const double x1 = start_x + i * O.segW;
    auto* seg = new TBox(x1, C.y - O.segH/2.0, x1 + O.segW, C.y + O.segH/2.0);
    seg->SetFillStyle(0);
    seg->SetLineColor(kRed);
    seg->SetLineWidth(1);
    seg->SetLineStyle(1);
    seg->Draw("SAME");

    const int seg_num = (i < center_idx) ? (i - center_idx) : (i - center_idx + 1);
    if (O.labelEvery>0 && seg_num!=0 && (std::abs(seg_num) % O.labelEvery == 0)) {
      auto* t = new TLatex(x1 + 0.5*O.segW, C.y + O.segH/2.0 + (std::isfinite(O.segH)? (O.segH*0.07) : 8.0),
                           Form("%+d", seg_num));
      t->SetTextAlign(22);
      t->SetTextColor(kRed);
      t->SetTextSize(0.03);
      t->Draw("SAME");
    }
  }
  // 중앙 0 라벨
  auto* tz = new TLatex(C.x, C.y + O.segH/2.0 + (std::isfinite(O.segH)? (O.segH*0.07) : 8.0), "0");
  tz->SetTextAlign(22);
  tz->SetTextColor(kRed);
  tz->SetTextSize(0.03);
  tz->Draw("SAME");
}

// =============================
// ==========  MAIN  ===========
// =============================
void VP_view(const char* filepath="../E45_VP5.root")
{
  // 파일/트리
  TFile* f = TFile::Open(filepath);
  if (!f || f->IsZombie()) { std::cerr << "[Error] cannot open: " << filepath << "\n"; return; }
  auto* tree = (TTree*)f->Get("g4hyptpc");
  if (!tree) { std::cerr << "[Error] TTree 'g4hyptpc' not found.\n"; return; }

  // Reader
  TTreeReader reader(tree);
  TTreeReaderValue<std::vector<TParticle>> VP(reader,  "VP");
  TTreeReaderValue<std::vector<TParticle>> BH2(reader, "BH2");

  // 히스토그램
  std::map<int, TH2D*> H;
  for (int vp : {1,2,6,7,8}) {
    const char* pre = (vp==1 || vp==2) ? "Upstream (BVH_U-like)" : "Downstream (BVH_D-like)";
    H[vp] = BookVPHist(vp, pre);
  }
  const auto& Rb = kBH2Range;
  TH2D* h_bh2 = new TH2D("h_bh2", "BH2 Hit Distribution;X [mm];Y [mm]",
                         Rb.nxbins, Rb.xmin, Rb.xmax, Rb.nybins, Rb.ymin, Rb.ymax);

  // 오버레이 박스 (VP별)
  std::map<int, Box> BoxMap;
  for (int vp : {1,2,6,7,8}) BoxMap[vp] = MakeBox(kOverlay.at(vp), kCenter.at(vp));

  // 카운트
  Long64_t nBH2 = 0;
  Long64_t N_in_overlay[9] = {0}; // 이벤트 단위 카운트 (1..8 사용)

  // overlay∩overlay 교차 (분모는 기본적으로 nBH2)
  Long64_t N_vp2_vp6 = 0, N_vp2_vp7 = 0, N_vp2_vp8 = 0;
  Long64_t N_vp1_vp6 = 0, N_vp1_vp7 = 0, N_vp1_vp8 = 0;

  // ★ 추가: 2) BH2 && VP1(overlay in) 분모, 3) BH2 && VP2(overlay in) 분모
  Long64_t nBH2_and_VP1 = 0;
  Long64_t nBH2_and_VP2 = 0;

  // 루프
  while (reader.Next()) {
    const auto& vps = *VP;
    const auto& b   = *BH2;

    const bool passBH2 = !b.empty();
    if (passBH2) ++nBH2;
    for (const auto& h : b) h_bh2->Fill(h.Vx(), h.Vy());

    // (이벤트 단위) 플래그
    bool inBoxVP1=false, inBoxVP2=false, inBoxVP6=false, inBoxVP7=false, inBoxVP8=false;

    // 히스토그램 채움용 (z-슬라이스)
    bool hitVP1=false, hitVP2=false, hitVP6=false, hitVP7=false, hitVP8=false;

    for (const auto& p : vps) {
      const double x = p.Vx(), y = p.Vy(), z = p.Vz();
      auto inZ   = [&](int vp)->bool { const auto& Z = kVPZ.at(vp); return std::abs(z - Z.z) < Z.dz; };
      auto inBox = [&](int vp)->bool { const auto& bx = BoxMap.at(vp); return (x>bx.xL && x<bx.xR && y>bx.yL && y<bx.yR); };

      // z-슬라이스에 들어오면 히스토그램 채움
      if (inZ(1)) { H[1]->Fill(x,y); hitVP1 = true; if (inBox(1)) inBoxVP1 = true; }
      if (inZ(2)) { H[2]->Fill(x,y); hitVP2 = true; if (inBox(2)) inBoxVP2 = true; }
      if (inZ(6)) { H[6]->Fill(x,y); hitVP6 = true; if (inBox(6)) inBoxVP6 = true; }
      if (inZ(7)) { H[7]->Fill(x,y); hitVP7 = true; if (inBox(7)) inBoxVP7 = true; }
      if (inZ(8)) { H[8]->Fill(x,y); hitVP8 = true; if (inBox(8)) inBoxVP8 = true; }
    }

    // 한 이벤트당 한 번만 overlay-in 카운트
    if (passBH2 && inBoxVP1) ++N_in_overlay[1];
    if (passBH2 && inBoxVP2) ++N_in_overlay[2];
    if (passBH2 && inBoxVP6) ++N_in_overlay[6];
    if (passBH2 && inBoxVP7) ++N_in_overlay[7];
    if (passBH2 && inBoxVP8) ++N_in_overlay[8];

    // 2),3)용 분모 카운트
    if (passBH2 && inBoxVP1) ++nBH2_and_VP1;
    if (passBH2 && inBoxVP2) ++nBH2_and_VP2;

    // overlay∩overlay 교차 (분모는 nBH2, 또는 아래 2)/3)에서 별도 분모 사용)
    if (passBH2) {
      if (inBoxVP2 && inBoxVP6) ++N_vp2_vp6;
      if (inBoxVP2 && inBoxVP7) ++N_vp2_vp7;
      if (inBoxVP2 && inBoxVP8) ++N_vp2_vp8;
      if (inBoxVP1 && inBoxVP6) ++N_vp1_vp6;
      if (inBoxVP1 && inBoxVP7) ++N_vp1_vp7;
      if (inBoxVP1 && inBoxVP8) ++N_vp1_vp8;
    }
  }

  auto pct = [&](Long64_t n, Long64_t den)->double { return (den>0) ? 100.0*double(n)/double(den) : 0.0; };

  // ===================== 출력 1) 기존 요약 =====================
  std::cout << "================ Summary (1) ==================\n";
  std::cout << "1) ROOT file           : " << filepath << "\n";
  std::cout << "2) Total events(info)  : " << kTotalEventsInfo << "\n";
  std::cout << "3) BH2 passing events  : " << nBH2 << "\n";
  std::cout << "-----------------------------------------------\n";
  std::cout << "Overlays (BH2-normalized)\n";
  std::cout << Form("   VP1 overlay in  : %lld (%.4f %%)\n", N_in_overlay[1], pct(N_in_overlay[1], nBH2));
  std::cout << Form("   VP2 overlay in  : %lld (%.4f %%)\n", N_in_overlay[2], pct(N_in_overlay[2], nBH2));
  std::cout << Form("   VP6 overlay in  : %lld (%.4f %%)\n", N_in_overlay[6], pct(N_in_overlay[6], nBH2));
  std::cout << Form("   VP7 overlay in  : %lld (%.4f %%)\n", N_in_overlay[7], pct(N_in_overlay[7], nBH2));
  std::cout << Form("   VP8 overlay in  : %lld (%.4f %%)\n", N_in_overlay[8], pct(N_in_overlay[8], nBH2));
  std::cout << "-----------------------------------------------\n";
  std::cout << "overlay∩overlay (BH2-normalized)\n";
  std::cout << Form("4) VP2&VP6 : N = %lld, frac = %.6f (%.4f %%)\n", N_vp2_vp6, (nBH2? double(N_vp2_vp6)/nBH2:0.0), pct(N_vp2_vp6, nBH2));
  std::cout << Form("5) VP2&VP7 : N = %lld, frac = %.6f (%.4f %%)\n", N_vp2_vp7, (nBH2? double(N_vp2_vp7)/nBH2:0.0), pct(N_vp2_vp7, nBH2));
  std::cout << Form("6) VP2&VP8 : N = %lld, frac = %.6f (%.4f %%)\n", N_vp2_vp8, (nBH2? double(N_vp2_vp8)/nBH2:0.0), pct(N_vp2_vp8, nBH2));
  std::cout << Form("7) VP1&VP6 : N = %lld, frac = %.6f (%.4f %%)\n", N_vp1_vp6, (nBH2? double(N_vp1_vp6)/nBH2:0.0), pct(N_vp1_vp6, nBH2));
  std::cout << Form("8) VP1&VP7 : N = %lld, frac = %.6f (%.4f %%)\n", N_vp1_vp7, (nBH2? double(N_vp1_vp7)/nBH2:0.0), pct(N_vp1_vp7, nBH2));
  std::cout << Form("9) VP1&VP8 : N = %lld, frac = %.6f (%.4f %%)\n", N_vp1_vp8, (nBH2? double(N_vp1_vp8)/nBH2:0.0), pct(N_vp1_vp8, nBH2));
  std::cout << "================================================\n\n";

  // ===================== 출력 2) BH2 && VP1 분모 =====================
  std::cout << "================ Summary (2)  BH2 && VP1 baseline ==================\n";
  std::cout << "BH2 && VP1 (overlay-in) events: " << nBH2_and_VP1 << "\n";
  // BH2 대비 감소 수/비율
  Long64_t loss_VP1 = (nBH2 > nBH2_and_VP1) ? (nBH2 - nBH2_and_VP1) : 0;
  std::cout << "Reduction vs BH2-only         : " << loss_VP1
            << " (" << pct(loss_VP1, nBH2) << " % of BH2)\n";
  // VP1 분모로 downstream 교차 재계산
  std::cout << "overlay∩overlay (normalized by BH2&&VP1)\n";
  std::cout << Form("VP1&VP6 : %lld (%.4f %%)\n", N_vp1_vp6, pct(N_vp1_vp6, nBH2_and_VP1));
  std::cout << Form("VP1&VP7 : %lld (%.4f %%)\n", N_vp1_vp7, pct(N_vp1_vp7, nBH2_and_VP1));
  std::cout << Form("VP1&VP8 : %lld (%.4f %%)\n", N_vp1_vp8, pct(N_vp1_vp8, nBH2_and_VP1));
  std::cout << "====================================================================\n\n";

  // ===================== 출력 3) BH2 && VP2 분모 =====================
  std::cout << "================ Summary (3)  BH2 && VP2 baseline ==================\n";
  std::cout << "BH2 && VP2 (overlay-in) events: " << nBH2_and_VP2 << "\n";
  Long64_t loss_VP2 = (nBH2 > nBH2_and_VP2) ? (nBH2 - nBH2_and_VP2) : 0;
  std::cout << "Reduction vs BH2-only         : " << loss_VP2
            << " (" << pct(loss_VP2, nBH2) << " % of BH2)\n";
  std::cout << "overlay∩overlay (normalized by BH2&&VP2)\n";
  std::cout << Form("VP2&VP6 : %lld (%.4f %%)\n", N_vp2_vp6, pct(N_vp2_vp6, nBH2_and_VP2));
  std::cout << Form("VP2&VP7 : %lld (%.4f %%)\n", N_vp2_vp7, pct(N_vp2_vp7, nBH2_and_VP2));
  std::cout << Form("VP2&VP8 : %lld (%.4f %%)\n", N_vp2_vp8, pct(N_vp2_vp8, nBH2_and_VP2));
  std::cout << "====================================================================\n";

  // ====== 그림 ======
  TCanvas* cU = new TCanvas("c_vpU", "VP1-2 (overlay)", 1200, 600);
  cU->Divide(2,1);
  { int i=1; for (int vp : {1,2}) { cU->cd(i++); H[vp]->Draw("COLZ"); DrawOverlay(vp); } }

  TCanvas* cD = new TCanvas("c_vpD", "VP6-8 (overlay)", 1200, 900);
  cD->Divide(2,2);
  { int pad=1; for (int vp : {6,7,8}) { cD->cd(pad++); H[vp]->Draw("COLZ"); DrawOverlay(vp); } }

  TCanvas* cBH2 = new TCanvas("c_bh2", "BH2 2D", 800, 600);
  h_bh2->Draw("COLZ");

  std::cout << "[Info] Done.\n";
}

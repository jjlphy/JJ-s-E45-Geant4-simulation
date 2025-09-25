// VP_view.C  (Reader-based robust version: no TTreeFormula; works for vector<TParticle>)
// usage:
//   root -l
//   .L VP_view.C+
//   VP_view("../E45_VP4.root");

#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TBox.h"
#include "TLine.h"
#include "TLatex.h"
#include "TString.h"
#include "TMath.h"
#include "TParticle.h"
#include <map>
#include <vector>
#include <unordered_set>
#include <iostream>

// ======================= [사용자 설정 파라미터] =======================
// (1) 히스토그램 범위(축)
static double kU_xmin = -200, kU_xmax =  200, kU_ymin = -100, kU_ymax =  100;  // VP1,2
static double kD_xmin =  100, kD_xmax =  800, kD_ymin = -200, kD_ymax =  200;  // VP6-8
static int    kNbinsXU = 400, kNbinsYU = 200;
static int    kNbinsXD = 700, kNbinsYD = 400;

// (BH2 2D 히스토그램 범위 — 필요시 조정)
static double kBH2_xmin = -100, kBH2_xmax = 200;
static double kBH2_ymin = -100, kBH2_ymax = 100;
static int    kNbinsXBH2 = 150, kNbinsYBH2 = 100;

// (2) z-평면 위치(파일에 맞춤): |Vz - z0| < dz
static std::map<int,double> kVPZ = {
  {1, -1000.0}, {2, -900.0},
  {6,   900.0}, {7,  1200.0}, {8,  1500.0}
};
static double kVPdz = 1.0; // mm

// (3) 오버레이(격자) 파라미터
struct OverlayParams {
  int    nSeg;        // 세그먼트 개수
  double segW;        // X-폭 [mm]
  double segH;        // Y-길이 [mm]
  int    labelEvery;  // 라벨 간격 (0=표시 안함)
  bool   showZero;    // 중앙(0) 라벨 표시
  int    lineColor;   // kRed 등
  int    lineWidth;   // 선 굵기
  int    lineStyle;   // 1=실선
  double textSize;    // 라벨 글자 크기
  double textYOffset; // 라벨 y 오프셋
};

// --- BVH_U (2×10×60 mm), 세그 15개 ---
static OverlayParams kU = {
  /*nSeg*/      15,
  /*segW*/      10.0,
  /*segH*/      60.0,
  /*labelEvery*/5,
  /*showZero*/  true,
  /*lineColor*/ kRed,
  /*lineWidth*/ 1,
  /*lineStyle*/ 1,
  /*textSize*/  0.03,
  /*textYOffset*/ 8.0
};

// --- BVH_D (2×10×140 mm), 세그 54개 ---
static OverlayParams kD = {
  /*nSeg*/      54,
  /*segW*/      10.0,
  /*segH*/      140.0,
  /*labelEvery*/10,
  /*showZero*/  true,
  /*lineColor*/ kRed,
  /*lineWidth*/ 1,
  /*lineStyle*/ 1,
  /*textSize*/  0.03,
  /*textYOffset*/ 12.0
};

// (4) ★ VP별 오버레이 중심 위치 ★
static std::map<int, double> kCenterX = {
  {1, -30.0}, {2, -30.0},
  {6, 340.0}, {7, 450.0}, {8, 600.0}
};
static std::map<int, double> kCenterY = {
  {1, 0.0}, {2, 0.0},
  {6, 0.0}, {7, 0.0}, {8, 0.0}
};

// (5) total events (정보 출력용)
static const Long64_t kTotalEvents = 1000000; // 요청값 표기

// ====================================================================
// 오버레이 그리기
static void draw_overlay(const OverlayParams& P, double centerX=0.0, double centerY=0.0)
{
  const double totalW  = P.nSeg * P.segW;
  const double start_x = centerX - 0.5 * totalW;
  const int    center_idx = P.nSeg / 2;

  for (int i=0; i<P.nSeg; ++i) {
    const double x1 = start_x + i * P.segW;
    TBox* seg = new TBox(x1, centerY - P.segH/2.0, x1 + P.segW, centerY + P.segH/2.0);
    seg->SetFillStyle(0);
    seg->SetLineColor(P.lineColor);
    seg->SetLineWidth(P.lineWidth);
    seg->SetLineStyle(P.lineStyle);
    seg->Draw("SAME");

    int seg_num = (i < center_idx) ? (i - center_idx) : (i - center_idx + 1);
    if (P.labelEvery>0 && seg_num!=0 && (std::abs(seg_num) % P.labelEvery == 0)) {
      TLatex* t = new TLatex(x1 + 0.5*P.segW, centerY + P.segH/2.0 + P.textYOffset,
                             Form("%+d", seg_num));
      t->SetTextAlign(22);
      t->SetTextColor(P.lineColor);
      t->SetTextSize(P.textSize);
      t->Draw("SAME");
    }
  }

  if (P.showZero) {
    TLatex* tz = new TLatex(centerX, centerY + P.segH/2.0 + P.textYOffset, "0");
    tz->SetTextAlign(22);
    tz->SetTextColor(P.lineColor);
    tz->SetTextSize(P.textSize);
    tz->Draw("SAME");
  }
}

static TH2D* book_vp_hist(int vp, const char* titlePrefix)
{
  const bool isU = (vp==1 || vp==2);
  TString name  = Form("h_vp%d", vp);
  TString title = Form("%s (VP%d, z=%.0f mm);X [mm];Y [mm]", titlePrefix, vp, kVPZ[vp]);
  if (isU) {
    return new TH2D(name, title, kNbinsXU, kU_xmin, kU_xmax,
                                kNbinsYU, kU_ymin, kU_ymax);
  } else {
    return new TH2D(name, title, kNbinsXD, kD_xmin, kD_xmax,
                                kNbinsYD, kD_ymin, kD_ymax);
  }
}

// 오버레이 박스 경계 계산
struct Box { double xL, xR, yL, yR; };
static Box MakeBox(const OverlayParams& P, double cx, double cy) {
  const double totalW = P.nSeg * P.segW;
  return { cx - 0.5*totalW, cx + 0.5*totalW, cy - 0.5*P.segH, cy + 0.5*P.segH };
}

// ------------------- [메인] -------------------
void VP_view(const char* filepath="../E45_VP4.root")
{
  // 파일/트리
  TFile* f = TFile::Open(filepath);
  if (!f || f->IsZombie()) {
    std::cerr << "[Error] cannot open file: " << filepath << std::endl;
    return;
  }
  TTree* tree = (TTree*)f->Get("g4hyptpc");
  if (!tree) {
    std::cerr << "[Error] TTree 'g4hyptpc' not found." << std::endl;
    return;
  }

  // Reader 준비
  TTreeReader reader(tree);
  TTreeReaderValue<std::vector<TParticle>> vp(reader, "VP");
  TTreeReaderValue<std::vector<TParticle>> bh2(reader, "BH2");

  // 히스토/카운트
  std::map<int, TH2D*> H;
  for (int vpID : {1,2,6,7,8}) {
    H[vpID] = book_vp_hist(vpID, (vpID==1||vpID==2) ? "Upstream (BVH_U-like)" : "Downstream (BVH_D-like)");
  }
  TH2D* h_bh2 = new TH2D("h_bh2", "BH2 Hit Distribution; X [mm]; Y [mm]",
                         kNbinsXBH2, kBH2_xmin, kBH2_xmax,
                         kNbinsYBH2, kBH2_ymin, kBH2_ymax);

  // 오버레이 박스 (VP별)
  std::map<int, Box> BoxMap = {
    {1, MakeBox(kU, kCenterX[1], kCenterY[1])},
    {2, MakeBox(kU, kCenterX[2], kCenterY[2])},
    {6, MakeBox(kD, kCenterX[6], kCenterY[6])},
    {7, MakeBox(kD, kCenterX[7], kCenterY[7])},
    {8, MakeBox(kD, kCenterX[8], kCenterY[8])},
  };

  // 카운트 변수
  Long64_t nTotal = tree->GetEntries();
  Long64_t nBH2   = 0;

  Long64_t N_in_overlay[9] = {0}; // index by vpID (1..8 used)
  // VP 교차 카운트 (분모는 BH2)
  Long64_t N_vp2_vp6 = 0, N_vp2_vp7 = 0, N_vp2_vp8 = 0;
  Long64_t N_vp1_vp6 = 0, N_vp1_vp7 = 0, N_vp1_vp8 = 0;

  // 메인 루프
  Long64_t iev = 0;
  while (reader.Next()) {
    ++iev;
    const auto& V = *vp;
    const auto& B = *bh2;

    bool passBH2 = !B.empty();
    if (passBH2) ++nBH2;

    // BH2 2D 히스토그램
    for (const auto& hit : B) {
      h_bh2->Fill(hit.Vx(), hit.Vy());
    }

    // 이 이벤트가 각 VP(z 슬라이스)에 히트했는지 flag
    bool hitVP1=false, hitVP2=false, hitVP6=false, hitVP7=false, hitVP8=false;

    // VP 히스토/오버레이 카운트
    for (const auto& p : V) {
      const double vx = p.Vx(), vy = p.Vy(), vz = p.Vz();

      auto inZ = [&](int vpID)->bool {
        return std::abs(vz - kVPZ[vpID]) < kVPdz;
      };
      auto inBox = [&](int vpID)->bool {
        const auto& bx = BoxMap[vpID];
        return (vx > bx.xL && vx < bx.xR && vy > bx.yL && vy < bx.yR);
      };

      // Upstream (VP1,2)
      if (inZ(1)) { H[1]->Fill(vx, vy); hitVP1 = true; if (passBH2 && inBox(1)) ++N_in_overlay[1]; }
      if (inZ(2)) { H[2]->Fill(vx, vy); hitVP2 = true; if (passBH2 && inBox(2)) ++N_in_overlay[2]; }
      // Downstream (VP6,7,8)
      if (inZ(6)) { H[6]->Fill(vx, vy); hitVP6 = true; if (passBH2 && inBox(6)) ++N_in_overlay[6]; }
      if (inZ(7)) { H[7]->Fill(vx, vy); hitVP7 = true; if (passBH2 && inBox(7)) ++N_in_overlay[7]; }
      if (inZ(8)) { H[8]->Fill(vx, vy); hitVP8 = true; if (passBH2 && inBox(8)) ++N_in_overlay[8]; }
    }

    // 교차 카운트 (분모는 BH2 이벤트)
    if (passBH2) {
      if (hitVP2 && hitVP6) ++N_vp2_vp6;
      if (hitVP2 && hitVP7) ++N_vp2_vp7;
      if (hitVP2 && hitVP8) ++N_vp2_vp8;
      if (hitVP1 && hitVP6) ++N_vp1_vp6;
      if (hitVP1 && hitVP7) ++N_vp1_vp7;
      if (hitVP1 && hitVP8) ++N_vp1_vp8;
    }
  }

  // 출력 (요청 포맷)
  auto pct = [&](Long64_t n)->double { return (nBH2>0) ? 100.0*double(n)/double(nBH2) : 0.0; };

  std::cout << "================ Summary ==================\n";
  std::cout << "1) ROOT file           : " << filepath << "\n";
  std::cout << "2) Total events        : " << kTotalEvents << " (info)\n";
  std::cout << "3) BH2 passing events  : " << nBH2 << "\n";
  std::cout << "-------------------------------------------\n";
  std::cout << "Overlays (BH2-normalized, for reference)\n";
  std::cout << Form("   VP1 overlay in  : %lld (%.4f %%)\n", N_in_overlay[1], pct(N_in_overlay[1]));
  std::cout << Form("   VP2 overlay in  : %lld (%.4f %%)\n", N_in_overlay[2], pct(N_in_overlay[2]));
  std::cout << Form("   VP6 overlay in  : %lld (%.4f %%)\n", N_in_overlay[6], pct(N_in_overlay[6]));
  std::cout << Form("   VP7 overlay in  : %lld (%.4f %%)\n", N_in_overlay[7], pct(N_in_overlay[7]));
  std::cout << Form("   VP8 overlay in  : %lld (%.4f %%)\n", N_in_overlay[8], pct(N_in_overlay[8]));
  std::cout << "-------------------------------------------\n";
  std::cout << "4) VP2(-900) & VP6(+900)    : N = " << N_vp2_vp6 << ",  frac = "
            << (nBH2>0 ? double(N_vp2_vp6)/nBH2 : 0.0) << " (" << pct(N_vp2_vp6) << " %)\n";
  std::cout << "5) VP2(-900) & VP7(+1200)   : N = " << N_vp2_vp7 << ",  frac = "
            << (nBH2>0 ? double(N_vp2_vp7)/nBH2 : 0.0) << " (" << pct(N_vp2_vp7) << " %)\n";
  std::cout << "6) VP2(-900) & VP8(+1500)   : N = " << N_vp2_vp8 << ",  frac = "
            << (nBH2>0 ? double(N_vp2_vp8)/nBH2 : 0.0) << " (" << pct(N_vp2_vp8) << " %)\n";
  std::cout << "7) VP1(-1000) & VP6(+900)   : N = " << N_vp1_vp6 << ",  frac = "
            << (nBH2>0 ? double(N_vp1_vp6)/nBH2 : 0.0) << " (" << pct(N_vp1_vp6) << " %)\n";
  std::cout << "8) VP1(-1000) & VP7(+1200)  : N = " << N_vp1_vp7 << ",  frac = "
            << (nBH2>0 ? double(N_vp1_vp7)/nBH2 : 0.0) << " (" << pct(N_vp1_vp7) << " %)\n";
  std::cout << "9) VP1(-1000) & VP8(+1500)  : N = " << N_vp1_vp8 << ",  frac = "
            << (nBH2>0 ? double(N_vp1_vp8)/nBH2 : 0.0) << " (" << pct(N_vp1_vp8) << " %)\n";
  std::cout << "===========================================\n";

  //  그림 그리기
  TCanvas* cU = new TCanvas("c_vpU", "VP1-2 (BVH_U overlay)", 1200, 600);
  cU->Divide(2,1);
  { int vpID=1; cU->cd(1); H[vpID]->Draw("COLZ"); draw_overlay(kU, kCenterX[vpID], kCenterY[vpID]); }
  { int vpID=2; cU->cd(2); H[vpID]->Draw("COLZ"); draw_overlay(kU, kCenterX[vpID], kCenterY[vpID]); }

  TCanvas* cD = new TCanvas("c_vpD", "VP6-8 (BVH_D overlay)", 1200, 900);
  cD->Divide(2,2);
  { int pad=1; for (int vpID: {6,7,8}) { cD->cd(pad++); H[vpID]->Draw("COLZ"); draw_overlay(kD, kCenterX[vpID], kCenterY[vpID]); } }

  TCanvas* cBH2 = new TCanvas("c_bh2", "BH2 2D", 800, 600);
  h_bh2->Draw("COLZ");

  std::cout << "[Info] Done. Reader-based workflow completed.\n";
}

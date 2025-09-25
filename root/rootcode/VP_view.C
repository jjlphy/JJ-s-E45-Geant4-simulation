// VP_view.C  (updated)
// 사용법:
//   root -l
//   .L VP_view.C+
//   VP_view("../E45_with_2pi.root");
//
// 기능 요약
//  - VP1,2 (Upstream): BVH_U 오버레이 (nSeg=15, 2×10×60 mm)
//  - VP6,7,8 (Downstream): BVH_D 오버레이 (nSeg=54, 2×10×140 mm)
//  - z 슬라이스: |VP.Vz() - z0| < kVPdz

#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TBox.h"
#include "TLine.h"
#include "TLatex.h"
#include "TString.h"
#include "TMath.h"
#include <map>
#include <utility>
#include <iostream>

// ======================= [사용자 설정 파라미터] =======================
// (1) 히스토그램 범위(축) — 필요 시 조정
static double kU_xmin = -300, kU_xmax =  300, kU_ymin = -200, kU_ymax =  200;  // VP1,2용
static double kD_xmin = -910, kD_xmax =  910, kD_ymin = -310, kD_ymax =  310;  // VP6-8용
static int    kNbinsXU = 200, kNbinsYU = 200;  // Upstream bin
static int    kNbinsXD = 400, kNbinsYD = 200;  // Downstream bin

// (2) z-평면 위치(단위 mm)와 선택 두께(|Vz - z0| < dz)  ← 실제 파일에 맞춤
static std::map<int,double> kVPZ = {
  {1, -1000.0}, {2, -900.0},
  {6,   900.0}, {7,  1200.0}, {8,  1500.0}
};
static double kVPdz = 1.0; // mm

// (3) 오버레이(격자) 파라미터
struct OverlayParams {
  int    nSeg;        // 세그먼트 개수
  double segW;        // 세그먼트 가로 폭 [mm]
  double segH;        // 세그먼트 세로 길이 [mm]
  int    labelEvery;  // 라벨 간격 (0이면 표시 안 함)
  bool   showZero;    // 중앙(0) 라벨 표시
  int    lineColor;   // kRed 등
  int    lineWidth;   // 선 굵기
  int    lineStyle;   // 1=실선, 2=점선 등
  double textSize;    // 라벨 글자 크기
  double textYOffset; // 라벨 y 오프셋
};

// --- BVH_U (2×10×60 mm), 세그 15개 ---
static OverlayParams kU = {
  /*nSeg*/      15,
  /*segW*/      10.0,    // X-폭
  /*segH*/      60.0,    // Y-길이(세그먼트 길이=60 mm)
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

// ====================================================================

static void draw_overlay(const OverlayParams& P)
{
  const double totalW  = P.nSeg * P.segW;
  const double start_x = -0.5 * totalW;
  const int center_idx = P.nSeg / 2;

  for (int i=0; i<P.nSeg; ++i) {
    const double x1 = start_x + i * P.segW;
    TBox* seg = new TBox(x1, -P.segH/2.0, x1 + P.segW, P.segH/2.0);
    seg->SetFillStyle(0);
    seg->SetLineColor(P.lineColor);
    seg->SetLineWidth(P.lineWidth);
    seg->SetLineStyle(P.lineStyle);
    seg->Draw("SAME");

    // 중앙 0은 비우고 양쪽 ±번호
    int seg_num = (i < center_idx) ? (i - center_idx) : (i - center_idx + 1);

    if (P.labelEvery>0 && seg_num!=0 && (std::abs(seg_num) % P.labelEvery == 0)) {
      TLatex* t = new TLatex(x1 + 0.5*P.segW, P.segH/2.0 + P.textYOffset,
                             Form("%+d", seg_num));
      t->SetTextAlign(22);
      t->SetTextColor(P.lineColor);
      t->SetTextSize(P.textSize);
      t->Draw("SAME");
    }
  }

  if (P.showZero) {
    TLatex* tz = new TLatex(0, P.segH/2.0 + P.textYOffset, "0");
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

void VP_view(const char* filepath="../E45_with_2pi.root")
{
  // 0) 파일/트리 열기
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

  // 1) 히스토그램 준비 (VP1,2,6,7,8)
  std::map<int, TH2D*> H;
  for (int vp : {1,2,6,7,8}) {
    H[vp] = book_vp_hist(vp, (vp==1||vp==2) ? "Upstream (BVH_U-like)" : "Downstream (BVH_D-like)");
  }

  // 2) 데이터 채우기: |Vz - z0| < kVPdz
  for (auto& kv : H) {
    int vp = kv.first;
    const double z0 = kVPZ[vp];
    TString sel  = Form("abs(VP.Vz() - (%.6f)) < %.6f", z0, kVPdz);
    TString targ = Form("VP.Vy():VP.Vx() >> %s", kv.second->GetName());
    tree->Draw(targ, sel, "goff");
  }

  // 3) 캔버스: U(1,2) 1창, D(6,7,8) 1창
  TCanvas* cU = new TCanvas("c_vpU", "VP1-2 (BVH_U overlay)", 1200, 600);
  cU->Divide(2,1);
  TCanvas* cD = new TCanvas("c_vpD", "VP6-8 (BVH_D overlay)", 1200, 900);
  cD->Divide(2,2); // 3장만 쓰고 1칸은 비게 둠

  // 4) 그리기 + 오버레이
  // Upstream
  for (int i=0; i<2; ++i) {
    int vp = (i==0 ? 1 : 2);
    cU->cd(i+1);
    H[vp]->Draw("COLZ");
    draw_overlay(kU);  // BVH_U 스타일 (15개, 60 mm)
  }

  // Downstream (VP6,7,8)
  int pad=1;
  for (int vp : {6,7,8}) {
    cD->cd(pad++);
    H[vp]->Draw("COLZ");
    draw_overlay(kD);  // BVH_D 스타일 (54개, 140 mm)
  }

  std::cout << "[Info] Done. You can tune overlay by editing kU/kD at the top." << std::endl;
}

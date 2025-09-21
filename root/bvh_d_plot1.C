// bvh_d_vetorate.C
// Usage:
//   root -l
//   .L bvh_d_vetorate.C+
//   draw_bvhd();

#include "TCanvas.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TStyle.h"

void draw_bvhd() {
  gStyle->SetOptStat(0);

  // -----------------------------
  // 데이터 입력
  // -----------------------------
  // 1) 세그먼트 개수 변화
  double seg_x[]    = {32, 34, 36, 38, 40, 50};
  double seg_rate[] = {97.986, 97.952, 98.027, 98.000, 98.072, 98.206};
  int nseg = sizeof(seg_x)/sizeof(seg_x[0]);

  // 2) 세그먼트 길이 변화
  double len_x[]    = {140, 200, 250};
  double len_rate[] = {97.986, 99.085, 99.223};
  int nlen = sizeof(len_x)/sizeof(len_x[0]);

  // -----------------------------
  // 캔버스 분할 (2개 그래프)
  // -----------------------------
  TCanvas *c1 = new TCanvas("c1","BVH_D optimization",900,400);
  c1->Divide(2,1);

  // -----------------------------
  // 그래프 1: 세그먼트 개수
  // -----------------------------
  c1->cd(1);
  auto gr1 = new TGraph(nseg, seg_x, seg_rate);
  gr1->SetTitle("BVH_D Segments vs Veto rate;# Segments;Veto rate [%]");
  gr1->SetLineColor(kBlue);
  gr1->SetLineWidth(2);
  gr1->SetMarkerColor(kBlue);
  gr1->SetMarkerStyle(20);
  gr1->SetMarkerSize(1.2);
  gr1->GetYaxis()->SetRangeUser(97.8, 98.3);
  gr1->Draw("ALP");

  // -----------------------------
  // 그래프 2: 세그먼트 길이
  // -----------------------------
  c1->cd(2);
  auto gr2 = new TGraph(nlen, len_x, len_rate);
  gr2->SetTitle("BVH_D Length vs Veto rate;Length [mm];Veto rate [%]");
  gr2->SetLineColor(kRed+1);
  gr2->SetLineWidth(2);
  gr2->SetMarkerColor(kRed+1);
  gr2->SetMarkerStyle(21);
  gr2->SetMarkerSize(1.2);
  gr2->GetYaxis()->SetRangeUser(97.8, 99.4);
  gr2->Draw("ALP");

  c1->Update();
}

#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TBox.h"
#include "TLine.h"
#include "TParticle.h"
#include "TStyle.h"
#include <vector>
#include <iostream>

//--- 각 검출기별 오버레이 그리기 함수 ---

void draw_bh2_overlay(double center_x) {
    const double seg_w = 14.0;
    const double seg_h = 100.0;
    const int n_segs = 15;
    const double total_w = n_segs * seg_w;
    const double start_x = center_x - total_w / 2.0;

    // 전체 외곽선
    TBox *outline = new TBox(start_x, -seg_h/2.0, start_x + total_w, seg_h/2.0);
    outline->SetFillStyle(0);
    outline->SetLineColor(kRed);
    outline->SetLineWidth(2);
    outline->Draw("SAME");

    // 세그먼트 구분선
    for (int i = 1; i < n_segs; ++i) {
        double x_pos = start_x + i * seg_w;
        TLine *line = new TLine(x_pos, -seg_h/2.0, x_pos, seg_h/2.0);
        line->SetLineColor(kRed);
        line->SetLineStyle(2);
        line->Draw("SAME");
    }
}

void draw_bvh_u_overlay(double center_x) {
    const double seg_w = 10.0;
    const double seg_h = 140.0; // BVH_U의 Y길이를 140mm로 가정
    const int n_segs = 26;
    const double total_w = n_segs * seg_w;
    const double start_x = center_x - total_w / 2.0;

    TBox *outline = new TBox(start_x, -seg_h/2.0, start_x + total_w, seg_h/2.0);
    outline->SetFillStyle(0);
    outline->SetLineColor(kRed);
    outline->SetLineWidth(2);
    outline->Draw("SAME");

    for (int i = 1; i < n_segs; ++i) {
        double x_pos = start_x + i * seg_w;
        TLine *line = new TLine(x_pos, -seg_h/2.0, x_pos, seg_h/2.0);
        line->SetLineColor(kRed);
        line->SetLineStyle(2);
        line->Draw("SAME");
    }
}

void draw_bvh_d_overlay(double center_x) {
    const double seg_w = 10.0;
    const double seg_h = 140.0;
    const int n_segs = 40;
    const double total_w = n_segs * seg_w;
    const double start_x = center_x - total_w / 2.0;

    TBox *outline = new TBox(start_x, -seg_h/2.0, start_x + total_w, seg_h/2.0);
    outline->SetFillStyle(0);
    outline->SetLineColor(kRed);
    outline->SetLineWidth(2);
    outline->Draw("SAME");

    for (int i = 1; i < n_segs; ++i) {
        double x_pos = start_x + i * seg_w;
        TLine *line = new TLine(x_pos, -seg_h/2.0, x_pos, seg_h/2.0);
        line->SetLineColor(kRed);
        line->SetLineStyle(2);
        line->Draw("SAME");
    }
}


// --- 메인 함수 ---
void view_detectors_BVH() {
    gStyle->SetOptStat(1111);
    gStyle->SetPalette(kBird);

    // 1. 파일 및 TTree 열기
    TFile *f = new TFile("E45_BVH2.root"); // 분석할 파일 이름
    if (!f || f->IsZombie()) { return; }
    TTree *tree = (TTree*)f->Get("g4hyptpc");
    if (!tree) { return; }

    // 2. 브랜치 읽기 설정
    std::vector<TParticle> *bh2_hits = nullptr;
    std::vector<TParticle> *bvh_u_hits = nullptr;
    std::vector<TParticle> *bvh_d_hits = nullptr;

    tree->SetBranchAddress("BH2", &bh2_hits);
    tree->SetBranchAddress("BVH_U", &bvh_u_hits);
    tree->SetBranchAddress("BVH_D", &bvh_d_hits);

    // 3. 히스토그램 정의 (이전 결과 기반으로 범위 설정)
    TH2D *h_hit_bh2   = new TH2D("h_hit_bh2",   "BH2 Hit Distribution;X [mm];Y [mm]", 250, -100, 150, 200, -100, 100);
    TH2D *h_hit_bvh_u = new TH2D("h_hit_bvh_u", "BVH_U Hit Distribution;X [mm];Y [mm]", 300, -150, 150, 200, -150, 150);
    TH2D *h_hit_bvh_d = new TH2D("h_hit_bvh_d", "BVH_D Hit Distribution;X [mm];Y [mm]", 500, 150, 650, 200, -150, 150);
    
    // 4. 이벤트 루프 및 히스토그램 채우기
    long long n_entries = tree->GetEntries();
    for (long long i = 0; i < n_entries; ++i) {
        tree->GetEntry(i);

        for (const auto& hit : *bh2_hits) {
            h_hit_bh2->Fill(hit.Vx(), hit.Vy());
        }
        for (const auto& hit : *bvh_u_hits) {
            h_hit_bvh_u->Fill(hit.Vx(), hit.Vy());
        }
        for (const auto& hit : *bvh_d_hits) {
            h_hit_bvh_d->Fill(hit.Vx(), hit.Vy());
        }
    }

    // 5. 시각화
    TCanvas *c_hits = new TCanvas("c_hits", "2D Hit Distributions with Overlays", 1800, 600);
    c_hits->Divide(3, 1);

    // BH2 그리기
    c_hits->cd(1);
    h_hit_bh2->Draw("COLZ");
    draw_bh2_overlay(35.0); // BH2 중심 X = 35.0 mm

    // BVH_U 그리기
    c_hits->cd(2);
    h_hit_bvh_u->Draw("COLZ");
    draw_bvh_u_overlay(30.0); // BVH_U 중심 X = 30.0 mm

    // BVH_D 그리기
    c_hits->cd(3);
    h_hit_bvh_d->Draw("COLZ");
    draw_bvh_d_overlay(400.0); // BVH_D 중심 X = 400.0 mm
}
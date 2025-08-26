#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TParticle.h>
#include <TPolyLine.h>
#include <TLine.h> // 세그먼트 선을 위해 포함
#include <TLegend.h>
#include <TMath.h>
#include <vector>
#include <iostream>
#include <string>

#pragma link C++ class vector<TParticle>+;

// ## 최종 버전: 범용 플롯 함수 ##
void plot_detector(TTree* tree, const std::string& detName, int nSegs, double segWidth, double segHeight, double rotation_deg = 0.0) {
    std::cout << "처리 중인 검출기: " << detName << "..." << std::endl;

    double totalWidth = nSegs * segWidth;
    double x_min = -totalWidth / 2.0;
    double x_max = totalWidth / 2.0;
    double y_min = -segHeight / 2.0;
    double y_max = segHeight / 2.0;

    double hist_x_min = x_min - 50, hist_x_max = x_max + 50;
    double hist_y_min = y_min - 50, hist_y_max = y_max + 50;
    if (rotation_deg != 0.0) {
        hist_x_min = -200; hist_x_max = 200;
        hist_y_min = -200; hist_y_max = 200;
    }

    std::string histName = "h_" + detName;
    std::string histTitle = detName + " Hit Distribution;X [mm];Y [mm];Hits";
    TH2F *h_hitmap = new TH2F(histName.c_str(), histTitle.c_str(), 400, hist_x_min, hist_x_max, 400, hist_y_min, hist_y_max);

    std::vector<TParticle> *hits = nullptr;
    tree->SetBranchAddress(detName.c_str(), &hits);

    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        if (hits) {
            for (const TParticle& particle : *hits) {
                h_hitmap->Fill(particle.Vx(), particle.Vy());
            }
        }
    }
    tree->ResetBranchAddresses();

    std::string canvasName = "c_" + detName;
    TCanvas *c = new TCanvas(canvasName.c_str(), (detName + " Hit Map").c_str(), 800, 800);
    h_hitmap->SetStats(0);
    h_hitmap->Draw("COLZ");

    // --- 외곽선 및 세그먼트 선 그리기 ---
    double angle_rad = rotation_deg * TMath::Pi() / 180.0;
    double cos_a = TMath::Cos(angle_rad);
    double sin_a = TMath::Sin(angle_rad);

    // 1. 외곽선 그리기
    double corners_x[5] = {x_min, x_max, x_max, x_min, x_min};
    double corners_y[5] = {y_min, y_min, y_max, y_max, y_min};
    double rotated_x[5], rotated_y[5];
    for (int i = 0; i < 5; ++i) {
        rotated_x[i] = corners_x[i] * cos_a - corners_y[i] * sin_a;
        rotated_y[i] = corners_x[i] * sin_a + corners_y[i] * cos_a;
    }
    TPolyLine *outline = new TPolyLine(5, rotated_x, rotated_y);
    outline->SetLineColor(kRed);
    outline->SetLineWidth(2);
    outline->SetFillStyle(0);
    outline->Draw("SAME");

    // 2. 내부 세그먼트 선 그리기
    for (int i = 1; i < nSegs; ++i) {
        // 로컬 좌표계 기준 세그먼트 선의 시작과 끝
        double line_x_local = x_min + i * segWidth;
        double line_y_start_local = y_min;
        double line_y_end_local = y_max;

        // 회전 변환 적용
        double x_start_rot = line_x_local * cos_a - line_y_start_local * sin_a;
        double y_start_rot = line_x_local * sin_a + line_y_start_local * cos_a;
        double x_end_rot = line_x_local * cos_a - line_y_end_local * sin_a;
        double y_end_rot = line_x_local * sin_a + line_y_end_local * cos_a;
        
        TLine* seg_line = new TLine(x_start_rot, y_start_rot, x_end_rot, y_end_rot);
        seg_line->SetLineColor(kRed);
        seg_line->SetLineStyle(2); // 점선으로 표시
        seg_line->Draw("SAME");
    }
}


// ## 메인 함수 ##
void plot_all_individually() {
    TFile *f = gFile;
    if (!f || f->IsZombie()) { /* ... */ return; }
    TTree *tree = (TTree*)f->Get("g4hyptpc");
    if (!tree) { /* ... */ return; }
    
    plot_detector(tree, "T0", 5, 32.0, 160.0, 45.0);
    plot_detector(tree, "BVH_U", 10, 10.0, 60.0);
    plot_detector(tree, "BH2", 12, 14.0, 100.0);
    plot_detector(tree, "BVH_D", 60, 10.0, 140.0);
}
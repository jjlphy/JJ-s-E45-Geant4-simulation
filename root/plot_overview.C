#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TParticle.h>
#include <TPolyLine.h>
#include <TMath.h>
#include <vector>
#include <iostream>
#include <string>

#pragma link C++ class vector<TParticle>+;

// ## 헬퍼 함수: 검출기 외곽선을 그리는 함수 ##
void draw_detector_outlines() {
    // --- T0 (기울어진 검출기) ---
    double t0_w = 5 * 32.0; // 160mm
    double t0_h = 160.0;
    double angle_rad = 45.0 * TMath::Pi() / 180.0;
    double cos_a = TMath::Cos(angle_rad);
    double sin_a = TMath::Sin(angle_rad);

    double corners_x[5] = {-t0_w/2, t0_w/2, t0_w/2, -t0_w/2, -t0_w/2};
    double corners_y[5] = {-t0_h/2, -t0_h/2, t0_h/2, t0_h/2, -t0_h/2};
    
    double rotated_corners_x[5], rotated_corners_y[5];
    for (int i = 0; i < 5; ++i) {
        rotated_corners_x[i] = corners_x[i] * cos_a - corners_y[i] * sin_a;
        rotated_corners_y[i] = corners_x[i] * sin_a + corners_y[i] * cos_a;
    }
    
    TPolyLine *t0_outline = new TPolyLine(5, rotated_corners_x, rotated_corners_y);
    t0_outline->SetLineColor(kGreen);
    t0_outline->SetLineWidth(2);
    t0_outline->SetFillStyle(0);
    t0_outline->Draw("SAME");

    // --- 다른 검출기들 ---
    TPolyLine *bvhu_outline = new TPolyLine(5, new double[5]{-50, 50, 50, -50, -50}, new double[5]{-30, -30, 30, 30, -30});
    bvhu_outline->SetLineColor(kCyan); bvhu_outline->SetLineWidth(2); bvhu_outline->SetFillStyle(0);
    bvhu_outline->Draw("SAME");
    
    TPolyLine *bh2_outline = new TPolyLine(5, new double[5]{-84, 84, 84, -84, -84}, new double[5]{-50, -50, 50, 50, -50});
    bh2_outline->SetLineColor(kMagenta); bh2_outline->SetLineWidth(2); bh2_outline->SetFillStyle(0);
    bh2_outline->Draw("SAME");

    TPolyLine *bvhd_outline = new TPolyLine(5, new double[5]{-300, 300, 300, -300, -300}, new double[5]{-70, -70, 70, 70, -70});
    bvhd_outline->SetLineColor(kRed); bvhd_outline->SetLineWidth(2); bvhd_outline->SetFillStyle(0);
    bvhd_outline->Draw("SAME");
}

// ## 메인 함수 ##
void plot_overview() {
    TFile *f = gFile;
    if (!f || f->IsZombie()) {
        std::cout << "오류: ROOT 파일이 열려있지 않습니다." << std::endl;
        return;
    }
    
    // ## 수정된 부분: TTree 이름을 g4hyptpc로 변경 ##
    TTree *tree = (TTree*)f->Get("g4hyptpc");
    if (!tree) {
        std::cout << "오류: TTree 'g4hyptpc'를 찾을 수 없습니다." << std::endl;
        return;
    }

    TH2F *h_overview = new TH2F("h_overview", "Overall Hit Distribution;X [mm];Y [mm];Hits",
                                800, -400, 400,
                                300, -150, 150);

    std::vector<TParticle> *t0_hits = nullptr, *bvhu_hits = nullptr, *bh2_hits = nullptr, *bvhd_hits = nullptr;

    tree->SetBranchAddress("T0", &t0_hits);
    tree->SetBranchAddress("BVH_U", &bvhu_hits);
    tree->SetBranchAddress("BH2", &bh2_hits);
    tree->SetBranchAddress("BVH_D", &bvhd_hits);

    Long64_t nEntries = tree->GetEntries();
    std::cout << "총 " << nEntries << "개의 이벤트를 처리합니다..." << std::endl;
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        if (t0_hits)   { for (const auto& p : *t0_hits)   h_overview->Fill(p.Vx(), p.Vy()); }
        if (bvhu_hits) { for (const auto& p : *bvhu_hits) h_overview->Fill(p.Vx(), p.Vy()); }
        if (bh2_hits)  { for (const auto& p : *bh2_hits)  h_overview->Fill(p.Vx(), p.Vy()); }
        if (bvhd_hits) { for (const auto& p : *bvhd_hits) h_overview->Fill(p.Vx(), p.Vy()); }
    }
    std::cout << "이벤트 처리 완료!" << std::endl;
    
    TCanvas *c_overview = new TCanvas("c_overview", "Overall Detector Hit Map", 1000, 700);
    h_overview->SetStats(0);
    h_overview->Draw("COLZ");

    draw_detector_outlines();
}
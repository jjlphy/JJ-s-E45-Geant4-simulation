#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TParticle.h>
#include <TLine.h>
#include <vector>
#include <iostream>

// ## 추가된 부분 ##
// ROOT 컴파일러에게 vector<TParticle>의 Dictionary를 강제로 생성하도록 지시합니다.
#pragma link C++ class vector<TParticle>+;

void plot_bvh_d() {
    TFile *f = gFile;
    if (!f || f->IsZombie()) {
        std::cout << "오류: ROOT 파일이 열려있지 않습니다." << std::endl;
        return;
    }

    TTree *tree = (TTree*)f->Get("g4hyptpc");
    if (!tree) {
        std::cout << "오류: TTree 'g4hyptpc'를 찾을 수 없습니다." << std::endl;
        return;
    }

    TH2F *h_hitmap_xy = new TH2F("h_hitmap_xy", "BVH_D Hit Distribution;X [mm];Y [mm];Hits",
                                600, -300.0, 300.0,
                                140, -70.0, 70.0);

    std::vector<TParticle> *bvh_hits = nullptr;
    // SetBranchAddress 호출 전에 경고가 발생하므로, 이 부분이 중요합니다.
    tree->SetBranchAddress("BVH_D", &bvh_hits);

    Long64_t nEntries = tree->GetEntries();
    std::cout << "총 " << nEntries << "개의 이벤트를 처리합니다..." << std::endl;

    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);

        // bvh_hits 포인터가 유효한지 매번 확인 (안전장치)
        if (bvh_hits) {
            for (const TParticle& particle : *bvh_hits) {
                double hit_x = particle.Vx();
                double hit_y = particle.Vy();
                h_hitmap_xy->Fill(hit_x, hit_y);
            }
        }
    }
    std::cout << "이벤트 처리 완료!" << std::endl;

    TCanvas *c1 = new TCanvas("c1", "BVH_D Physical Hit Map", 900, 500);
    c1->cd();
    h_hitmap_xy->SetStats(0);
    h_hitmap_xy->Draw("COLZ");

    for (int i = 0; i <= 60; ++i) {
        double x_pos = -300.0 + (i * 10.0);
        TLine *line = new TLine(x_pos, -70.0, x_pos, 70.0);
        line->SetLineColor(kRed);
        line->SetLineStyle(1);
        line->Draw("SAME");
    }
    
    TLine* box_top = new TLine(-300.0, 70.0, 300.0, 70.0);
    TLine* box_bottom = new TLine(-300.0, -70.0, 300.0, -70.0);
    box_top->SetLineColor(kRed);
    box_bottom->SetLineColor(kRed);
    box_top->Draw("SAME");
}

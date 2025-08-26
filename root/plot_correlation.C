#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TParticle.h>
#include <vector>
#include <iostream>
#include <string>

#pragma link C++ class vector<TParticle>+;

void plot_correlation() {
    TFile *f = gFile;
    if (!f || f->IsZombie()) { /* ... */ return; }
    TTree *tree = (TTree*)f->Get("g4hyptpc");
    if (!tree) { /* ... */ return; }

    // 3개의 2D 히스토그램 생성
    TH2F *h_T0_vs_BH2 = new TH2F("h_T0_vs_BH2", "T0 vs BH2 Hit Correlation;BH2 Segment ID;T0 Segment ID", 12, 0.5, 12.5, 5, 0.5, 5.5);
    TH2F *h_BVHU_vs_BH2 = new TH2F("h_BVHU_vs_BH2", "BVH_U vs BH2 Hit Correlation;BH2 Segment ID;BVH_U Segment ID", 12, 0.5, 12.5, 10, 0.5, 10.5);
    TH2F *h_BVHU_vs_BVHD = new TH2F("h_BVHU_vs_BVHD", "BVH_U vs BVH_D Hit Correlation;BVH_D Segment ID;BVH_U Segment ID", 60, 0.5, 60.5, 10, 0.5, 10.5);


    // 필요한 모든 브랜치에 변수 연결
    std::vector<TParticle> *t0_hits = nullptr;
    std::vector<TParticle> *bh2_hits = nullptr;
    std::vector<TParticle> *bvhu_hits = nullptr;
    std::vector<TParticle> *bvhd_hits = nullptr;

    tree->SetBranchAddress("T0", &t0_hits);
    tree->SetBranchAddress("BH2", &bh2_hits);
    tree->SetBranchAddress("BVH_U", &bvhu_hits);
    tree->SetBranchAddress("BVH_D", &bvhd_hits);

    Long64_t nEntries = tree->GetEntries();
    std::cout << "총 " << nEntries << "개의 이벤트를 분석합니다..." << std::endl;

    for (long long i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);

        // T0 vs BH2
        if (t0_hits && t0_hits->size() > 0 && bh2_hits && bh2_hits->size() > 0) {
            h_T0_vs_BH2->Fill(bh2_hits->at(0).GetMother(1) + 1, t0_hits->at(0).GetMother(1) + 1);
        }

        // BVH_U vs BH2
        if (bvhu_hits && bvhu_hits->size() > 0 && bh2_hits && bh2_hits->size() > 0) {
            h_BVHU_vs_BH2->Fill(bh2_hits->at(0).GetMother(1) + 1, bvhu_hits->at(0).GetMother(1) + 1);
        }

        // ## 추가된 부분: BVH_U vs BVH_D ##
        if (bvhu_hits && bvhu_hits->size() > 0 && bvhd_hits && bvhd_hits->size() > 0) {
            h_BVHU_vs_BVHD->Fill(bvhd_hits->at(0).GetMother(1) + 1, bvhu_hits->at(0).GetMother(1) + 1);
        }
    }
    std::cout << "분석 완료!" << std::endl;

    // 결과 그리기
    TCanvas *c1 = new TCanvas("c1", "T0 vs BH2 Correlation", 700, 500);
    h_T0_vs_BH2->Draw("COLZ");

    TCanvas *c2 = new TCanvas("c2", "BVH_U vs BH2 Correlation", 700, 500);
    h_BVHU_vs_BH2->Draw("COLZ");

    TCanvas *c3 = new TCanvas("c3", "BVH_U vs BVH_D Correlation", 900, 500);
    h_BVHU_vs_BVHD->Draw("COLZ");
}
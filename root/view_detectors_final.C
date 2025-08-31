// view_detectors.C (라이브러리 로드 최종본)
#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TParticle.h"
#include "TSystem.h" // gSystem을 사용하기 위해 필요
#include <vector>

void view_detectors_final() {
    // [가장 중요한 수정]
    // TParticle 클래스의 설명서(Dictionary)가 포함된 라이브러리를 명시적으로 로드합니다.
    gSystem->Load("libPhysics");

    TFile *f = new TFile("E45_test_VP.root");
    if (!f || f->IsZombie()) { return; }
    TTree *tree = (TTree*)f->Get("g4hyptpc");
    if (!tree) { return; }

    // Branch와 연결할 C++ 변수(포인터)를 선언합니다.
    std::vector<TParticle> *vp_hits = nullptr;
    std::vector<TParticle> *t0_hits = nullptr;
    std::vector<TParticle> *bh2_hits = nullptr;
    tree->SetBranchAddress("VP", &vp_hits);
    tree->SetBranchAddress("T0", &t0_hits);
    tree->SetBranchAddress("BH2", &bh2_hits);

    // 히스토그램들을 생성합니다.
    TH2D *h_vp4 = new TH2D("h_vp4", "VP4 at z=-760 mm; X [mm]; Y [mm]", 200, -150, 150, 200, -150, 150);
    TH2D *h_vp5 = new TH2D("h_vp5", "VP5 at z=470 mm; X [mm]; Y [mm]", 200, -400, 400, 200, -400, 400);
    TH2D *h_vp6 = new TH2D("h_vp6", "VP6 at z=490 mm; X [mm]; Y [mm]", 200, -400, 400, 200, -400, 400);
    TH2D *h_t0  = new TH2D("h_t0",  "T0 Hit Distribution; X [mm]; Y [mm]", 200, -150, 150, 200, -150, 150);
    TH2D *h_bh2 = new TH2D("h_bh2", "BH2 Hit Distribution; X [mm]; Y [mm]", 200, -150, 150, 200, -150, 150);

    // 이벤트 루프를 돌면서 히스토그램을 직접 채웁니다.
    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; i++) {
        tree->GetEntry(i);
        
        if (vp_hits) { // 포인터가 유효한지 확인
            for (const auto& hit : *vp_hits) {
                double z = hit.Vz();
                if (abs(z - (-760.0)) < 1.0) h_vp4->Fill(hit.Vx(), hit.Vy());
                else if (abs(z - 470.0) < 1.0) h_vp5->Fill(hit.Vx(), hit.Vy());
                else if (abs(z - 490.0) < 1.0) h_vp6->Fill(hit.Vx(), hit.Vy());
            }
        }
        if (t0_hits) {
            for (const auto& hit : *t0_hits) h_t0->Fill(hit.Vx(), hit.Vy());
        }
        if (bh2_hits) {
            for (const auto& hit : *bh2_hits) h_bh2->Fill(hit.Vx(), hit.Vy());
        }
    }

    // 캔버스를 생성하고 결과를 그립니다.
    TCanvas *c1 = new TCanvas("c1", "Detector Hit Distributions", 1200, 1800);
    c1->Divide(2, 3);
    
    c1->cd(1); h_vp4->Draw("COLZ");
    c1->cd(2); h_bh2->Draw("COLZ");
    c1->cd(3); h_vp5->Draw("COLZ");
    c1->cd(4); h_vp6->Draw("COLZ");
    c1->cd(5); h_t0->Draw("COLZ");
}
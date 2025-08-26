#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TParticle.h>
#include <TLine.h>
#include <vector>
#include <iostream>
#include <string>

// vector<TParticle>의 Dictionary를 강제로 생성하도록 지시합니다.
#pragma link C++ class vector<TParticle>+;

// ## 핵심: 재사용 가능한 히스토그램 생성 함수 ##d
// 어떤 검출기든 이름, 세그먼트 수, 크기 정보만 주면 그려줍니다.
void plot_detector_hitmap(TTree* tree, const std::string& detName, int nSegs, double segWidth, double segHeight) {
    std::cout << "처리 중인 검출기: " << detName << "..." << std::endl;

    // 1. 검출기 물리적 크기 계산
    double totalWidth = nSegs * segWidth;
    double x_min = -totalWidth / 2.0;
    double x_max = totalWidth / 2.0;
    double y_min = -segHeight / 2.0;
    double y_max = segHeight / 2.0;

    // 2. 히스토그램 생성
    std::string histName = "h_" + detName;
    std::string histTitle = detName + " Hit Distribution;X [mm];Y [mm];Hits";
    TH2F *h_hitmap = new TH2F(histName.c_str(), histTitle.c_str(),
                              (int)totalWidth, x_min, x_max,   // 1mm 간격으로 bin 설정
                              (int)segHeight, y_min, y_max);

    // 3. 브랜치 연결 및 이벤트 루프
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
    
    // 이 브랜치 주소는 다음 검출기를 위해 초기화해주는 것이 좋습니다.
    tree->ResetBranchAddresses();

    // 4. 캔버스에 그리기
    std::string canvasName = "c_" + detName;
    TCanvas *c = new TCanvas(canvasName.c_str(), (detName + " Hit Map").c_str(), 900, 500);
    c->cd();
    h_hitmap->SetStats(0);
    h_hitmap->Draw("COLZ");

    // 5. 경계선 그리기
    for (int i = 0; i <= nSegs; ++i) {
        double x_pos = x_min + (i * segWidth);
        TLine *line = new TLine(x_pos, y_min, x_pos, y_max);
        line->SetLineColor(kRed);
        line->Draw("SAME");
    }
    TLine* box_top = new TLine(x_min, y_max, x_max, y_max);
    TLine* box_bottom = new TLine(x_min, y_min, x_max, y_min);
    box_top->SetLineColor(kRed);
    box_bottom->SetLineColor(kRed);
    box_top->Draw("SAME");
    box_bottom->Draw("SAME");

    std::cout << detName << " 플롯 생성 완료." << std::endl;
}


// ## 메인 함수: 이 함수를 호출하면 됩니다 ##
void run_all_plots() {
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

    // 1. BVH_U 플롯 생성
    // 이름: "BVH_U", 세그먼트 수: 10, 세그먼트 폭: 10mm, 세그먼트 높이: 60mm
    plot_detector_hitmap(tree, "BVH_U", 10, 10.0, 60.0);

    // 2. BH2 플롯 생성
    // 이름: "BH2", 세그먼트 수: 12, 세그먼트 폭: 14mm, 세그먼트 높이: 100mm
    plot_detector_hitmap(tree, "BH2", 12, 14.0, 100.0);
}
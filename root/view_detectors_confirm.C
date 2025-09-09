#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TBox.h"
#include "TLine.h"
#include "TLatex.h"
#include "TMath.h"
#include "TParticle.h"
#include <vector>
#include <map>
#include <iostream>

//--- 오버레이 그리기 함수 (설정 변수를 인자로 받도록 수정) ---
void draw_bvh_overlay(double center_x, int n_segs, double seg_h, double seg_w = 10.0) {
    double total_w = n_segs * seg_w;
    double start_x = center_x - total_w / 2.0;
    
    for (int i = 0; i < n_segs; ++i) {
        double x1 = start_x + i * seg_w;
        TBox *seg = new TBox(x1, -seg_h/2.0, x1 + seg_w, seg_h/2.0);
        seg->SetFillStyle(0); seg->SetLineColor(kRed); seg->SetLineWidth(1);
        seg->Draw("SAME");
    }
    TLatex *label_zero = new TLatex(center_x, seg_h/2.0 + 15, "0");
    label_zero->SetTextAlign(22); label_zero->SetTextColor(kRed); label_zero->SetTextSize(0.03);
    label_zero->Draw("SAME");
}


// --- 메인 함수 ---
void view_detectors_confirm() {
    // ======[ 1. 설계 파라미터 설정 (이 부분만 수정하세요) ]======
    // BVH_U 설정
    double bvh_u_center_x = 40.0;    // mm, 오버레이 X축 중심
    int    bvh_u_n_segs   = 22;      // 세그먼트 개수
    double bvh_u_seg_h    = 140.0;   // mm, 세그먼트 높이 (60 또는 140)

    // BVH_D 설정
    double bvh_d_center_x = 350.0;   // mm, 오버레이 X축 중심
    int    bvh_d_n_segs   = 32;      // 세그먼트 개수
    double bvh_d_seg_h    = 140.0;   // mm, 세그먼트 높이 (60 또는 140)
    // ==========================================================

    TFile *f = new TFile("E45_newprofile15_980.root");
    if (!f || f->IsZombie()) { return; }
    TTree *tree = (TTree*)f->Get("g4hyptpc");
    if (!tree) { return; }

    std::map<int, double> vp_z_positions = {
        {1, -1000.0}, {2, -900.0}, {3, -800.0}, {4, -760.0},
        {5, 840.0},   {6, 900.0},  {7, 950.0},  {8, 1000.0}
    };
    
    // --- 오버레이 경계 자동 계산 ---
    double bvh_u_w = 10.0;
    double bvh_u_x_min = bvh_u_center_x - (bvh_u_n_segs * bvh_u_w) / 2.0;
    double bvh_u_x_max = bvh_u_center_x + (bvh_u_n_segs * bvh_u_w) / 2.0;
    double bvh_u_y_min = -bvh_u_seg_h / 2.0;
    double bvh_u_y_max =  bvh_u_seg_h / 2.0;

    double bvh_d_w = 10.0;
    double bvh_d_x_min = bvh_d_center_x - (bvh_d_n_segs * bvh_d_w) / 2.0;
    double bvh_d_x_max = bvh_d_center_x + (bvh_d_n_segs * bvh_d_w) / 2.0;
    double bvh_d_y_min = -bvh_d_seg_h / 2.0;
    double bvh_d_y_max =  bvh_d_seg_h / 2.0;

    // --- TTree 브랜치 및 카운터 설정 ---
    std::vector<TParticle> *vp_hits = nullptr;
    tree->SetBranchAddress("VP", &vp_hits);

    std::map<int, long long> total_hits;
    std::map<int, long long> outside_hits;
    for (int i = 1; i <= 8; ++i) { total_hits[i] = 0; outside_hits[i] = 0; }

    // --- 이벤트 루프 ---
    for (Long64_t i = 0; i < tree->GetEntries(); ++i) {
        tree->GetEntry(i);
        for (const auto& hit : *vp_hits) {
            double hitX = hit.Vx(); double hitY = hit.Vy(); double hitZ = hit.Vz();
            for (auto const& [vp_num, z_pos] : vp_z_positions) {
                if (TMath::Abs(hitZ - z_pos) < 1.0) {
                    total_hits[vp_num]++;
                    if (vp_num <= 4) { // BVH_U 영역
                        if (hitX < bvh_u_x_min || hitX > bvh_u_x_max || hitY < bvh_u_y_min || hitY > bvh_u_y_max) {
                            outside_hits[vp_num]++;
                        }
                    } else { // BVH_D 영역
                        if (hitX < bvh_d_x_min || hitX > bvh_d_x_max || hitY < bvh_d_y_min || hitY > bvh_d_y_max) {
                            outside_hits[vp_num]++;
                        }
                    }
                    break;
                }
            }
        }
    }

    // --- 결과 출력 ---
    std::cout << "=====================================================" << std::endl;
    std::cout << "      Detector Acceptance Analysis Results" << std::endl;
    std::cout << "=====================================================" << std::endl;
    std::cout << "--------- BVH_U Region (Upstream) ---------" << std::endl;
    for (int i = 1; i <= 4; ++i) {
        double percent_outside = (total_hits[i] > 0) ? (100.0 * outside_hits[i] / total_hits[i]) : 0;
        printf("[VP%d @ Z=%.1f mm]\n", i, vp_z_positions.at(i));
        printf(" - Total Hits        : %lld\n", total_hits[i]);
        printf(" - Hits Outside Overlay: %lld\n", outside_hits[i]);
        printf(" - Outside Percentage: %.2f %%\n\n", percent_outside);
    }
    std::cout << "--------- BVH_D Region (Downstream) ---------" << std::endl;
    for (int i = 5; i <= 8; ++i) {
        double percent_outside = (total_hits[i] > 0) ? (100.0 * outside_hits[i] / total_hits[i]) : 0;
        printf("[VP%d @ Z=%.1f mm]\n", i, vp_z_positions.at(i));
        printf(" - Total Hits        : %lld\n", total_hits[i]);
        printf(" - Hits Outside Overlay: %lld\n", outside_hits[i]);
        printf(" - Outside Percentage: %.2f %%\n\n", percent_outside);
    }
    std::cout << "=====================================================" << std::endl;

    // --- 히스토그램 생성 및 시각화 ---
    TH2D *h_vp1 = new TH2D("h_vp1", Form("BVH_U region(VP1, z=%.0f); X [mm]; Y [mm]", vp_z_positions[1]), 200, -300, 300, 200, -200, 200);
    TH2D *h_vp2 = new TH2D("h_vp2", Form("BVH_U region(VP2, z=%.0f); X [mm]; Y [mm]", vp_z_positions[2]), 200, -300, 300, 200, -200, 200);
    TH2D *h_vp3 = new TH2D("h_vp3", Form("BVH_U region(VP3, z=%.0f); X [mm]; Y [mm]", vp_z_positions[3]), 200, -300, 300, 200, -200, 200);
    TH2D *h_vp4 = new TH2D("h_vp4", Form("BVH_U region(VP4, z=%.0f); X [mm]; Y [mm]", vp_z_positions[4]), 200, -300, 300, 200, -200, 200);
    TH2D *h_vp5 = new TH2D("h_vp5", Form("BVH_D region(VP5, z=%.0f); X [mm]; Y [mm]", vp_z_positions[5]), 400, 0, 800, 200, -200, 200);
    TH2D *h_vp6 = new TH2D("h_vp6", Form("BVH_D region(VP6, z=%.0f); X [mm]; Y [mm]", vp_z_positions[6]), 400, 0, 800, 200, -200, 200);
    TH2D *h_vp7 = new TH2D("h_vp7", Form("BVH_D region(VP7, z=%.0f); X [mm]; Y [mm]", vp_z_positions[7]), 400, 0, 800, 200, -200, 200);
    TH2D *h_vp8 = new TH2D("h_vp8", Form("BVH_D region(VP8, z=%.0f); X [mm]; Y [mm]", vp_z_positions[8]), 400, 0, 800, 200, -200, 200);
    
    std::map<int, TH2D*> hist_map = { {1,h_vp1}, {2,h_vp2}, {3,h_vp3}, {4,h_vp4}, {5,h_vp5}, {6,h_vp6}, {7,h_vp7}, {8,h_vp8} };
    for(auto const& [vp_num, z_pos] : vp_z_positions) {
        TString selection = Form("abs(VP.Vz() - (%.1f)) < 1.0", z_pos);
        tree->Draw(Form("VP.Vy():VP.Vx() >> h_vp%d", vp_num), selection, "goff");
    }

    TCanvas *c_bvh_u = new TCanvas("c_bvh_u", "BVH U Region (VP1-4)", 1200, 1200);
    c_bvh_u->Divide(2, 2);
    for (int i = 1; i <= 4; ++i) {
        c_bvh_u->cd(i);
        hist_map[i]->Draw("COLZ");
        draw_bvh_overlay(bvh_u_center_x, bvh_u_n_segs, bvh_u_seg_h);
    }

    TCanvas *c_bvh_d = new TCanvas("c_bvh_d", "BVH D Region (VP5-8)", 1200, 1200);
    c_bvh_d->Divide(2, 2);
    for (int i = 5; i <= 8; ++i) {
        c_bvh_d->cd(i - 4);
        hist_map[i]->Draw("COLZ");
        draw_bvh_overlay(bvh_d_center_x, bvh_d_n_segs, bvh_d_seg_h);
    }
}
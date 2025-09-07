// view_detectors_final_with_all_labels_v2.C

#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TBox.h"
#include "TLine.h"
#include "TLatex.h"
#include "TMath.h"
#include <map>

//--- 오버레이 그리기 함수 (코드를 깔끔하게 하기 위해 분리) ---

// BVH_U (VP1, 2, 3, 4) 오버레이 그리기 함수
void draw_bvh_u_overlay() {
    double seg_w = 10.0, seg_h = 60.0; int n_segs = 20;
    double total_w = n_segs * seg_w, start_x = -total_w / 2.0;
    int center_idx = n_segs / 2;
    for (int i = 0; i < n_segs; ++i) {
        double x1 = start_x + i * seg_w;
        TBox *seg = new TBox(x1, -seg_h/2.0, x1 + seg_w, seg_h/2.0);
        seg->SetFillStyle(0); seg->SetLineColor(kRed); seg->SetLineWidth(1);
        seg->Draw("SAME");

        int seg_num = (i < center_idx) ? (i - center_idx) : (i - center_idx + 1);
        if (abs(seg_num) == 5 || abs(seg_num) == 10) {
            TLatex *label = new TLatex(x1 + seg_w/2.0, seg_h/2.0 + 10, Form("%+d", seg_num));
            label->SetTextAlign(22); label->SetTextColor(kRed); label->SetTextSize(0.03);
            label->Draw("SAME");
        }
    }
    TLatex *label_zero = new TLatex(0, seg_h/2.0 + 10, "0");
    label_zero->SetTextAlign(22); label_zero->SetTextColor(kRed); label_zero->SetTextSize(0.03);
    label_zero->Draw("SAME");
}

// BVH_D (VP5, 6, 7, 8) 오버레이 그리기 함수
void draw_bvh_d_overlay() {
    double seg_w = 10.0, seg_h = 140.0; int n_segs = 70;
    double total_w = n_segs * seg_w, start_x = -total_w / 2.0;
    int center_idx = n_segs / 2;
    for (int i = 0; i < n_segs; ++i) {
        double x1 = start_x + i * seg_w;
        TBox *seg = new TBox(x1, -seg_h/2.0, x1 + seg_w, seg_h/2.0);
        seg->SetFillStyle(0); seg->SetLineColor(kRed); seg->SetLineWidth(1);
        seg->Draw("SAME");
        int seg_num = (i < center_idx) ? (i - center_idx) : (i - center_idx + 1);
        if (seg_num != 0 && abs(seg_num) % 10 == 0) {
            TLatex *label = new TLatex(x1 + seg_w/2.0, seg_h/2.0 + 15, Form("%+d", seg_num));
            label->SetTextAlign(22); label->SetTextColor(kRed); label->SetTextSize(0.03);
            label->Draw("SAME");
        }
    }
    TLatex *label_zero = new TLatex(0, seg_h/2.0 + 15, "0");
    label_zero->SetTextAlign(22); label_zero->SetTextColor(kRed); label_zero->SetTextSize(0.03);
    label_zero->Draw("SAME");
}


// --- 메인 함수 ---
void view_detectors() {
    // 0) 파일/트리 열기
    TFile *f = new TFile("E45_newprofile4_980.root");
    if (!f || f->IsZombie()) { return; }
    TTree *tree = (TTree*)f->Get("g4hyptpc");
    if (!tree) { return; }

    // 1) 히스토그램 생성
    // Z 위치 정보 (DCGeomParam 파일 기준)
    std::map<int, double> vp_z_positions = {
        {1, -1000.0}, {2, -900.0}, {3, -800.0}, {4, -760.0},
        {5, 840.0},   {6, 900.0},  {7, 950.0},  {8, 1000.0}
    };

    // DetSize 파일 기준 VP 크기 (X, Y full length)
    // VP1-4: 800x800, VP5-8: 1000x1000
    // 히스토그램 범위는 detector 크기보다 약간 넓게 설정
    TH2D *h_vp1 = new TH2D("h_vp1", Form("BVH_U region(VP1, z=%.0f); X [mm]; Y [mm]", vp_z_positions[1]), 200, -410, 410, 200, -410, 410);
    TH2D *h_vp2 = new TH2D("h_vp2", Form("BVH_U region(VP2, z=%.0f); X [mm]; Y [mm]", vp_z_positions[2]), 200, -410, 410, 200, -410, 410);
    TH2D *h_vp3 = new TH2D("h_vp3", Form("BVH_U region(VP3, z=%.0f); X [mm]; Y [mm]", vp_z_positions[3]), 200, -410, 410, 200, -410, 410);
    TH2D *h_vp4 = new TH2D("h_vp4", Form("BVH_U region(VP4, z=%.0f); X [mm]; Y [mm]", vp_z_positions[4]), 200, -410, 410, 200, -410, 410);

    TH2D *h_vp5 = new TH2D("h_vp5", Form("BVH_D region(VP5, z=%.0f); X [mm]; Y [mm]", vp_z_positions[5]), 200, -510, 510, 200, -510, 510);
    TH2D *h_vp6 = new TH2D("h_vp6", Form("BVH_D region(VP6, z=%.0f); X [mm]; Y [mm]", vp_z_positions[6]), 200, -510, 510, 200, -510, 510);
    TH2D *h_vp7 = new TH2D("h_vp7", Form("BVH_D region(VP7, z=%.0f); X [mm]; Y [mm]", vp_z_positions[7]), 200, -510, 510, 200, -510, 510);
    TH2D *h_vp8 = new TH2D("h_vp8", Form("BVH_D region(VP8, z=%.0f); X [mm]; Y [mm]", vp_z_positions[8]), 200, -510, 510, 200, -510, 510);

    TH2D *h_t0  = new TH2D("h_t0",  "T0 Hit Distribution; X [mm]; Y [mm]", 200, -150, 150, 200, -150, 150);
    TH2D *h_bh2 = new TH2D("h_bh2", "BH2 Hit Distribution; X [mm]; Y [mm]", 200, -150, 150, 200, -150, 150);
    
    // 히스토그램 포인터 맵 (나중에 그리기 편하게)
    std::map<int, TH2D*> hist_map = {
        {1, h_vp1}, {2, h_vp2}, {3, h_vp3}, {4, h_vp4},
        {5, h_vp5}, {6, h_vp6}, {7, h_vp7}, {8, h_vp8}
    };

    // 2) 데이터 채우기
    for(auto const& [vp_num, z_pos] : vp_z_positions) {
        TString hist_name = Form("h_vp%d", vp_num);
        TString selection = Form("abs(VP.Vz() - (%.1f)) < 1.0", z_pos);
        tree->Draw(Form("VP.Vy():VP.Vx() >> %s", hist_name.Data()), selection, "goff");
    }
    tree->Draw("T0.Vy():T0.Vx() >> h_t0",  "", "goff");
    tree->Draw("BH2.Vy():BH2.Vx() >> h_bh2", "", "goff");

    // 3) 캔버스 생성 및 그리기 (총 10개 plot, 2x5)
    TCanvas *c1 = new TCanvas("c1", "Detector Hit Distributions", 1000, 2000);
    c1->Divide(2, 5);

    // VP 1-4 그리기 (BVH_U)
    for (int i = 1; i <= 4; ++i) {
        c1->cd(i);
        hist_map[i]->Draw("COLZ");
        draw_bvh_u_overlay();
    }

    // VP 5-8 그리기 (BVH_D)
    for (int i = 5; i <= 8; ++i) {
        c1->cd(i);
        hist_map[i]->Draw("COLZ");
        draw_bvh_d_overlay();
    }
    
    // BH2 그리기
    c1->cd(9);
    h_bh2->Draw("COLZ");
    {
        const double segW   = 14.0; const int nSeg = 15; const double cx = 12.0;
        const double halfW  = 0.5 * segW * nSeg; const double xL = cx - halfW;
        const double xR     = cx + halfW; const double halfY = 100.0/2.0;
        TBox *box_bh2 = new TBox(xL, -halfY, xR, halfY);
        box_bh2->SetFillStyle(0); box_bh2->SetLineColor(kRed); box_bh2->SetLineWidth(2);
        box_bh2->Draw("SAME");
        for (int i = 0; i <= nSeg; ++i) {
            double x = xL + i * segW;
            TLine *L = new TLine(x, -halfY, x, halfY);
            L->SetLineColor(kRed); L->SetLineStyle(3); L->Draw("SAME");
        }
        TLine *Lc = new TLine(cx, -halfY, cx, halfY);
        Lc->SetLineColor(kRed); Lc->SetLineStyle(2); Lc->Draw("SAME");
        TLatex *labL = new TLatex(xL, halfY+6, "-91");
        TLatex *labR = new TLatex(xR, halfY+6, "+63");
        labL->SetTextAlign(23); labR->SetTextAlign(23);
        labL->SetTextColor(kRed); labR->SetTextColor(kRed);
        labL->SetTextSize(0.03); labR->SetTextSize(0.03);
        labL->Draw("SAME"); labR->Draw("SAME");
    }

    // T0 그리기
    c1->cd(10);
    h_t0->Draw("COLZ");
    {
        double x_t0[4], y_t0[4], angle=45.0*TMath::Pi()/180.0, w=80, h=80;
        x_t0[0] = -w*cos(angle)-h*sin(angle); y_t0[0] = -w*sin(angle)+h*cos(angle);
        x_t0[1] =  w*cos(angle)-h*sin(angle); y_t0[1] =  w*sin(angle)+h*cos(angle);
        x_t0[2] =  w*cos(angle)+h*sin(angle); y_t0[2] =  w*sin(angle)-h*cos(angle);
        x_t0[3] = -w*cos(angle)+h*sin(angle); y_t0[3] = -w*sin(angle)-h*cos(angle);
        for(int i=0; i<4; ++i) {
            TLine *line = new TLine(x_t0[i], y_t0[i], x_t0[(i+1)%4], y_t0[(i+1)%4]);
            line->SetLineColor(kRed); line->SetLineWidth(2);
            line->Draw("SAME");
        }
    }
}
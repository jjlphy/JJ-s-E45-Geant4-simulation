// view_detectors_final_with_all_labels.C
#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TBox.h"
#include "TLine.h"
#include "TLatex.h"
#include "TMath.h"

void view_detectors() {
    TFile *f = new TFile("../test.root");
    if (!f || f->IsZombie()) { return; }
    TTree *tree = (TTree*)f->Get("g4hyptpc");
    if (!tree) { return; }

    // 1. 히스토그램 생성
    TH2D *h_vp4 = new TH2D("h_vp4", "BVH_U region(VP4); X [mm]; Y [mm]", 200, -150, 150, 200, -150, 150);
    TH2D *h_vp5 = new TH2D("h_vp5", "BVH_D region, close(VP5); X [mm]; Y [mm]", 200, -400, 400, 200, -400, 400);
    TH2D *h_vp6 = new TH2D("h_vp6", "BVH_D region, far(VP6); X [mm]; Y [mm]", 200, -400, 400, 200, -400, 400);
    TH2D *h_t0  = new TH2D("h_t0",  "T0 Hit Distribution; X [mm]; Y [mm]", 200, -150, 150, 200, -150, 150);
    TH2D *h_bh2 = new TH2D("h_bh2", "BH2 Hit Distribution; X [mm]; Y [mm]", 200, -150, 150, 200, -150, 150);

    // 2. 데이터 채우기
    tree->Draw("VP.Vy():VP.Vx() >> h_vp4", "abs(VP.Vz() - (-760.0)) < 1.0", "goff");
    tree->Draw("VP.Vy():VP.Vx() >> h_vp5", "abs(VP.Vz() - 470.0) < 1.0", "goff");
    tree->Draw("VP.Vy():VP.Vx() >> h_vp6", "abs(VP.Vz() - 490.0) < 1.0", "goff");
    tree->Draw("T0.Vy():T0.Vx() >> h_t0", "", "goff");
    tree->Draw("BH2.Vy():BH2.Vx() >> h_bh2", "", "goff");
    
    // 3. 캔버스 생성 및 그리기
    TCanvas *c1 = new TCanvas("c1", "Detector Hit Distributions", 1200, 1800);
    c1->Divide(2, 3);

    // --- BVH_U region(VP4) at z=-760 ---
    c1->cd(1);
    h_vp4->Draw("COLZ");
    double seg_w4 = 10.0, seg_h4 = 60.0; int n_segs4 = 20;
    double total_w4 = n_segs4 * seg_w4, start_x4 = -total_w4 / 2.0;
    int center_idx4 = n_segs4 / 2;
    for (int i = 0; i < n_segs4; ++i) {
        double x1 = start_x4 + i * seg_w4;
        TBox *seg = new TBox(x1, -seg_h4/2.0, x1 + seg_w4, seg_h4/2.0);
        seg->SetFillStyle(0); seg->SetLineColor(kRed); seg->SetLineWidth(1);
        seg->Draw("SAME");
        
        // [추가된 라벨링]
        int seg_num = (i < center_idx4) ? (i - center_idx4) : (i - center_idx4 + 1);
        if (abs(seg_num) == 5 || abs(seg_num) == 10) {
            TLatex *label = new TLatex(x1 + seg_w4/2.0, seg_h4/2.0 + 10, Form("%+d", seg_num));
            label->SetTextAlign(22); label->SetTextColor(kRed); label->SetTextSize(0.03);
            label->Draw("SAME");
        }
    }
    TLatex *label_zero4 = new TLatex(0, seg_h4/2.0 + 10, "0");
    label_zero4->SetTextAlign(22); label_zero4->SetTextColor(kRed); label_zero4->SetTextSize(0.03);
    label_zero4->Draw("SAME");


    // --- BH2 Hit Distribution ---
    c1->cd(2);
    h_bh2->Draw("COLZ");
    TBox *box_bh2 = new TBox(-84, -50, 84, 50);
    box_bh2->SetFillStyle(0); box_bh2->SetLineColor(kRed); box_bh2->SetLineWidth(2);
    box_bh2->Draw("SAME");

    // --- BVH_D region, close(VP5) at z=470 ---
    c1->cd(3);
    h_vp5->Draw("COLZ");
    double seg_w5 = 10.0, seg_h5 = 140.0; int n_segs5 = 70;
    double total_w5 = n_segs5 * seg_w5, start_x5 = -total_w5 / 2.0;
    int center_idx5 = n_segs5 / 2;
    for (int i = 0; i < n_segs5; ++i) {
        double x1 = start_x5 + i * seg_w5;
        TBox *seg = new TBox(x1, -seg_h5/2.0, x1 + seg_w5, seg_h5/2.0);
        seg->SetFillStyle(0); seg->SetLineColor(kRed); seg->SetLineWidth(1);
        seg->Draw("SAME");
        int seg_num = (i < center_idx5) ? (i - center_idx5) : (i - center_idx5 + 1);
        if (seg_num != 0 && abs(seg_num) % 10 == 0) {
            TLatex *label = new TLatex(x1 + seg_w5/2.0, seg_h5/2.0 + 15, Form("%+d", seg_num));
            label->SetTextAlign(22); label->SetTextColor(kRed); label->SetTextSize(0.03);
            label->Draw("SAME");
        }
    }
    TLatex *label_zero5 = new TLatex(0, seg_h5/2.0 + 15, "0");
    label_zero5->SetTextAlign(22); label_zero5->SetTextColor(kRed); label_zero5->SetTextSize(0.03);
    label_zero5->Draw("SAME");


    // --- BVH_D region, far(VP6) at z=490 ---
    c1->cd(4);
    h_vp6->Draw("COLZ");
    double seg_w6 = 10.0, seg_h6 = 140.0; int n_segs6 = 70;
    double total_w6 = n_segs6 * seg_w6, start_x6 = -total_w6 / 2.0;
    int center_idx6 = n_segs6 / 2;
    for (int i = 0; i < n_segs6; ++i) {
        double x1 = start_x6 + i * seg_w6;
        TBox *seg = new TBox(x1, -seg_h6/2.0, x1 + seg_w6, seg_h6/2.0);
        seg->SetFillStyle(0); seg->SetLineColor(kRed); seg->SetLineWidth(1);
        seg->Draw("SAME");
        int seg_num = (i < center_idx6) ? (i - center_idx6) : (i - center_idx6 + 1);
        if (seg_num != 0 && abs(seg_num) % 10 == 0) {
            TLatex *label = new TLatex(x1 + seg_w6/2.0, seg_h6/2.0 + 15, Form("%+d", seg_num));
            label->SetTextAlign(22); label->SetTextColor(kRed); label->SetTextSize(0.03);
            label->Draw("SAME");
        }
    }
    TLatex *label_zero6 = new TLatex(0, seg_h6/2.0 + 15, "0");
    label_zero6->SetTextAlign(22); label_zero6->SetTextColor(kRed); label_zero6->SetTextSize(0.03);
    label_zero6->Draw("SAME");
    
    // --- T0 Hit Distribution ---
    c1->cd(5);
    h_t0->Draw("COLZ");
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
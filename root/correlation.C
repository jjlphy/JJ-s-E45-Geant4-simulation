#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TParticle.h"
#include "TStyle.h"
#include <vector>
#include <iostream>

void correlation() {
    // 스타일 설정
    gStyle->SetOptStat(1111);
    gStyle->SetPalette(57); // kBird

    // 1. 파일 및 TTree 열기
    TFile *f = new TFile("E45_BVH1.root"); // 분석할 파일 이름
    if (!f || f->IsZombie()) { return; }
    TTree *tree = (TTree*)f->Get("g4hyptpc");
    if (!tree) { return; }

    // 2. 브랜치 읽기 설정
    std::vector<TParticle> *bh2_hits = nullptr;
    std::vector<TParticle> *bvh_u_hits = nullptr;
    std::vector<TParticle> *bvh_d_hits = nullptr;
    std::vector<TParticle> *htof_hits = nullptr;

    tree->SetBranchAddress("BH2", &bh2_hits);
    tree->SetBranchAddress("BVH_U", &bvh_u_hits);
    tree->SetBranchAddress("BVH_D", &bvh_d_hits);
    tree->SetBranchAddress("HTOF", &htof_hits);

    // 3. 히스토그램 정의
    // 1D edep 스펙트럼
    TH1D *h_edep_bh2 = new TH1D("h_edep_bh2", "BH2 Energy Deposit;Edep [MeV];Entries", 200, 0, 4.0);
    TH1D *h_edep_bvh_u = new TH1D("h_edep_bvh_u", "BVH_U Energy Deposit;Edep [MeV];Entries", 200, 0, 2.0);
    TH1D *h_edep_bvh_d = new TH1D("h_edep_bvh_d", "BVH_D Energy Deposit;Edep [MeV];Entries", 200, 0, 2.0);
    TH1D *h_edep_htof = new TH1D("h_edep_htof", "HTOF Energy Deposit;Edep [MeV];Entries", 200, 0, 10.0);
    
    // 2D Hit 분포
    TH2D *h_hit_bh2   = new TH2D("h_hit_bh2",   "BH2 Hit Distribution;X [mm];Y [mm]", 200, -100, 150, 200, -100, 100);
    TH2D *h_hit_bvh_u = new TH2D("h_hit_bvh_u", "BVH_U Hit Distribution;X [mm];Y [mm]", 200, -150, 150, 200, -150, 150);
    TH2D *h_hit_bvh_d = new TH2D("h_hit_bvh_d", "BVH_D Hit Distribution;X [mm];Y [mm]", 400, 150, 650, 200, -150, 150);

    // 2D 상관관계 히스토그램
    TH2D *h_corr_bh2_bvh_u = new TH2D("h_corr_bh2_bvh_u", "BH2 vs BVH_U (HTOF multiplicity #geq 2);BH2 Seg ID;BVH_U Seg ID", 15, -0.5, 14.5, 26, -0.5, 25.5);
    TH2D *h_corr_bh2_bvh_d = new TH2D("h_corr_bh2_bvh_d", "BH2 vs BVH_D (HTOF multiplicity #geq 2);BH2 Seg ID;BVH_D Seg ID", 15, -0.5, 14.5, 40, -0.5, 39.5);
    TH2D *h_corr_bvh_u_bvh_d = new TH2D("h_corr_bvh_u_bvh_d", "BVH_U vs BVH_D (HTOF multiplicity #geq 2);BVH_U Seg ID;BVH_D Seg ID", 26, -0.5, 25.5, 40, -0.5, 39.5);
    
    // 4. 이벤트 루프 및 히스토그램 채우기
    long long n_entries = tree->GetEntries();
    for (long long i = 0; i < n_entries; ++i) {
        tree->GetEntry(i);

        // 0.1 MIP 기준 edep 제한 값 (MeV)
        const double edep_cut_bh2 = 0.1; // 1.0 MeV * 0.1
        const double edep_cut_bvh = 0.04; // 0.4 MeV * 0.1
        const double edep_cut_htof = 0.2; // 2.0 MeV * 0.1

        int n_bh2_hits = 0, n_bvh_u_hits = 0, n_bvh_d_hits = 0, n_htof_hits = 0;
        int bh2_id = -1, bvh_u_id = -1, bvh_d_id = -1;
        
        for (const auto& hit : *bh2_hits) {
            h_edep_bh2->Fill(hit.GetWeight());
            if (hit.GetWeight() > edep_cut_bh2) {
                n_bh2_hits++;
                bh2_id = hit.GetMother(0);
                h_hit_bh2->Fill(hit.Vx(), hit.Vy());
            }
        }
        for (const auto& hit : *bvh_u_hits) {
            h_edep_bvh_u->Fill(hit.GetWeight());
            if (hit.GetWeight() > edep_cut_bvh) {
                n_bvh_u_hits++;
                bvh_u_id = hit.GetMother(0);
                h_hit_bvh_u->Fill(hit.Vx(), hit.Vy());
            }
        }
        for (const auto& hit : *bvh_d_hits) {
            h_edep_bvh_d->Fill(hit.GetWeight());
            if (hit.GetWeight() > edep_cut_bvh) {
                n_bvh_d_hits++;
                bvh_d_id = hit.GetMother(0);
                h_hit_bvh_d->Fill(hit.Vx(), hit.Vy());
            }
        }
        for (const auto& hit : *htof_hits) {
            h_edep_htof->Fill(hit.GetWeight());
            if (hit.GetWeight() > edep_cut_htof) {
                n_htof_hits++;
            }
        }
        
        if (n_htof_hits >= 2) {
            if (n_bh2_hits == 1 && n_bvh_u_hits == 1 && n_bvh_d_hits == 1) {
                h_corr_bh2_bvh_u->Fill(bh2_id, bvh_u_id);
                h_corr_bh2_bvh_d->Fill(bh2_id, bvh_d_id);
                h_corr_bvh_u_bvh_d->Fill(bvh_u_id, bvh_d_id);
            }
        }
    }

    // 5. 시각화
    TCanvas *c_edep = new TCanvas("c_edep", "Energy Deposition Spectra", 1200, 800);
    c_edep->Divide(2, 2);
    c_edep->cd(1); gPad->SetLogy(); h_edep_bh2->Draw();
    c_edep->cd(2); gPad->SetLogy(); h_edep_bvh_u->Draw();
    c_edep->cd(3); gPad->SetLogy(); h_edep_bvh_d->Draw();
    c_edep->cd(4); gPad->SetLogy(); h_edep_htof->Draw();
    
    TCanvas *c_hits = new TCanvas("c_hits", "2D Hit Distributions", 1800, 600);
    c_hits->Divide(3, 1);
    c_hits->cd(1); h_hit_bh2->Draw("COLZ");
    c_hits->cd(2); h_hit_bvh_u->Draw("COLZ");
    c_hits->cd(3); h_hit_bvh_d->Draw("COLZ");

    TCanvas *c_corr = new TCanvas("c_corr", "Correlations (HTOF multiplicity >= 2)", 1500, 500);
    c_corr->Divide(3, 1);
    c_corr->cd(1); h_corr_bh2_bvh_u->Draw("COLZ");
    c_corr->cd(2); h_corr_bh2_bvh_d->Draw("COLZ");
    c_corr->cd(3); h_corr_bvh_u_bvh_d->Draw("COLZ");
}
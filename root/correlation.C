// correlation_labels.C
// 1) No HTOF cut 3개 (캔버스 c_no)
// 2) HTOF multiplicity >= multCut 3개 (캔버스 c_cut)
// 각 2D bin에 "전체 이벤트 대비 %"(작은 글씨) 오버레이.
// 사용법:
// root -l
// root [0] .L gen.dict.C+                         // 필요시 딕셔너리
// root [1] .L correlation_labels.C+
// root [2] correlation_labels("E45_BVH2.root");   // 실행

#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TParticle.h"
#include "TStyle.h"
#include "TLine.h"
#include "TLatex.h"
#include <vector>
#include <unordered_set>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <string>

#define SEG(P)   (P.GetMother(0))   // 세그 인덱스(로그로 확인됨)
#define EDEP(P)  (P.GetWeight())    // HTOF 에너지(Weight)

static void DrawSegGrid(TH2* h, int color=kGray+2, int style=3, int width=1){
    if(!h) return;
    auto xa = h->GetXaxis();
    auto ya = h->GetYaxis();
    double x1 = xa->GetXmin(), x2 = xa->GetXmax();
    double y1 = ya->GetXmin(), y2 = ya->GetXmax();
    int xlo = (int)std::floor(x1+0.5), xhi = (int)std::ceil(x2-0.5);
    int ylo = (int)std::floor(y1+0.5), yhi = (int)std::ceil(y2-0.5);
    for(int xi=xlo; xi<=xhi; ++xi){
        double xv = xi + 0.5;
        TLine *lv = new TLine(xv, y1, xv, y2);
        lv->SetLineColor(color); lv->SetLineStyle(style); lv->SetLineWidth(width); lv->Draw("same");
    }
    for(int yi=ylo; yi<=yhi; ++yi){
        double yv = yi + 0.5;
        TLine *lh = new TLine(x1, yv, x2, yv);
        lh->SetLineColor(color); lh->SetLineStyle(style); lh->SetLineWidth(width); lh->Draw("same");
    }
}

// bin 중앙에 작은 글씨로 % 텍스트를 찍는다(0.1% 미만 생략 가능)
static void DrawPercentLabels(TH2* h, const std::vector<std::vector<int>>& occCount,
                              Long64_t totalEvents, double minPercentToDraw=0.1,
                              double textSize=0.018, int color=kBlack)
{
    if(!h || totalEvents<=0) return;
    TLatex tx; tx.SetNDC(false); tx.SetTextSize(textSize); tx.SetTextColor(color);

    int nx = h->GetNbinsX();
    int ny = h->GetNbinsY();
    for(int ix=1; ix<=nx; ++ix){
        for(int iy=1; iy<=ny; ++iy){
            int cnt = occCount[ix-1][iy-1];
            if(cnt<=0) continue;
            double pct = 100.0 * (double)cnt / (double)totalEvents;
            if(pct < minPercentToDraw) continue; // 너무 작은 건 생략해서 가독성 유지
            double x = h->GetXaxis()->GetBinCenter(ix);
            double y = h->GetYaxis()->GetBinCenter(iy);
            tx.DrawLatex(x, y, Form("%.1f%%", pct));
        }
    }
}

void correlation(const char* filename="E45_BVH2.root",
                        int multCut=2,            // HTOF multiplicity 기준
                        bool useEnergyCut=true,   // true: Edep>edepCut 후 개수 셈
                        double edepCut=0.2,       // HTOF 에너지 컷
                        double minPercentToDraw=0.1, // 이 % 미만 라벨 생략(가독성)
                        bool drawGrid=true,       // 세그 경계 점선
                        bool verbose=true)
{
    gStyle->SetOptStat(1111);
    gStyle->SetPalette(kBird);

    TFile* f = TFile::Open(filename, "READ");
    if(!f || f->IsZombie()){ std::cerr<<"[ERR] open "<<filename<<" failed\n"; return; }
    TTree* tree = dynamic_cast<TTree*>(f->Get("g4hyptpc"));
    if(!tree){ std::cerr<<"[ERR] TTree g4hyptpc not found\n"; return; }

    std::vector<TParticle> *bh2=nullptr, *bvhu=nullptr, *bvhd=nullptr, *htof=nullptr;
    tree->SetBranchAddress("BH2",   &bh2);
    tree->SetBranchAddress("BVH_U", &bvhu);
    tree->SetBranchAddress("BVH_D", &bvhd);
    tree->SetBranchAddress("HTOF",  &htof);

    const Long64_t N = tree->GetEntries();
    if(N==0){ std::cerr<<"[WARN] 0 entries\n"; return; }

    // 프리스캔으로 실제 세그 범위 확인 → binning 고정
    int loBH2=+1000000000, hiBH2=-1000000000;
    int loU  =+1000000000, hiU  =-1000000000;
    int loD  =+1000000000, hiD  =-1000000000;
    Long64_t Nscan = std::min<Long64_t>(N, 50000);
    for(Long64_t i=0;i<Nscan;++i){
        tree->GetEntry(i);
        if(bh2)  for(const auto& p:*bh2){ int id=SEG(p); loBH2=std::min(loBH2,id); hiBH2=std::max(hiBH2,id); }
        if(bvhu) for(const auto& p:*bvhu){int id=SEG(p); loU  =std::min(loU,  id); hiU  =std::max(hiU,  id); }
        if(bvhd) for(const auto& p:*bvhd){int id=SEG(p); loD  =std::min(loD,  id); hiD  =std::max(hiD,  id); }
    }
    auto fix=[&](int& lo,int& hi){ if(lo>hi){ lo=0; hi=0; } };
    fix(loBH2,hiBH2); fix(loU,hiU); fix(loD,hiD);

    int nxBH2 = std::max(1, hiBH2-loBH2+1);
    int nyU   = std::max(1, hiU  -loU  +1);
    int nyD   = std::max(1, hiD  -loD  +1);

    double xloBH2 = loBH2-0.5, xhiBH2 = hiBH2+0.5;
    double yloU   = loU  -0.5, yhiU   = hiU  +0.5;
    double yloD   = loD  -0.5, yhiD   = hiD  +0.5;

    // 히스토그램(heatmap) 6개
    TH2D* h_bh2_u_no = new TH2D("h_bh2_u_no", "BH2 vs BVH_U (No HTOF cut);BH2 Seg ID;BVH_U Seg ID",
                                nxBH2,xloBH2,xhiBH2, nyU,yloU,yhiU);
    TH2D* h_bh2_d_no = new TH2D("h_bh2_d_no", "BH2 vs BVH_D (No HTOF cut);BH2 Seg ID;BVH_D Seg ID",
                                nxBH2,xloBH2,xhiBH2, nyD,yloD,yhiD);
    TH2D* h_u_d_no   = new TH2D("h_u_d_no",   "BVH_U vs BVH_D (No HTOF cut);BVH_U Seg ID;BVH_D Seg ID",
                                nyU,yloU,yhiU, nyD,yloD,yhiD);

    TH2D* h_bh2_u_ct = new TH2D("h_bh2_u_ct", Form("BH2 vs BVH_U (HTOF #geq %d);BH2 Seg ID;BVH_U Seg ID",multCut),
                                nxBH2,xloBH2,xhiBH2, nyU,yloU,yhiU);
    TH2D* h_bh2_d_ct = new TH2D("h_bh2_d_ct", Form("BH2 vs BVH_D (HTOF #geq %d);BH2 Seg ID;BVH_D Seg ID",multCut),
                                nxBH2,xloBH2,xhiBH2, nyD,yloD,yhiD);
    TH2D* h_u_d_ct   = new TH2D("h_u_d_ct",   Form("BVH_U vs BVH_D (HTOF #geq %d);BVH_U Seg ID;BVH_D Seg ID",multCut),
                                nyU,yloU,yhiU, nyD,yloD,yhiD);

    // ===== 이벤트 단위 bin 점유 카운트(퍼센트용) =====
    // same shape as histograms: [ix][iy]에 "이 bin을 점유한 이벤트 수" 저장
    std::vector<std::vector<int>> occ_no_bu(nxBH2, std::vector<int>(nyU, 0));  // BH2 vs U (no cut)
    std::vector<std::vector<int>> occ_no_bd(nxBH2, std::vector<int>(nyD, 0));  // BH2 vs D (no cut)
    std::vector<std::vector<int>> occ_no_ud(nyU,   std::vector<int>(nyD, 0));  // U   vs D (no cut)

    std::vector<std::vector<int>> occ_ct_bu(nxBH2, std::vector<int>(nyU, 0));  // BH2 vs U (cut)
    std::vector<std::vector<int>> occ_ct_bd(nxBH2, std::vector<int>(nyD, 0));  // BH2 vs D (cut)
    std::vector<std::vector<int>> occ_ct_ud(nyU,   std::vector<int>(nyD, 0));  // U   vs D (cut)

    auto within = [](int ix,int nx,int iy,int ny)->bool{
        return (ix>=0 && ix<nx && iy>=0 && iy<ny);
    };

    Long64_t n_pass = 0;

    // ===== 메인 루프 =====
    for(Long64_t i=0;i<N;++i){
        tree->GetEntry(i);

        // --- No cut: 이 이벤트가 점유한 bin 집합(중복 제거) ---
        std::unordered_set<long long> seen_bu_no, seen_bd_no, seen_ud_no;
        if(bh2 && bvhu){
            for(const auto& a:*bh2){
                int ia = SEG(a)-loBH2;
                for(const auto& b:*bvhu){
                    int ib = SEG(b)-loU;
                    if(!within(ia,nxBH2, ib,nyU)) continue;
                    long long key = ((long long)ia<<32) | (unsigned long long)ib;
                    seen_bu_no.insert(key);
                    h_bh2_u_no->Fill(ia+loBH2, ib+loU); // heatmap(히트수 시각화)
                }
            }
        }
        if(bh2 && bvhd){
            for(const auto& a:*bh2){
                int ia = SEG(a)-loBH2;
                for(const auto& b:*bvhd){
                    int ib = SEG(b)-loD;
                    if(!within(ia,nxBH2, ib,nyD)) continue;
                    long long key = ((long long)ia<<32) | (unsigned long long)ib;
                    seen_bd_no.insert(key);
                    h_bh2_d_no->Fill(ia+loBH2, ib+loD);
                }
            }
        }
        if(bvhu && bvhd){
            for(const auto& a:*bvhu){
                int ia = SEG(a)-loU;
                for(const auto& b:*bvhd){
                    int ib = SEG(b)-loD;
                    if(!within(ia,nyU, ib,nyD)) continue;
                    long long key = ((long long)ia<<32) | (unsigned long long)ib;
                    seen_ud_no.insert(key);
                    h_u_d_no->Fill(ia+loU, ib+loD);
                }
            }
        }
        // 이벤트 단위 카운트(+1씩)
        for(auto key: seen_bu_no){ int ia= (int)(key>>32), ib=(int)(key&0xffffffff); occ_no_bu[ia][ib]++; }
        for(auto key: seen_bd_no){ int ia= (int)(key>>32), ib=(int)(key&0xffffffff); occ_no_bd[ia][ib]++; }
        for(auto key: seen_ud_no){ int ia= (int)(key>>32), ib=(int)(key&0xffffffff); occ_no_ud[ia][ib]++; }

        // --- HTOF ≥ multCut 조건 판단 ---
        int mult = 0;
        if(htof){
            if(useEnergyCut){ for(const auto& h : *htof) if(EDEP(h)>edepCut) ++mult; }
            else            { mult = (int)htof->size(); }
        }
        if(mult < multCut) continue;
        ++n_pass;

        // --- Cut 적용 버전(이벤트 단위 중복 제거) ---
        std::unordered_set<long long> seen_bu_ct, seen_bd_ct, seen_ud_ct;
        if(bh2 && bvhu){
            for(const auto& a:*bh2){
                int ia = SEG(a)-loBH2;
                for(const auto& b:*bvhu){
                    int ib = SEG(b)-loU;
                    if(!within(ia,nxBH2, ib,nyU)) continue;
                    long long key = ((long long)ia<<32) | (unsigned long long)ib;
                    seen_bu_ct.insert(key);
                    h_bh2_u_ct->Fill(ia+loBH2, ib+loU);
                }
            }
        }
        if(bh2 && bvhd){
            for(const auto& a:*bh2){
                int ia = SEG(a)-loBH2;
                for(const auto& b:*bvhd){
                    int ib = SEG(b)-loD;
                    if(!within(ia,nxBH2, ib,nyD)) continue;
                    long long key = ((long long)ia<<32) | (unsigned long long)ib;
                    seen_bd_ct.insert(key);
                    h_bh2_d_ct->Fill(ia+loBH2, ib+loD);
                }
            }
        }
        if(bvhu && bvhd){
            for(const auto& a:*bvhu){
                int ia = SEG(a)-loU;
                for(const auto& b:*bvhd){
                    int ib = SEG(b)-loD;
                    if(!within(ia,nyU, ib,nyD)) continue;
                    long long key = ((long long)ia<<32) | (unsigned long long)ib;
                    seen_ud_ct.insert(key);
                    h_u_d_ct->Fill(ia+loU, ib+loD);
                }
            }
        }
        for(auto key: seen_bu_ct){ int ia= (int)(key>>32), ib=(int)(key&0xffffffff); occ_ct_bu[ia][ib]++; }
        for(auto key: seen_bd_ct){ int ia= (int)(key>>32), ib=(int)(key&0xffffffff); occ_ct_bd[ia][ib]++; }
        for(auto key: seen_ud_ct){ int ia= (int)(key>>32), ib=(int)(key&0xffffffff); occ_ct_ud[ia][ib]++; }
    }

    if(verbose){
        std::cout << "\n=====================================================\n";
        std::cout << " File                      : " << filename << "\n";
        std::cout << " Total Events              : " << N << "\n";
        std::cout << " Events with HTOF >= " << multCut
                  << (useEnergyCut? Form(" (Edep>%.3f)", edepCut):" (count only)")
                  << " : " << n_pass << "\n";
        std::cout << " Labels show: (events_bin / TotalEvents="<<N<<") * 100\n";
        std::cout << "=====================================================\n";
    }

    // ---- 캔버스 1: No HTOF cut (3개) + bin별 % 라벨 ----
    TCanvas* c_no = new TCanvas("c_no", "No HTOF cut (event % labels)", 1800, 600);
    c_no->Divide(3,1);
    c_no->cd(1); h_bh2_u_no->Draw("COLZ"); if(drawGrid) DrawSegGrid(h_bh2_u_no);
               DrawPercentLabels(h_bh2_u_no, occ_no_bu, N, minPercentToDraw, 0.018, kBlack);
    c_no->cd(2); h_bh2_d_no->Draw("COLZ"); if(drawGrid) DrawSegGrid(h_bh2_d_no);
               DrawPercentLabels(h_bh2_d_no, occ_no_bd, N, minPercentToDraw, 0.018, kBlack);
    c_no->cd(3); h_u_d_no  ->Draw("COLZ"); if(drawGrid) DrawSegGrid(h_u_d_no);
               DrawPercentLabels(h_u_d_no,   occ_no_ud, N, minPercentToDraw, 0.018, kBlack);
    c_no->Update();

    // ---- 캔버스 2: HTOF >= multCut (3개) + bin별 % 라벨 ----
    TCanvas* c_ct = new TCanvas("c_ct", Form("HTOF #geq %d (event %% labels)", multCut), 1800, 600);
    c_ct->Divide(3,1);
    c_ct->cd(1); h_bh2_u_ct->Draw("COLZ"); if(drawGrid) DrawSegGrid(h_bh2_u_ct);
               DrawPercentLabels(h_bh2_u_ct, occ_ct_bu, N, minPercentToDraw, 0.018, kRed+1);
    c_ct->cd(2); h_bh2_d_ct->Draw("COLZ"); if(drawGrid) DrawSegGrid(h_bh2_d_ct);
               DrawPercentLabels(h_bh2_d_ct, occ_ct_bd, N, minPercentToDraw, 0.018, kRed+1);
    c_ct->cd(3); h_u_d_ct  ->Draw("COLZ"); if(drawGrid) DrawSegGrid(h_u_d_ct);
               DrawPercentLabels(h_u_d_ct,   occ_ct_ud, N, minPercentToDraw, 0.018, kRed+1);
    c_ct->Update();

    // 저장 원하면:
    // c_no->SaveAs("noHTOF_with_labels.pdf");
    // c_ct->SaveAs("withHTOF_with_labels.pdf");
}

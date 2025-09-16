// correlation_BH2.C  (compact axes & labels + per-detector Edep thresholds)
// root -l
// .L correlation_BH2.C+
// correlation_BH2("E45_BVH2.root",
//   /*dedupPerEvent=*/true,
//   /*htofEdepCut=*/0.0, /*htofMultCut=*/0,
//   /*labelMinPct=*/0.2,
//   /*labelPctOfBH2=*/true,
//   /*logz=*/false, /*save=*/false,
//   /*axisLabelSize=*/0.035, /*axisTitleSize=*/0.050, /*overlayTextSize=*/0.020,
//   /*thrBH2=*/0.00, /*thrBVHU=*/0.00, /*thrBVHD=*/0.00); // ← 새로 추가된 3개(기본 0)

#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TParticle.h"
#include "TLatex.h"
#include "TString.h"
#include "TPaletteAxis.h"
#include <vector>
#include <unordered_set>
#include <algorithm>
#include <iostream>

static const int N_BH2  = 15; // 0..14
static const int N_BVHU = 22; // 0..21
static const int N_BVHD = 32; // 0..31

// ---------- overlay ----------
static void drawCellLabels(TH2* h, double denom,
                           double minPct=0.2, double textSize=0.020)
{
  if(!h || denom<=0) return;
  TLatex tx; tx.SetTextAlign(22); tx.SetTextFont(42); tx.SetTextSize(textSize);
  const int nx=h->GetNbinsX(), ny=h->GetNbinsY();
  for(int ix=1; ix<=nx; ++ix){
    for(int iy=1; iy<=ny; ++iy){
      const double n=h->GetBinContent(ix,iy);
      if(n<=0) continue;
      const double pct=100.0*n/denom;
      if(pct<minPct) continue;
      tx.DrawLatex(h->GetXaxis()->GetBinCenter(ix),
                   h->GetYaxis()->GetBinCenter(iy),
                   Form("%.1f%% (%.0f)", pct, n));
    }
  }
}

// ---------- axis labels ----------
static void setSegmentLabels(TH2D* h, int nU, int nD,
                             double axisLabelSize, double axisTitleSize)
{
  for(int i=1;i<=nU;++i) h->GetXaxis()->SetBinLabel(i, Form("%d", i-1));
  for(int j=1;j<=nD;++j) h->GetYaxis()->SetBinLabel(j, Form("%d", j-1));
  h->GetXaxis()->LabelsOption("v");      // vertical (겹침 방지)

  h->GetXaxis()->SetLabelSize(axisLabelSize);
  h->GetYaxis()->SetLabelSize(axisLabelSize);
  h->GetXaxis()->SetTitleSize(axisTitleSize);
  h->GetYaxis()->SetTitleSize(axisTitleSize);
  h->GetXaxis()->SetTitleOffset(1.0);
  h->GetYaxis()->SetTitleOffset(1.0);
  h->GetZaxis()->SetLabelSize(0.030);    // 팔레트 라벨 기본 크기
}

// ---------- histo factory ----------
static TH2D* makeUD(const char* name, const char* title,
                    int nU, int nD,
                    double axisLabelSize, double axisTitleSize)
{
  TH2D* h = new TH2D(name, title,
                     nU, -0.5, nU-0.5,
                     nD, -0.5, nD-0.5);
  h->GetXaxis()->SetTitle("BVH_U Seg");
  h->GetYaxis()->SetTitle("BVH_D Seg");
  setSegmentLabels(h, nU, nD, axisLabelSize, axisTitleSize);
  return h;
}

// ---------- unique segments per event (no-threshold; 유지용) ----------
static std::vector<int> uniqSegs(const std::vector<TParticle>* v, int nmax)
{
  std::unordered_set<int> s;
  if(v){
    for(const auto& p:*v){
      int id = p.GetMother(1);
      if(0<=id && id<nmax) s.insert(id);
    }
  }
  return std::vector<int>(s.begin(), s.end());
}

// ---------- unique segments with Edep threshold (새로 추가) ----------
static std::vector<int> uniqSegsThr(const std::vector<TParticle>* v, int nmax, double thr)
{
  std::unordered_set<int> s;
  if(v){
    for(const auto& p:*v){
      if(p.GetWeight() <= thr) continue;         // 임계값 적용
      int id = p.GetMother(1);
      if(0<=id && id<nmax) s.insert(id);
    }
  }
  return std::vector<int>(s.begin(), s.end());
}

void correlation_BH2(const char* fname="E45_BVH%.root",
                        bool   dedupPerEvent   = true,
                        double htofEdepCut    = 0.0,
                        int    htofMultCut    = 0,
                        double labelMinPct    = 0.2,
                        bool   labelPctOfBH2  = true,
                        bool   logz           = false,
                        bool   save           = false,
                        double axisLabelSize  = 0.035,
                        double axisTitleSize  = 0.050,
                        double overlayTextSize= 0.020,
                        // ===== 여기부터 새 인자(기본 0.0 = 필터 없음) =====
                        double thrBH2         = 0.0,
                        double thrBVHU        = 0.0,
                        double thrBVHD        = 0.0)
{
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kBird);

  TFile* f = TFile::Open(fname,"READ");
  if(!f || f->IsZombie()){ std::cerr<<"[ERR] cannot open "<<fname<<"\n"; return; }
  TTree* tr = (TTree*)f->Get("g4hyptpc");
  if(!tr){ std::cerr<<"[ERR] no TTree g4hyptpc\n"; return; }

  std::vector<TParticle> *bh2=nullptr, *bvhu=nullptr, *bvhd=nullptr, *htof=nullptr;
  tr->SetBranchAddress("BH2",   &bh2);
  tr->SetBranchAddress("BVH_U", &bvhu);
  tr->SetBranchAddress("BVH_D", &bvhd);
  tr->SetBranchAddress("HTOF",  &htof);

  // histos per BH2 seg
  std::vector<TH2D*> hud(N_BH2, nullptr);
  for(int h=0; h<N_BH2; ++h){
    hud[h] = makeUD(Form("h_ud_bh2_%d",h),
                    Form("BH2=%d",h),
                    N_BVHU, N_BVHD, axisLabelSize, axisTitleSize);
  }

  const Long64_t N = tr->GetEntries();
  Long64_t N_pass = 0;
  std::vector<Long64_t> N_bh2_seg(N_BH2, 0); // denominator

  for(Long64_t i=0;i<N;++i){
    tr->GetEntry(i);

    bool use_evt = true;
    if(htofMultCut>0){
      int mHTOF=0; if(htof) for(const auto& h:*htof) if(h.GetWeight()>htofEdepCut) ++mHTOF;
      use_evt = (mHTOF>=htofMultCut);
      if(use_evt) ++N_pass;
    }
    if(!use_evt) continue;

    // --- 세그먼트 수집 (임계값 적용) ---
    std::vector<int> U, D, H;
    if(dedupPerEvent){
      U = uniqSegsThr(bvhu, N_BVHU, thrBVHU);
      D = uniqSegsThr(bvhd, N_BVHD, thrBVHD);
      H = uniqSegsThr(bh2,  N_BH2,  thrBH2 );
    }else{
      if(bvhu) for(const auto& p:*bvhu){ if(p.GetWeight()<=thrBVHU) continue; int id=p.GetMother(1); if(0<=id&&id<N_BVHU) U.push_back(id); }
      if(bvhd) for(const auto& p:*bvhd){ if(p.GetWeight()<=thrBVHD) continue; int id=p.GetMother(1); if(0<=id&&id<N_BVHD) D.push_back(id); }
      if(bh2)  for(const auto& p:*bh2 ){ if(p.GetWeight()<=thrBH2 ) continue; int id=p.GetMother(1); if(0<=id&&id<N_BH2 ) H.push_back(id); }
    }

    // 분모: 해당 BH2 세그에 에너지-컷을 통과한 "이벤트" 개수
    // (중복 제거된 H 기준)
    for(int h: H) if(0<=h && h<N_BH2) ++N_bh2_seg[h];

    // 채우기
    for(int h: H){
      if(h<0||h>=N_BH2) continue;
      for(int u: U)
        for(int d: D)
          hud[h]->Fill(u,d);
    }
  }

  // 5 canvases × 3 panels
  const double denomAll = (double)((htofMultCut>0)?N_pass:N);
  for(int g=0; g<5; ++g){
    TCanvas* c = new TCanvas(Form("c_grp_%d",g),
                             Form("BVH_U vs BVH_D per BH2 (group %d)",g),
                             2100, 700);
    c->Divide(3,1, 0.002, 0.002);

    for(int j=0;j<3;++j){
      int h = g*3 + j;
      c->cd(j+1);
      if(logz) gPad->SetLogz();
      gPad->SetGrid(1,1);
      gPad->SetLeftMargin(0.13);
      gPad->SetRightMargin(0.12);
      gPad->SetBottomMargin(0.18);

      hud[h]->Draw("COLZ");

      // shrink palette labels a bit
      if(auto pal = (TPaletteAxis*)hud[h]->GetListOfFunctions()->FindObject("palette")){
        pal->SetLabelSize(0.030);
        pal->SetTitleSize(0.035);
      }

      double denom = labelPctOfBH2 ? (double)std::max<Long64_t>(1, N_bh2_seg[h])
                                   : std::max(1.0, denomAll);
      drawCellLabels(hud[h], denom, labelMinPct, overlayTextSize);
    }

    if(save){
      c->SaveAs(Form("bh2_slices_grp%d.png",g));
      c->SaveAs(Form("bh2_slices_grp%d.pdf",g));
    }
  }

  // summary
  std::cout << "\n==================== Summary ====================\n";
  std::cout << "File                : " << fname << "\n";
  std::cout << "Total events        : " << N << "\n";
  if(htofMultCut>0)
    std::cout << "HTOF mult >= " << htofMultCut << " : " << N_pass
              << " (" << (100.0*N_pass/std::max<Long64_t>(1,N)) << " %)\n";
  std::cout << "Per-BH2 events (denominator):\n";
  for(int h=0; h<N_BH2; ++h) std::cout << "  BH2 " << h << " : " << N_bh2_seg[h] << "\n";
  std::cout << "Label basis         : "
            << (labelPctOfBH2 ? "per-BH2 events" : ((htofMultCut>0)?"events passing cut":"all events"))
            << "\n";
  std::cout << "Edep thr (BH2/U/D)  : " << thrBH2 << " / " << thrBVHU << " / " << thrBVHD << "\n";
  std::cout << "=================================================\n";
}

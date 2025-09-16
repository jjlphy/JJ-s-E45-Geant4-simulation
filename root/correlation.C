// correlation.C
// root -l
// .L correlation.C+
// correlation("E45_BVH2.root", true, 0.2, 2, 0.2, true, true, false, false);

#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TParticle.h"
#include "TLatex.h"
#include "TString.h"
#include <vector>
#include <unordered_set>
#include <algorithm>
#include <iostream>

static const int N_BH2  = 15; // 0..14
static const int N_BVHU = 22; // 0..21 (UPSTREAM)
static const int N_BVHD = 32; // 0..31 (DOWNSTREAM)

// -------- cell labels -------------------------------------------------------
static void drawCellLabels(TH2* h, double denom,
                           double minPct=0.2, double textSize=0.018)
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

// -------- 2D histo factory (모든 bin 라벨 + 'none' 추가) --------------------
static TH2D* mk2(const char* name,const char* title,
                 int nx,int ny, bool withNone,
                 const char* xlab, const char* ylab)
{
  const int    bx  = withNone ? nx+1 : nx;
  const int    by  = withNone ? ny+1 : ny;
  const double xlo = -0.5, xhi = withNone ? nx+0.5 : nx-0.5;
  const double ylo = -0.5, yhi = withNone ? ny+0.5 : ny-0.5;

  TH2D* h = new TH2D(name, Form("%s;%s;%s", title, xlab, ylab),
                     bx, xlo, xhi, by, ylo, yhi);

  // 👉 bin 라벨을 0..N-1로 모두 채움
  for(int i=1;i<=nx;++i) h->GetXaxis()->SetBinLabel(i, Form("%d", i-1));
  for(int j=1;j<=ny;++j) h->GetYaxis()->SetBinLabel(j, Form("%d", j-1));
  if(withNone){
    h->GetXaxis()->SetBinLabel(nx+1,"none");
    h->GetYaxis()->SetBinLabel(ny+1,"none");
  }

  // 보기 옵션
  h->GetXaxis()->LabelsOption("v"); // x축 라벨 세로 회전(겹침 방지)
  h->GetXaxis()->SetLabelSize(0.018);
  h->GetYaxis()->SetLabelSize(0.018);
  h->GetXaxis()->SetNdivisions(510);
  h->GetYaxis()->SetNdivisions(510);
  return h;
}

void correlation(const char* fname="E45_BVH%.root",
                 bool   dedupPerEvent   = true,
                 double htofEdepCut    = 0.2,
                 int    htofMultCut    = 2,
                 double labelMinPct    = 0.2,
                 bool   labelPctOfEvents = true,
                 bool   includeNoHit   = true,
                 bool   logz           = false,
                 bool   save           = false)
{
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kBird);

  // open
  TFile* f = TFile::Open(fname,"READ");
  if(!f || f->IsZombie()){ std::cerr<<"[ERR] cannot open "<<fname<<"\n"; return; }
  TTree* tr = (TTree*)f->Get("g4hyptpc");
  if(!tr){ std::cerr<<"[ERR] no TTree g4hyptpc\n"; return; }

  std::vector<TParticle> *bh2=nullptr, *bvhu=nullptr, *bvhd=nullptr, *htof=nullptr;
  tr->SetBranchAddress("BH2",   &bh2);
  tr->SetBranchAddress("BVH_U", &bvhu);
  tr->SetBranchAddress("BVH_D", &bvhd);
  tr->SetBranchAddress("HTOF",  &htof);

  // histos
  TH2D *h_bh2_u_nocut = mk2("h_bh2_u_nocut","BH2 vs BVH_U (no HTOF cut)",
                             N_BH2, N_BVHU, includeNoHit, "BH2 Seg", "BVH_U Seg");
  TH2D *h_bh2_d_nocut = mk2("h_bh2_d_nocut","BH2 vs BVH_D (no HTOF cut)",
                             N_BH2, N_BVHD, includeNoHit, "BH2 Seg", "BVH_D Seg");
  TH2D *h_u_d_nocut   = mk2("h_u_d_nocut",  "BVH_U vs BVH_D (no HTOF cut)",
                             N_BVHU, N_BVHD, includeNoHit, "BVH_U Seg", "BVH_D Seg");

  TH2D *h_bh2_u_ct = mk2("h_bh2_u_ct",Form("BH2 vs BVH_U (HTOF multiplicity #geq %d)",htofMultCut),
                          N_BH2, N_BVHU, includeNoHit, "BH2 Seg", "BVH_U Seg");
  TH2D *h_bh2_d_ct = mk2("h_bh2_d_ct",Form("BH2 vs BVH_D (HTOF multiplicity #geq %d)",htofMultCut),
                          N_BH2, N_BVHD, includeNoHit, "BH2 Seg", "BVH_D Seg");
  TH2D *h_u_d_ct   = mk2("h_u_d_ct",  Form("BVH_U vs BVH_D (HTOF multiplicity #geq %d)",htofMultCut),
                          N_BVHU, N_BVHD, includeNoHit, "BVH_U Seg", "BVH_D Seg");

  const int BH2_NONE  = N_BH2;
  const int BVHU_NONE = N_BVHU;
  const int BVHD_NONE = N_BVHD;

  const Long64_t N = tr->GetEntries();
  Long64_t N_pass = 0;

  for(Long64_t i=0;i<N;++i){
    tr->GetEntry(i);

    // HTOF mult
    int mHTOF=0;
    if(htof) for(const auto& h:*htof) if(h.GetWeight()>htofEdepCut) ++mHTOF;
    const bool passCT = (mHTOF>=htofMultCut);
    if(passCT) ++N_pass;

    // collect segs
    std::vector<int> vBH2, vU, vD;
    if(dedupPerEvent){
      std::unordered_set<int> sBH2,sU,sD;
      if(bh2)  for(const auto& p:*bh2){ int id=p.GetMother(1); if(0<=id && id<N_BH2)  sBH2.insert(id); }
      if(bvhu) for(const auto& p:*bvhu){int id=p.GetMother(1); if(0<=id && id<N_BVHU) sU.insert(id);  }
      if(bvhd) for(const auto& p:*bvhd){int id=p.GetMother(1); if(0<=id && id<N_BVHD) sD.insert(id);  }
      vBH2.assign(sBH2.begin(), sBH2.end());
      vU.assign(sU.begin(), sU.end());
      vD.assign(sD.begin(), sD.end());
    }else{
      if(bh2)  for(const auto& p:*bh2){ int id=p.GetMother(1); if(0<=id && id<N_BH2)  vBH2.push_back(id); }
      if(bvhu) for(const auto& p:*bvhu){int id=p.GetMother(1); if(0<=id && id<N_BVHU) vU.push_back(id);  }
      if(bvhd) for(const auto& p:*bvhd){int id=p.GetMother(1); if(0<=id && id<N_BVHD) vD.push_back(id);  }
    }

    // no-hit bucket
    if(includeNoHit){
      if(vBH2.empty()) vBH2.push_back(BH2_NONE);
      if(vU.empty())   vU.push_back(BVHU_NONE);
      if(vD.empty())   vD.push_back(BVHD_NONE);
    }

    // fill
    for(int a: vBH2){
      for(int b: vU){
        h_bh2_u_nocut->Fill(a,b);
        if(passCT) h_bh2_u_ct->Fill(a,b);
      }
      for(int c: vD){
        h_bh2_d_nocut->Fill(a,c);
        if(passCT) h_bh2_d_ct->Fill(a,c);
      }
    }
    for(int b: vU){
      for(int c: vD){
        h_u_d_nocut->Fill(b,c);
        if(passCT) h_u_d_ct->Fill(b,c);
      }
    }
  }

  auto draw_triple = [&](TH2D* h1, TH2D* h2, TH2D* h3,
                         double denomEvents,
                         const char* cname, const char* ctitle){
    TCanvas* c = new TCanvas(cname, ctitle, 1800, 600);
    c->Divide(3,1);
    if(logz){ c->cd(1)->SetLogz(); c->cd(2)->SetLogz(); c->cd(3)->SetLogz(); }

    c->cd(1); gPad->SetGrid(1,1); h1->Draw("COLZ");
    drawCellLabels(h1, labelPctOfEvents ? denomEvents : h1->Integral(), labelMinPct);
    c->cd(2); gPad->SetGrid(1,1); h2->Draw("COLZ");
    drawCellLabels(h2, labelPctOfEvents ? denomEvents : h2->Integral(), labelMinPct);
    c->cd(3); gPad->SetGrid(1,1); h3->Draw("COLZ");
    drawCellLabels(h3, labelPctOfEvents ? denomEvents : h3->Integral(), labelMinPct);
    return c;
  };

  TCanvas* c_nc = draw_triple(h_bh2_u_nocut, h_bh2_d_nocut, h_u_d_nocut,
                              (double)N, "c_corr_nocut", "BH2/BVH correlations (no HTOF cut)");
  TCanvas* c_ct = draw_triple(h_bh2_u_ct,    h_bh2_d_ct,    h_u_d_ct,
                              (double)N_pass, "c_corr_cut",
                              Form("BH2/BVH correlations (HTOF mult #geq %d)", htofMultCut));

  std::cout << "\n=====================================================\n";
  std::cout << " File                        : " << fname << "\n";
  std::cout << " Total Events                : " << N << "\n";
  std::cout << " HTOF mult >= " << htofMultCut << " events : " << N_pass
            << "  (" << (100.0 * N_pass / std::max<Long64_t>(1, N)) << " %)\n";
  std::cout << " Label basis                 : "
            << (labelPctOfEvents ? "percent of events" : "percent of hist counts") << "\n";
  std::cout << " Include 'no hit' bin        : " << (includeNoHit?"true":"false") << "\n";
  std::cout << "=====================================================\n";

  if(save){
    c_nc->SaveAs("corr_nocut.png");  c_nc->SaveAs("corr_nocut.pdf");
    c_ct->SaveAs("corr_htofcut.png"); c_ct->SaveAs("corr_htofcut.pdf");
  }
}

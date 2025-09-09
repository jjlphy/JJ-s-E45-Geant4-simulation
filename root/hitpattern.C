// hitpattern.C
// Usage:
//   root -l
//   .L hitpattern.C+
//   hitpattern("E45_BVH2.root", /*distinctPerEvent=*/false, /*save=*/true);

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TParticle.h"
#include <vector>
#include <unordered_set>
#include <iostream>
#include <algorithm>

// Segment counts (fixed binning 0..N-1)
static const int N_BH2  = 15; // 0..14
static const int N_BVHU = 26; // 0..25
static const int N_BVHD = 40; // 0..39

// Pretty helper
static void set_style() {
  gStyle->SetOptStat(1111);
  gStyle->SetTitleFont(42, ""); gStyle->SetLabelFont(42, "XY");
  gStyle->SetTitleSize(0.05, "XY"); gStyle->SetLabelSize(0.04, "XY");
}

static void printUF(TH1* h) {
  const int nb = h->GetNbinsX();
  std::cout << "  underflow=" << h->GetBinContent(0)
            << ", overflow=" << h->GetBinContent(nb+1) << "\n";
}

// main
void hitpattern(const char* fname="E45_BVH3.root",
                bool distinctPerEvent=false,   // true: 이벤트당 같은 세그는 1회만 카운트
                bool save=false,               // true: png/pdf 저장
                Long64_t maxEvents=-1)         // -1: 전 이벤트
{
  set_style();

  // Open tree
  TFile* f = TFile::Open(fname,"READ");
  if(!f || f->IsZombie()){ std::cerr<<"[ERR] cannot open "<<fname<<"\n"; return; }
  TTree* tr = (TTree*)f->Get("g4hyptpc");
  if(!tr){ std::cerr<<"[ERR] no TTree g4hyptpc\n"; return; }

  // Branches
  std::vector<TParticle> *bh2=nullptr, *bvhu=nullptr, *bvhd=nullptr;
  tr->SetBranchAddress("BH2",   &bh2);
  tr->SetBranchAddress("BVH_U", &bvhu);
  tr->SetBranchAddress("BVH_D", &bvhd);

  // Histos (exactly 0..N-1)
  TH1D* h_bh2  = new TH1D("hBH2",  "BH2 HitPattern;BH2 SegID;counts",      N_BH2,  -0.5, N_BH2-0.5);
  TH1D* h_bvhu = new TH1D("hBVHU", "BVH_U HitPattern;BVH_U SegID;counts",  N_BVHU, -0.5, N_BVHU-0.5);
  TH1D* h_bvhd = new TH1D("hBVHD", "BVH_D HitPattern;BVH_D SegID;counts",  N_BVHD, -0.5, N_BVHD-0.5);

  h_bh2 ->SetFillColorAlpha(kAzure+1, 0.65);
  h_bvhu->SetFillColorAlpha(kSpring+5,0.65);
  h_bvhd->SetFillColorAlpha(kOrange+7,0.65);

  // For quick sanity
  int minU=+1e9, maxU=-1e9, minD=+1e9, maxD=-1e9;

  const Long64_t Ntot = tr->GetEntries();
  const Long64_t Nuse = (maxEvents<0 ? Ntot : std::min(maxEvents, Ntot));

  for(Long64_t i=0;i<Nuse;++i){
    tr->GetEntry(i);

    if(!distinctPerEvent){
      // Per-hit counting (duplicates allowed)
      if(bh2)  for(const auto& p:*bh2){
        int id = p.GetMother(1);                // <<<< segment ID lives in Mother(1)
        if(0<=id && id<N_BH2) h_bh2->Fill(id);
      }
      if(bvhu) for(const auto& p:*bvhu){
        int id = p.GetMother(1);
        minU = std::min(minU,id); maxU = std::max(maxU,id);
        if(0<=id && id<N_BVHU) h_bvhu->Fill(id);
      }
      if(bvhd) for(const auto& p:*bvhd){
        int id = p.GetMother(1);
        minD = std::min(minD,id); maxD = std::max(maxD,id);
        if(0<=id && id<N_BVHD) h_bvhd->Fill(id);
      }
    } else {
      // Per-event counting (deduplicate within an event)
      std::unordered_set<int> sBH2, sU, sD;

      if(bh2)  for(const auto& p:*bh2){
        int id = p.GetMother(1);
        if(0<=id && id<N_BH2) sBH2.insert(id);     // dedup here
      }
      if(bvhu) for(const auto& p:*bvhu){
        int id = p.GetMother(1);
        minU = std::min(minU,id); maxU = std::max(maxU,id);
        if(0<=id && id<N_BVHU) sU.insert(id);      // dedup here
      }
      if(bvhd) for(const auto& p:*bvhd){
        int id = p.GetMother(1);
        minD = std::min(minD,id); maxD = std::max(maxD,id);
        if(0<=id && id<N_BVHD) sD.insert(id);      // dedup here
      }

      for(int v: sBH2) h_bh2 ->Fill(v);
      for(int v: sU)   h_bvhu->Fill(v);
      for(int v: sD)   h_bvhd->Fill(v);
    }
  }

  // Console summary
  std::cout << "\n===== HitPattern summary ("<<fname<<") =====\n";
  std::cout << "Events used: " << Nuse << " / " << Ntot << "\n";
  std::cout << "BH2   integral=" << h_bh2 ->Integral() << " ";  printUF(h_bh2);
  std::cout << "BVH_U integral=" << h_bvhu->Integral() << " ";  printUF(h_bvhu);
  std::cout << "BVH_D integral=" << h_bvhd->Integral() << " ";  printUF(h_bvhd);
  if(minU<=maxU) std::cout << "Observed BVH_U seg range (Mother(1)): " << minU << ".." << maxU << "\n";
  if(minD<=maxD) std::cout << "Observed BVH_D seg range (Mother(1)): " << minD << ".." << maxD << "\n";
  std::cout << "distinctPerEvent=" << (distinctPerEvent?"true":"false") << "\n";

  // Draw
  TCanvas* c = new TCanvas("c_hitpat","HitPatterns",1500,450);
  c->Divide(3,1);
  c->cd(1); h_bh2 ->Draw("hist");
  c->cd(2); h_bvhu->Draw("hist");
  c->cd(3); h_bvhd->Draw("hist");
  c->Update();

  if(save){
    c->SaveAs("hitpatterns.png");
    c->SaveAs("hitpatterns.pdf");
  }
}

// hitpattern.C
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TParticle.h"
#include <vector>
#include <iostream>

static const int N_BH2  = 15; // 0..14
static const int N_BVHU = 22; // 0..21  <-- 수정
static const int N_BVHD = 32; // 0..31  <-- 수정

void hitpattern(const char* fname="../E45_with_SCH.root",
                bool DISTINCT_PER_EVENT=false,
                bool SAVE=false)
{
  gStyle->SetOptStat(1111);

  TFile* f = TFile::Open(fname,"READ");
  if(!f || f->IsZombie()){ std::cerr<<"[ERR] open "<<fname<<" failed\n"; return; }
  TTree* tr = (TTree*)f->Get("g4hyptpc");
  if(!tr){ std::cerr<<"[ERR] no TTree g4hyptpc\n"; return; }

  std::vector<TParticle> *bh2=nullptr, *bvhu=nullptr, *bvhd=nullptr;
  tr->SetBranchAddress("BH2",   &bh2);
  tr->SetBranchAddress("BVH_U", &bvhu);
  tr->SetBranchAddress("BVH_D", &bvhd);

  TH1D* h_bh2  = new TH1D("hBH2",  "BH2 HitPattern;BH2 SegID;counts",     N_BH2,  -0.5, N_BH2-0.5);
  TH1D* h_bvhu = new TH1D("hBVHU", "BVH_U HitPattern;BVH_U SegID;counts", N_BVHU, -0.5, N_BVHU-0.5);
  TH1D* h_bvhd = new TH1D("hBVHD", "BVH_D HitPattern;BVH_D SegID;counts", N_BVHD, -0.5, N_BVHD-0.5);

  const Long64_t N = tr->GetEntries();
  for(Long64_t i=0;i<N;++i){
    tr->GetEntry(i);

    if(!DISTINCT_PER_EVENT){
      if(bh2)  for(const auto& p:*bh2){  int id=p.GetMother(1); if(0<=id && id<N_BH2)  h_bh2->Fill(id); }
      if(bvhu) for(const auto& p:*bvhu){ int id=p.GetMother(1); if(0<=id && id<N_BVHU) h_bvhu->Fill(id); }
      if(bvhd) for(const auto& p:*bvhd){ int id=p.GetMother(1); if(0<=id && id<N_BVHD) h_bvhd->Fill(id); }
    }else{
      bool seenBH2[ N_BH2 ]  = {0};
      bool seenU [ N_BVHU ]  = {0};
      bool seenD [ N_BVHD ]  = {0};
      if(bh2)  for(const auto& p:*bh2){  int id=p.GetMother(1); if(0<=id && id<N_BH2)  seenBH2[id]=true; }
      if(bvhu) for(const auto& p:*bvhu){ int id=p.GetMother(1); if(0<=id && id<N_BVHU) seenU[id]=true;  }
      if(bvhd) for(const auto& p:*bvhd){ int id=p.GetMother(1); if(0<=id && id<N_BVHD) seenD[id]=true;  }
      for(int i=0;i<N_BH2; ++i)  if(seenBH2[i]) h_bh2->Fill(i);
      for(int i=0;i<N_BVHU;++i)  if(seenU[i])   h_bvhu->Fill(i);
      for(int i=0;i<N_BVHD;++i)  if(seenD[i])   h_bvhd->Fill(i);
    }
  }

  TCanvas* c = new TCanvas("c_hitpat","HitPatterns",1500,450);
  c->Divide(3,1);
  c->cd(1); h_bh2 ->SetFillColorAlpha(kAzure+1,0.6);  h_bh2 ->Draw("hist");
  c->cd(2); h_bvhu->SetFillColorAlpha(kSpring+5,0.6); h_bvhu->Draw("hist");
  c->cd(3); h_bvhd->SetFillColorAlpha(kOrange+7,0.6); h_bvhd->Draw("hist");
  if(SAVE){ c->SaveAs("hitpatterns.png"); c->SaveAs("hitpatterns.pdf"); }
}


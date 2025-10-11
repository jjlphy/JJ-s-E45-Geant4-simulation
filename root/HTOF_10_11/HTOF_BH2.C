// HTOF_BH2_Gated_Summary.C
// 1) Total events
// 2) BH2 pass
// 3) BH2 pass & HTOF >=1
// 4) BH2 pass & HTOF >=2
// 5) BH2-gated HTOF HitPattern (unique segments/event), bins: 0..33 + NONE
//
// Usage:
//   root -l
//   .L HTOF_BH2_Gated_Summary.C+
//   HTOF_BH2_Gated_Summary("E45.root","g4hyptpc", /*save=*/true);

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TSystem.h"
#include "TInterpreter.h"
#include "TParticle.h"
#include "TString.h"
#include <vector>
#include <set>
#include <iostream>
#include <iomanip>

static const int kNTiles   = 34;   // 0..33
static const int kBIN_NONE = 34;   // 34th bin: NONE

struct _DICT_ {
  _DICT_(){
    gSystem->Load("libPhysics");
    gInterpreter->GenerateDictionary("vector<TParticle>","TParticle.h;vector");
  }
} _dict_;

static bool HasBranch(TTree* tr, const char* bname){
  return tr && tr->GetBranch(bname);
}

static inline double pct(double num, double den){
  return (den>0) ? (100.0 * num / den) : 0.0;
}

void HTOF_BH2_Gated_Summary(const char* filename="E45.root",
                            const char* treename="g4hyptpc",
                            bool save=true)
{
  // open
  TFile* f = TFile::Open(filename,"READ");
  if(!f || f->IsZombie()){ std::cerr<<"[ERR] open "<<filename<<" failed\n"; return; }
  TTree* T = (TTree*)f->Get(treename);
  if(!T){ std::cerr<<"[ERR] tree "<<treename<<" not found\n"; return; }

  // branches
  if(!HasBranch(T,"BH2")){ std::cerr<<"[ERR] need branch 'BH2' (vector<TParticle>)\n"; return; }
  const bool hasHTOFvec  = HasBranch(T,"HTOF");
  const bool hasHTOFcopy = HasBranch(T,"HTOF_copyNo"); // optional

  if(!(hasHTOFvec || hasHTOFcopy)){
    std::cerr<<"[ERR] need 'HTOF' (vector<TParticle>) or 'HTOF_copyNo' (vector<int>)\n";
    return;
  }

  std::vector<TParticle>* BH2 = nullptr;
  std::vector<TParticle>* HTOF = nullptr;
  std::vector<int>* HTOF_copyNo = nullptr;

  T->SetBranchAddress("BH2",&BH2);
  if(hasHTOFvec)  T->SetBranchAddress("HTOF",&HTOF);
  if(hasHTOFcopy) T->SetBranchAddress("HTOF_copyNo",&HTOF_copyNo);

  // histogram: 0..33 + NONE
  TH1D* hHitPat = new TH1D("hHTOF_hitpattern_BH2",
     "BH2-gated HTOF HitPattern (unique segments/event);Tile ID;Events",
     kNTiles+1, -0.5, kNTiles+0.5);
  hHitPat->SetDirectory(nullptr);
  hHitPat->Sumw2();
  // label all bins
  auto ax = hHitPat->GetXaxis();
  for (int b = 1; b <= kNTiles; ++b) ax->SetBinLabel(b, Form("%d", b-1)); // 1->0, ..., 34->33
  ax->SetBinLabel(kNTiles+1, "NONE");
  ax->LabelsOption("v");    // vertical labels

  // counters
  const Long64_t Ntot = T->GetEntries();
  Long64_t N_bh2 = 0;
  Long64_t N_bh2_htof_ge1 = 0;
  Long64_t N_bh2_htof_ge2 = 0;

  // loop
  for(Long64_t ie=0; ie<Ntot; ++ie){
    T->GetEntry(ie);
    if(!BH2) continue;

    // BH2 gate
    if(BH2->empty()) continue;
    N_bh2++;

    // unique HTOF tiles for this event
    std::set<int> uniq;
    if(hasHTOFcopy && HTOF_copyNo){
      for(int cn : *HTOF_copyNo){
        if(0<=cn && cn<kNTiles) uniq.insert(cn);
      }
    }else if(hasHTOFvec && HTOF){
      for(const auto& p : *HTOF){
        int cn = p.GetStatusCode();
        if(0<=cn && cn<kNTiles) uniq.insert(cn);
      }
    }

    if(uniq.empty()){
      hHitPat->Fill(kBIN_NONE);      // NONE
    }else{
      N_bh2_htof_ge1++;
      if((int)uniq.size()>=2) N_bh2_htof_ge2++;
      for(int cn : uniq) hHitPat->Fill(cn);
    }
  }

  // report (각 항목에 대해: 전체 기준 %, BH2 기준 % 같이 표기)
  std::cout<<std::fixed<<std::setprecision(2);
  std::cout<<"========== BH2-gated HTOF Summary ==========\n";
  std::cout<<"Total events                      : "<<Ntot<<"\n";
  std::cout<<"BH2 pass                          : "<<N_bh2
           <<" ("<<pct(N_bh2,Ntot)<<" % of Total)\n";
  std::cout<<"BH2 pass & HTOF >= 1              : "<<N_bh2_htof_ge1
           <<" ("<<pct(N_bh2_htof_ge1,Ntot)<<" % of Total, "
           <<pct(N_bh2_htof_ge1,N_bh2)<<" % of BH2)\n";
  std::cout<<"BH2 pass & HTOF >= 2              : "<<N_bh2_htof_ge2
           <<" ("<<pct(N_bh2_htof_ge2,Ntot)<<" % of Total, "
           <<pct(N_bh2_htof_ge2,N_bh2)<<" % of BH2)\n";
  std::cout<<"============================================\n";

  // draw
  TCanvas* c = new TCanvas("cBH2_HTOF","BH2-gated HTOF Summary", 950, 520);
  c->SetBottomMargin(0.20);   // for vertical labels
  hHitPat->SetLineWidth(2);
  hHitPat->Draw("hist");

  if(save){
    c->SaveAs("BH2_gated_HTOF_HitPattern.png");
  }
}

// HTOF_HitPat_copynum.C
// - Assumes HTOFSD/HTOFHit stores per-hit copy number into TParticle::StatusCode()
// - Draws 34-bin hit pattern (0..33) and multiplicity per event
//
// Usage:
//   root -l
//   .L HTOF_HitPat_copynum.C+
//   HTOF_HitPat_copynum("E45.root","g4hyptpc", /*uniquePerEvent=*/true, /*save=*/true);

#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TSystem.h"
#include "TInterpreter.h"
#include "TParticle.h"
#include <vector>
#include <set>
#include <iostream>
#include <iomanip>

// 34 tiles: copy_no = 0..33
static const int kNTiles = 34;

// Optional: also show multiplicity distribution up to this cap
static const int kMultCap = 12;

// Prepare dictionary for vector<TParticle> if needed
struct _DICT_ {
  _DICT_(){
    gSystem->Load("libPhysics");
    gInterpreter->GenerateDictionary("vector<TParticle>","TParticle.h;vector");
  }
} _dict_;

bool HasBranch(TTree* tr, const char* bname){
  return tr && tr->GetBranch(bname);
}

void HTOF_HitPat_copynum(const char* filename="E45.root",
                         const char* treename="g4hyptpc",
                         bool uniquePerEvent=true,
                         bool save=true)
{
  // open file & tree
  TFile* f = TFile::Open(filename, "READ");
  if(!f || f->IsZombie()){ std::cerr << "[ERR] Cannot open: " << filename << "\n"; return; }
  TTree* T = (TTree*)f->Get(treename);
  if(!T){ std::cerr << "[ERR] Tree not found: " << treename << "\n"; return; }

  // Prefer vector<int> HTOF_copyNo if exists; else use vector<TParticle> HTOF (StatusCode)
  const bool has_vecint = HasBranch(T, "HTOF_copyNo");
  const bool has_vecpar = HasBranch(T, "HTOF");

  if(!(has_vecint || has_vecpar)){
    std::cerr << "[ERR] Need branch 'HTOF_copyNo' (vector<int>) or 'HTOF' (vector<TParticle>)\n";
    T->Print();
    return;
  }

  std::vector<int>* v_copy = nullptr;
  std::vector<TParticle>* v_part = nullptr;
  if(has_vecint) T->SetBranchAddress("HTOF_copyNo", &v_copy);
  if(has_vecpar) T->SetBranchAddress("HTOF", &v_part);

  // Histograms
  TH1D* hTile = new TH1D("hHTOF_tile34",
                 "HTOF Hit Pattern (copy_no);Tile ID (0..33);Counts",
                 kNTiles, -0.5, kNTiles-0.5);

  TH1D* hMult = new TH1D("hHTOF_multiplicity",
                 "HTOF Multiplicity per Event;# unique tiles per event;Events",
                 kMultCap+1, -0.5, kMultCap+0.5);

  // Loop
  const Long64_t N = T->GetEntries();
  Long64_t n_evt = 0, n_evt_ge2 = 0;
  for(Long64_t ie=0; ie<N; ++ie){
    T->GetEntry(ie);

    std::set<int> uniq; // unique tiles for this event

    if(has_vecint && v_copy){
      for(int cn : *v_copy){
        if(cn>=0 && cn<kNTiles){
          if(uniquePerEvent){ uniq.insert(cn); } else { hTile->Fill(cn); }
        }
      }
    }
    else if(has_vecpar && v_part){
      for(const auto& p : *v_part){
        int cn = p.GetStatusCode(); // filled in HTOFSD patch
        if(cn>=0 && cn<kNTiles){
          if(uniquePerEvent){ uniq.insert(cn); } else { hTile->Fill(cn); }
        }
      }
    }

    if(uniquePerEvent){
      for(int cn : uniq) hTile->Fill(cn);
      int m = (int)uniq.size();
      if(m > kMultCap) m = kMultCap;
      hMult->Fill(m);
      if(!uniq.empty()){ n_evt++; if((int)uniq.size()>=2) n_evt_ge2++; }
    }
  }

  // Report
  if(uniquePerEvent){
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "[INFO] Events with >=1 HTOF hit : " << n_evt << "\n";
    std::cout << "[INFO] Multiplicity >= 2        : " << n_evt_ge2;
    if(n_evt>0) std::cout << "  (" << 100.0*double(n_evt_ge2)/double(n_evt) << " %)";
    std::cout << "\n";
  }else{
    std::cout << "[INFO] Filled all hits (non-unique) into hit-pattern.\n";
  }

  // Draw
  TCanvas* c1 = new TCanvas("cHTOF_tile34","HTOF 34-tile Hit Pattern",900,500);
  hTile->SetLineWidth(2);
  hTile->Draw("hist");

  TCanvas* c2 = nullptr;
  if(uniquePerEvent){
    c2 = new TCanvas("cHTOF_mult","HTOF Multiplicity",700,450);
    hMult->SetLineWidth(2);
    hMult->Draw("hist");
  }

  if(save){
    c1->SaveAs("HTOF_hitpattern_copyNo.png");
    if(c2) c2->SaveAs("HTOF_multiplicity_copyNo.png");
  }
}

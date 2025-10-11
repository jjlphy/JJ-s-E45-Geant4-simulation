// HTOF_BH2_Gated_Summary.C
// 1) 전체 이벤트 수
// 2) BH2 통과 이벤트 수
// 3) BH2 통과 중 HTOF >=1 이벤트 수
// 4) BH2 통과 중 HTOF >=2 이벤트 수
// 5) BH2 통과 이벤트의 HTOF 34-타일(0..33) + NONE(=34) 히트패턴 (유일 세그 기준)
//
// 사용법:
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
#include <vector>
#include <set>
#include <iostream>
#include <iomanip>

static const int kNTiles = 34;   // 0..33
static const int kBIN_NONE = 34; // 34번째 bin: NONE

struct _DICT_ {
  _DICT_(){
    gSystem->Load("libPhysics");
    gInterpreter->GenerateDictionary("vector<TParticle>","TParticle.h;vector");
  }
} _dict_;

static bool HasBranch(TTree* tr, const char* bname){
  return tr && tr->GetBranch(bname);
}

void HTOF_BH2_Gated_Summary(const char* filename="E45.root",
                            const char* treename="g4hyptpc",
                            bool save=true)
{
  // 파일/트리
  TFile* f = TFile::Open(filename,"READ");
  if(!f || f->IsZombie()){ std::cerr<<"[ERR] open "<<filename<<" failed\n"; return; }
  TTree* T = (TTree*)f->Get(treename);
  if(!T){ std::cerr<<"[ERR] tree "<<treename<<" not found\n"; return; }

  // 브랜치 준비
  const bool hasBH2 = HasBranch(T,"BH2");
  const bool hasHTOFvec = HasBranch(T,"HTOF");
  const bool hasHTOFcopy = HasBranch(T,"HTOF_copyNo"); // 있으면 우선 사용

  if(!hasBH2){ std::cerr<<"[ERR] need branch 'BH2' (vector<TParticle>)\n"; return; }
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

  // 히스토그램 (0..33 + NONE)
  TH1D* hHitPat = new TH1D("hHTOF_hitpattern_BH2",
     "BH2-gated HTOF HitPattern (unique segments/event);Tile ID;Events",
     kNTiles+1, -0.5, kNTiles+0.5);
  hHitPat->SetDirectory(nullptr);
  hHitPat->Sumw2();
  hHitPat->GetXaxis()->SetBinLabel(kNTiles+1, "NONE"); // bin index: 35th bin label

  // 카운터
  const Long64_t Ntot = T->GetEntries();
  Long64_t N_bh2=0, N_bh2_htof_ge1=0, N_bh2_htof_ge2=0;

  // 루프
  for(Long64_t ie=0; ie<Ntot; ++ie){
    T->GetEntry(ie);
    if(!BH2) continue;

    // 1) BH2 게이트
    if(BH2->empty()) continue;
    N_bh2++;

    // 2) 이 이벤트의 HTOF 유일 세그먼트 집합
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

    // 3) 통계/히스토그램
    if(uniq.empty()){
      hHitPat->Fill(kBIN_NONE); // NONE
    }else{
      N_bh2_htof_ge1++;
      if((int)uniq.size()>=2) N_bh2_htof_ge2++;
      for(int cn : uniq) hHitPat->Fill(cn); // 유일 세그 기준으로 채움
    }
  }

  // 리포트
  std::cout<<std::fixed<<std::setprecision(2);
  std::cout<<"========== BH2-gated HTOF Summary ==========\n";
  std::cout<<"Total events                      : "<<Ntot<<"\n";
  std::cout<<"BH2 pass                          : "<<N_bh2
           <<" ("<<(Ntot?100.0*double(N_bh2)/Ntot:0.0)<<" %)\n";
  std::cout<<"BH2 pass & HTOF >= 1              : "<<N_bh2_htof_ge1
           <<" ("<<(N_bh2?100.0*double(N_bh2_htof_ge1)/N_bh2:0.0)<<" % of BH2)\n";
  std::cout<<"BH2 pass & HTOF >= 2              : "<<N_bh2_htof_ge2
           <<" ("<<(N_bh2?100.0*double(N_bh2_htof_ge2)/N_bh2:0.0)<<" % of BH2)\n";
  std::cout<<"============================================\n";

  // 그림
  TCanvas* c = new TCanvas("cBH2_HTOF","BH2-gated HTOF Summary", 900, 500);
  hHitPat->SetLineWidth(2);
  hHitPat->Draw("hist");

  if(save){
    c->SaveAs("BH2_gated_HTOF_HitPattern.png");
  }
}

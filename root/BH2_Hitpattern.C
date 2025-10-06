// BH2_HitPattern_TP_v2.C  (ROOT 6.36 호환, vector<TParticle> 딕셔너리 자동 생성)
// - BH2 : std::vector<TParticle>
// - segField: "index" | "status" | "uniqueID" | "mother0" | "mother1" | "pdg"
// - valField: "E" | "Px" | "Py" | "Pz" | "Vx" | "Vy" | "Vz" | "Weight" | "Pt"
//root -l
//root [0] .L BH2_HitPattern_TP_v2.C+
//root [1] // 세그먼트=벡터 인덱스, 값=Energy(MeV), 컷=0.10
//root [1] BH2_HitPattern_TP_v2("../E45_with_SCH.root","g4hyptpc",0.10,"index","E",true);


#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <iomanip>
#include "TFile.h"
#include "TTree.h"
#include "TParticle.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TPaveText.h"
#include "TLatex.h"
#include "TSystem.h"
#include "TInterpreter.h"

static const int kNBH2 = 15;

static inline int getSegByField(const TParticle& p, const std::string& segField){
  if(segField=="status")    return p.GetStatusCode();
  if(segField=="uniqueID")  return p.GetUniqueID();
  if(segField=="mother0")   return p.GetMother(0);
  if(segField=="mother1")   return p.GetMother(1);
  if(segField=="pdg")       return p.GetPdgCode();
  // "index"는 외부에서 처리
  return p.GetStatusCode();
}

static inline double getVal(const TParticle& p, const std::string& valField){
  if(valField=="E")      return p.Energy();
  if(valField=="Px")     return p.Px();
  if(valField=="Py")     return p.Py();
  if(valField=="Pz")     return p.Pz();
  if(valField=="Vx")     return p.Vx();
  if(valField=="Vy")     return p.Vy();
  if(valField=="Vz")     return p.Vz();
  if(valField=="Weight") return p.GetWeight();
  if(valField=="Pt")     return p.Pt();
  return p.Energy();
}

void BH2_HitPattern_TP_v2(const char* filename,
                          const char* treename   = "g4hyptpc",
                          double cutValue        = 0.10,
                          const char* segFieldC  = "index",
                          const char* valFieldC  = "E",
                          bool normalizeToDenom  = true)
{
  // 딕셔너리 준비
  gSystem->Load("libEG");
  gInterpreter->GenerateDictionary("vector<TParticle>", "TParticle.h;vector");

  TFile* fin = TFile::Open(filename,"READ");
  if(!fin || fin->IsZombie()){ std::cerr<<"[ERR] open file\n"; return; }
  TTree* tr = (TTree*)fin->Get(treename);
  if(!tr){ std::cerr<<"[ERR] no tree\n"; fin->Close(); return; }

  std::vector<TParticle>* BH2 = nullptr;
  if(!tr->GetBranch("BH2")){ std::cerr<<"[ERR] no branch 'BH2'\n"; fin->Close(); return; }
  tr->SetBranchAddress("BH2",&BH2);

  std::string segField(segFieldC), valField(valFieldC);

  // --- 자동 판정: segField가 index가 아닐 때, 첫 이벤트에서 유효 세그가 하나도 안 잡히면 index로 전환
  bool autoSwitchedToIndex = false;
  if(segField!="index"){
    Long64_t nEnt = tr->GetEntries();
    if(nEnt>0){
      tr->GetEntry(0);
      bool anyValid = false;
      if(BH2){
        for(size_t i=0;i<BH2->size();++i){
          int seg = getSegByField(BH2->at(i), segField);
          if(0<=seg && seg<kNBH2){ anyValid = true; break; }
        }
      }
      if(!anyValid){
        std::cerr << "[WARN] segField='"<<segField<<"' 로는 유효 세그가 보이지 않습니다. 'index' 모드로 전환합니다.\n";
        segField = "index";
        autoSwitchedToIndex = true;
      }
    }
  }

  TH1F* h = new TH1F("hBH2","BH2 Hit Pattern;BH2 Segment ID;Counts", kNBH2, -0.5, kNBH2-0.5);
  h->SetStats(false);

  Long64_t nEnt = tr->GetEntries();
  Long64_t denom = 0;
  std::vector<Long64_t> segCounts(kNBH2,0);

  for(Long64_t ie=0; ie<nEnt; ++ie){
    tr->GetEntry(ie);
    if(!BH2) continue;

    bool pass = false;
    // 분모: 어떤 세그라도 cutValue 이상이면 BH2 pass
    for(size_t i=0;i<BH2->size();++i){
      int seg = (segField=="index") ? (int)i : getSegByField(BH2->at(i), segField);
      if(seg<0 || seg>=kNBH2) continue;
      double v = getVal(BH2->at(i), valField);
      if(v >= cutValue){ pass = true; break; }
    }
    if(!pass) continue;
    ++denom;

    // 세그 카운트
    for(size_t i=0;i<BH2->size();++i){
      int seg = (segField=="index") ? (int)i : getSegByField(BH2->at(i), segField);
      if(seg<0 || seg>=kNBH2) continue;
      double v = getVal(BH2->at(i), valField);
      if(v >= cutValue) ++segCounts[seg];
    }
  }

  for(int s=0;s<kNBH2;++s) h->SetBinContent(s+1, segCounts[s]);

  gStyle->SetOptStat(0);
  TCanvas* c = new TCanvas("cBH2TP","BH2 Hit Pattern (TParticle)", 900, 600);
  c->SetMargin(0.10,0.05,0.12,0.12);
  h->SetLineWidth(2);
  h->Draw("HIST");

  if(normalizeToDenom && denom>0){
    double maxCnt = h->GetMaximum();
    h->SetMaximum(maxCnt*1.25);
    TLatex tl; tl.SetTextSize(0.030);
    for(int s=0; s<kNBH2; ++s){
      double perc = 100.0 * segCounts[s] / (double)denom;
      double x = h->GetBinCenter(s+1);
      double y = segCounts[s] + 0.02*maxCnt;
      tl.DrawLatex(x-0.35, y, Form("%.2f%%", perc));
    }
  }

  TPaveText* pt = new TPaveText(0.60,0.60,0.94,0.94,"NDC");
  pt->SetFillColor(0);
  pt->SetTextAlign(12);
  pt->SetTextSize(0.032);
  pt->AddText(Form("File : %s", filename));
  pt->AddText(Form("Tree : %s", treename));
  pt->AddText("BH2 = vector<TParticle>");
  pt->AddText(Form("seg=%s%s, val=%s, cut=%.3f",
                   segField.c_str(),
                   (autoSwitchedToIndex?" (auto)":""),
                   valField.c_str(), cutValue));
  pt->AddText(Form("BH2-passing events: %lld (%.3f%%)",
                   denom, (nEnt>0? 100.0*denom/(double)nEnt:0.0)));
  pt->Draw();

  // 콘솔 요약
  std::cout << "========== BH2 Hit Pattern (vector<TParticle>) ==========\n";
  std::cout << "seg="<<segField<<(autoSwitchedToIndex?" (auto)":"")
            << ", val="<<valField<<", cut="<<cutValue<<"\n";
  std::cout << "Total entries: " << nEnt << "\n";
  std::cout << "BH2-passing(denom): " << denom
            << " (" << std::fixed<<std::setprecision(3)
            << (nEnt>0?100.0*denom/(double)nEnt:0.0) << "%)\n";
  std::cout << "---------------------------------------------------------\n";
  std::cout << " Seg | Counts";
  if(normalizeToDenom) std::cout << " | % of BH2-pass";
  std::cout << "\n";
  for(int s=0; s<kNBH2; ++s){
    double perc = (denom>0)? 100.0*segCounts[s]/(double)denom : 0.0;
    std::cout << std::setw(3)<<s<<" | "<<std::setw(8)<<segCounts[s];
    if(normalizeToDenom) std::cout << " | " << std::setw(8)<<std::fixed<<std::setprecision(3)<<perc;
    std::cout << "\n";
  }
  std::cout << "=========================================================\n";
}

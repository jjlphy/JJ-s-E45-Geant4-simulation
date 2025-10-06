// ProbeBH2.C  (ROOT 6.36 호환)
// BH2 브랜치(std::vector<TParticle>)의 필드 내용을 몇 이벤트 덤프
#include <iostream>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TParticle.h"

void ProbeBH2(const char* filename="../E45_with_SCH.root",
              const char* treename="g4hyptpc",
              int nEventsToScan=3)
{
  TFile* f = TFile::Open(filename, "READ");
  if(!f || f->IsZombie()){ std::cerr<<"[ERR] open file\n"; return; }
  TTree* T = (TTree*)f->Get(treename);
  if(!T){ std::cerr<<"[ERR] no tree\n"; return; }

  auto br = T->GetBranch("BH2");
  if(!br){ std::cerr<<"[ERR] no branch 'BH2'\n"; return; }
  std::cout << "[INFO] BH2 classname: " << br->GetClassName() << "\n";

  std::vector<TParticle>* v = nullptr;
  T->SetBranchAddress("BH2", &v);

  Long64_t nEnt = T->GetEntries();
  int nscan = std::min<Long64_t>(nEventsToScan, nEnt);
  for(int ie=0; ie<nscan; ++ie){
    T->GetEntry(ie);
    if(!v){ std::cout<<"[WARN] null vector at "<<ie<<"\n"; continue; }
    std::cout << "---- Event " << ie << " : size=" << v->size() << " ----\n";
    for(size_t i=0; i<v->size() && i<10; ++i){
      const TParticle& p = v->at(i);
      std::cout << "  ["<<i<<"]"
                << " PDG="      << p.GetPdgCode()
                << " status="   << p.GetStatusCode()
                << " mother=("  << p.GetMother(0) << "," << p.GetMother(1) << ")"
                << " uniqueID=" << p.GetUniqueID()
                << " E(MeV)="   << p.Energy()
                << " P=("       << p.Px() << "," << p.Py() << "," << p.Pz() << ")"
                << " Vtx=("     << p.Vx() << "," << p.Vy() << "," << p.Vz() << ")"
                << "\n";
    }
  }
  f->Close();
}


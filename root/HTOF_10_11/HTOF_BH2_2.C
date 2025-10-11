// HTOF_BH2_2.C
//
// (1) BH2-gated, HTOF tiles {0..5} vs BH2 segment(0..14) 2D 히스토
// (2) BH2-gated, HTOF(0..33) vs BH2(0..14) 2D 히스토
//
// 사용법:
//   root -l
//   .L HTOF_BH2_2.C+
//   HTOF_BH2_2("E45.root","g4hyptpc", /*save=*/true);

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH2I.h"
#include "TH1D.h"
#include "TSystem.h"
#include "TInterpreter.h"
#include "TParticle.h"
#include "TString.h"
#include <vector>
#include <set>
#include <cmath>
#include <iostream>

// ---------------------- geometry assumptions (E72) ----------------------
// BH2 위치/방향은 E72 파라미터(코드에서 RA/TA=0) 기준.
// BH2는 X 방향으로 15개 세그먼트, 각 폭 14 mm (전 구간 동일 폭).
// 월드좌표에서 BH2 중심은 (-10, 0, -560) mm 로 놓고 회전은 0으로 가정.
static const int    kNBH2Seg   = 15;
static const double kBH2_x0    = -10.0;  // [mm] BH2 center (world)
static const double kBH2_y0    =   0.0;
static const double kBH2_z0    = -560.0;
static const double kBH2_pitch = 14.0;   // [mm]
static const double kBH2_total = kNBH2Seg * kBH2_pitch; // 210 mm
// BH2 mother local에서 세그 i의 중심은 x = -total/2 + (i+0.5)*pitch
// 월드좌표와 동축 가정이므로 x' = x_world - kBH2_x0 로 사용.

static const int    kNHTOF = 34; // tiles 0..33

// vector<TParticle> 사전 준비
struct _DICT_ {
  _DICT_(){
    gSystem->Load("libPhysics");
    gInterpreter->GenerateDictionary("vector<TParticle>","TParticle.h;vector");
  }
} _dict_;

static bool HasBranch(TTree* tr, const char* bname){
  return tr && tr->GetBranch(bname);
}

// 월드좌표 (x,y,z) → BH2 세그먼트 인덱스(0..14), 실패 시 -1
static int MapBH2_WorldToSeg(double x, double y, double z)
{
  // (선택) z 근접성 체크를 하고 싶으면 아래 주석 해제 (허용 오차 10 mm 예시)
  // if (std::fabs(z - kBH2_z0) > 10.0) return -1;

  const double xloc = x - kBH2_x0; // BH2 중심 기준 x'
  const double xmin = -0.5*kBH2_total;           // -105
  const double xmax =  0.5*kBH2_total;           // +105
  if (xloc < xmin || xloc >= xmax) return -1;

  // 왼쪽 끝에서부터 pitch 간격으로 나눔
  const double rel = xloc - xmin;                // [0, total)
  int idx = int(std::floor(rel / kBH2_pitch));   // 0..14
  if (idx < 0) idx = 0;
  if (idx >= kNBH2Seg) idx = kNBH2Seg-1;
  return idx;
}

void HTOF_BH2_2(const char* filename="E45.root",
                const char* treename="g4hyptpc",
                bool save=true)
{
  // 파일/트리
  TFile* f = TFile::Open(filename,"READ");
  if(!f || f->IsZombie()){ std::cerr<<"[ERR] open "<<filename<<" failed\n"; return; }
  TTree* T = (TTree*)f->Get(treename);
  if(!T){ std::cerr<<"[ERR] tree "<<treename<<" not found\n"; return; }

  // 브랜치 확인
  if(!HasBranch(T,"BH2")){ std::cerr<<"[ERR] need branch 'BH2' (vector<TParticle>)\n"; return; }
  const bool hasHTOFvec  = HasBranch(T,"HTOF");
  const bool hasHTOFcopy = HasBranch(T,"HTOF_copyNo"); // 있으면 우선 사용

  if(!(hasHTOFvec || hasHTOFcopy)){
    std::cerr<<"[ERR] need 'HTOF' (vector<TParticle>) or 'HTOF_copyNo' (vector<int>)\n";
    return;
  }

  // 브랜치 포인터
  std::vector<TParticle>* BH2 = nullptr;
  std::vector<TParticle>* HTOF = nullptr;
  std::vector<int>* HTOF_copyNo = nullptr;

  T->SetBranchAddress("BH2",&BH2);
  if(hasHTOFvec)  T->SetBranchAddress("HTOF",&HTOF);
  if(hasHTOFcopy) T->SetBranchAddress("HTOF_copyNo",&HTOF_copyNo);

  // 히스토그램
  // (1) HTOF {0..5} vs BH2(0..14)
  TH2I* h05_vs_BH2 = new TH2I("hHTOF05_vs_BH2",
      "BH2-gated: HTOF tiles {0..5} vs BH2 segment;HTOF tile (0..5);BH2 seg (0..14)",
      6, -0.5, 5.5,   // x: 0..5
      kNBH2Seg, -0.5, kNBH2Seg-0.5); // y: 0..14
  h05_vs_BH2->SetDirectory(nullptr);

  // (2) 전체 HTOF(0..33) vs BH2(0..14)
  TH2I* hAll_vs_BH2 = new TH2I("hHTOF_vs_BH2",
      "BH2-gated: HTOF (0..33) vs BH2 (0..14);HTOF tile;BH2 seg",
      kNHTOF, -0.5, kNHTOF-0.5,
      kNBH2Seg, -0.5, kNBH2Seg-0.5);
  hAll_vs_BH2->SetDirectory(nullptr);

  // 루프
  const Long64_t N = T->GetEntries();
  for(Long64_t ie=0; ie<N; ++ie){
    T->GetEntry(ie);
    if(!BH2) continue;

    // BH2 게이트: BH2 비어있으면 스킵
    if(BH2->empty()) continue;

    // 이벤트 내 BH2 유일 세그 집합 (월드좌표 → seg 매핑)
    std::set<int> bh2Seg;
    for(const auto& p : *BH2){
      int s = MapBH2_WorldToSeg(p.Vx(), p.Vy(), p.Vz());
      if(0 <= s && s < kNBH2Seg) bh2Seg.insert(s);
    }
    if(bh2Seg.empty()) continue; // 매핑 실패 시 스킵

    // 이벤트 내 HTOF 유일 세그 집합
    std::set<int> htofSeg;
    if(hasHTOFcopy && HTOF_copyNo){
      for(int cn : *HTOF_copyNo){
        if(0 <= cn && cn < kNHTOF) htofSeg.insert(cn);
      }
    }else if(hasHTOFvec && HTOF){
      for(const auto& p : *HTOF){
        int cn = p.GetStatusCode(); // 0..33
        if(0 <= cn && cn < kNHTOF) htofSeg.insert(cn);
      }
    }

    // (1) HTOF {0..5} vs BH2
    for(int t : htofSeg){
      if(0 <= t && t <= 5){
        for(int s : bh2Seg){
          h05_vs_BH2->Fill(t, s);
        }
      }
    }

    // (2) 전체 HTOF vs BH2
    for(int t : htofSeg){
      for(int s : bh2Seg){
        hAll_vs_BH2->Fill(t, s);
      }
    }
  }

  // 그리기
  TCanvas* c1 = new TCanvas("cHTOF05_vs_BH2","HTOF{0..5} vs BH2", 800, 600);
  h05_vs_BH2->Draw("colz");
  c1->SetRightMargin(0.15);

  TCanvas* c2 = new TCanvas("cHTOF_vs_BH2","HTOF (0..33) vs BH2 (0..14)", 900, 650);
  hAll_vs_BH2->Draw("colz");
  c2->SetRightMargin(0.15);

  if(save){
    c1->SaveAs("HTOF05_vs_BH2.png");
    c2->SaveAs("HTOF_vs_BH2.png");
  }
}

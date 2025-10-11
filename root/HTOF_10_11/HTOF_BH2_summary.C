// HTOF_BH2_summary.C
//
// 사용법:
//   root -l
//   .L HTOF_BH2_summary.C+
//   HTOF_BH2_summary("E45.root","g4hyptpc", /*save=*/true);
//
// 기능:
//  <Section1> 전체 BH2 게이트 기준
//     - 전체 이벤트 수
//     - BH2 통과 이벤트 수 (전체 대비 %)
//     - BH2 & HTOF>=1 (전체 대비 %, BH2 대비 %)
//     - BH2 & HTOF>=2 (전체 대비 %, BH2 대비 %, BH2&HTOF>=1 대비 %)
//     - BH2 & HTOF==2 (전체 대비 %, BH2 대비 %, BH2&HTOF>=1 대비 %)
//     - HTOF==2인 이벤트가 (18,19),(19,20),...,(24,25) 어느 페어에 속하는지 요약
//
//  <Section2> BH2 세그먼트 3~10에 '하나라도' 닿은 이벤트로 동일 집계/출력
//  <Section3> BH2 세그먼트 4~9 에 '하나라도' 닿은 이벤트로 동일 집계/출력
//
//  히스토그램:
//    - HTOF 패턴(0..33) 1D: BH2(전체), BH2 seg 3-10, BH2 seg 4-9  → PNG 저장
//
// 브랜치 요구:
//    BH2 : vector<TParticle>   (Vx,Vy,Vz를 사용해서 BH2 세그먼트 매핑)
//    HTOF_copyNo : vector<int> (있으면 우선 사용)  또는
//    HTOF : vector<TParticle>  (StatusCode에 0..33 저장되어 있다고 가정)
//
// 주의:
//   - BH2 기하학 가정(E72 스타일): 중심(-10,0,-560)mm, X방향 15세그, pitch=14mm
//   - 필요하면 kBH2_x0 등 상수는 현 설정에 맞게 수정

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1I.h"
#include "TSystem.h"
#include "TInterpreter.h"
#include "TParticle.h"
#include "TString.h"
#include <vector>
#include <set>
#include <map>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <algorithm>

// ---------------------- geometry assumptions (E72-like) ----------------------
static const int    kNBH2Seg   = 15;     // 0..14
static const double kBH2_x0    = -10.0;  // [mm] world center
static const double kBH2_y0    =   0.0;
static const double kBH2_z0    = -560.0;
static const double kBH2_pitch = 14.0;   // [mm]
static const double kBH2_total = kNBH2Seg * kBH2_pitch; // 210 mm
// 세그 i 중심 (BH2 local x) = -total/2 + (i+0.5)*pitch (회전 0 가정)

static const int    kNHTOF     = 34;     // tiles 0..33

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

// 월드좌표 → BH2 세그 인덱스(0..14), 실패 시 -1
static int MapBH2_WorldToSeg(double x, double y, double z)
{
  // z 근접성 체크를 넣고 싶으면 활성화 (예: |z - kBH2_z0| <= 10 mm)
  // if (std::fabs(z - kBH2_z0) > 10.0) return -1;

  const double xloc = x - kBH2_x0;           // BH2 중심 기준 x'
  const double xmin = -0.5*kBH2_total;       // -105
  const double xmax =  0.5*kBH2_total;       // +105
  if (xloc < xmin || xloc >= xmax) return -1;

  const double rel = xloc - xmin;            // [0,total)
  int idx = int(std::floor(rel / kBH2_pitch)); // 0..14
  if (idx < 0) idx = 0;
  if (idx >= kNBH2Seg) idx = kNBH2Seg-1;
  return idx;
}

// 비율 출력용 헬퍼
static inline double pct(double num, double den){ return (den>0)? (100.0*num/den):0.0; }

// HTOF==2인 이벤트에서 (i,i+1) 페어 카운팅: i∈{18..24}
struct PairCount {
  std::map<std::pair<int,int>, long long> counts;
  long long others = 0;
  PairCount(){
    for(int i=18;i<=24;++i) counts[{i,i+1}] = 0;
  }
  void Fill(const std::set<int>& tiles){
    if((int)tiles.size()!=2){ return; }
    auto it = tiles.begin();
    int a = *it; ++it; int b = *it;
    if(a>b) std::swap(a,b);
    auto key = std::make_pair(a,b);
    auto f = counts.find(key);
    if(f!=counts.end()) f->second++;
    else others++;
  }
  void PrintSummary(const char* tag){
    std::cout << "  [" << tag << "] HTOF==2 pair summary:\n";
    long long total = 0;
    for(auto &kv : counts) total += kv.second;
    total += others;
    for(auto &kv : counts){
      std::cout << "    pair(" << kv.first.first << "," << kv.first.second << ") : "
                << kv.second << "\n";
    }
    std::cout << "    others : " << others << "\n";
    std::cout << "    total(HTOF==2) counted in summary : " << total << "\n";
  }
};

// 집계 박스(조건별로 한 번씩 사용)
struct StatsBox {
  // 전역 기준
  long long N_total = 0;

  // 게이트(조건)에 해당하는 BH2-선택 이벤트 수
  long long N_BH2sel = 0;

  // HTOF multiplicity 기준
  long long N_ge1 = 0;
  long long N_ge2 = 0;
  long long N_eq2 = 0;

  PairCount pairSummary;

  // HTOF 히스토(0..33) : 유니크 타일 기준으로 1 이벤트당 1카운트씩
  TH1I* hHTOF = nullptr;

  // 초기화
  void InitHist(const char* hname, const char* htitle){
    hHTOF = new TH1I(hname, htitle, kNHTOF, -0.5, kNHTOF-0.5);
    hHTOF->SetDirectory(nullptr);
  }

  // 결과 출력
  void Print(const char* title){
    std::cout << "\n<" << title << ">\n";
    std::cout << "  Total events                            : " << N_total << "\n";
    std::cout << "  BH2-selected events                     : " << N_BH2sel
              << "  (" << std::fixed << std::setprecision(3) << pct(N_BH2sel, N_total) << " % of Total)\n";

    std::cout << "  BH2 & HTOF>=1                           : " << N_ge1
              << "  (" << pct(N_ge1, N_total)  << " % of Total, "
              << pct(N_ge1, N_BH2sel)        << " % of BH2)\n";

    std::cout << "  BH2 & HTOF>=2                           : " << N_ge2
              << "  (" << pct(N_ge2, N_total)  << " % of Total, "
              << pct(N_ge2, N_BH2sel)        << " % of BH2, "
              << pct(N_ge2, N_ge1)           << " % of [BH2 & HTOF>=1])\n";

    std::cout << "  BH2 & HTOF==2                           : " << N_eq2
              << "  (" << pct(N_eq2, N_total)  << " % of Total, "
              << pct(N_eq2, N_BH2sel)        << " % of BH2, "
              << pct(N_eq2, N_ge1)           << " % of [BH2 & HTOF>=1])\n";

    pairSummary.PrintSummary(title);
  }
};

// BH2 세그 집합이 [lo..hi] 와 교집합이 있는지
static bool BH2HitsInRange(const std::set<int>& bh2Seg, int lo, int hi){
  for(int s : bh2Seg){
    if(lo <= s && s <= hi) return true;
  }
  return false;
}

void HTOF_BH2_summary(const char* filename="E45.root",
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
  const bool hasHTOFcopy = HasBranch(T,"HTOF_copyNo");

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

  // 섹션별 집계 박스 & 히스토
  StatsBox S1, S2, S3;
  S1.InitHist("hHTOF_BH2_all",   "HTOF pattern | BH2 (any seg);HTOF tile;Events");
  S2.InitHist("hHTOF_BH2_3_10",  "HTOF pattern | BH2 seg 3-10;HTOF tile;Events");
  S3.InitHist("hHTOF_BH2_4_9",   "HTOF pattern | BH2 seg 4-9;HTOF tile;Events");

  // 루프
  const Long64_t N = T->GetEntries();
  for(Long64_t ie=0; ie<N; ++ie){
    T->GetEntry(ie);
    S1.N_total++; S2.N_total++; S3.N_total++;

    // (1) BH2 세그 집합
    std::set<int> bh2Seg;
    if(BH2){
      for(const auto& p : *BH2){
        int s = MapBH2_WorldToSeg(p.Vx(), p.Vy(), p.Vz());
        if(0 <= s && s < kNBH2Seg) bh2Seg.insert(s);
      }
    }
    const bool passBH2_any   = !bh2Seg.empty();
    const bool passBH2_3_10  = BH2HitsInRange(bh2Seg, 3, 10);
    const bool passBH2_4_9   = BH2HitsInRange(bh2Seg, 4, 9);

    // (2) HTOF 세그 집합 (유니크)
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
    const int htMult = (int)htofSeg.size();

    // ---------- Section1: BH2(any) ----------
    if(passBH2_any){
      S1.N_BH2sel++;
      if(htMult>=1) S1.N_ge1++;
      if(htMult>=2) S1.N_ge2++;
      if(htMult==2){
        S1.N_eq2++;
        S1.pairSummary.Fill(htofSeg);
      }
      // HTOF 히스토그램(유니크 타일당 1카운트)
      for(int t : htofSeg) S1.hHTOF->Fill(t);
    }

    // ---------- Section2: BH2 seg 3~10 ----------
    if(passBH2_3_10){
      S2.N_BH2sel++;
      if(htMult>=1) S2.N_ge1++;
      if(htMult>=2) S2.N_ge2++;
      if(htMult==2){
        S2.N_eq2++;
        S2.pairSummary.Fill(htofSeg);
      }
      for(int t : htofSeg) S2.hHTOF->Fill(t);
    }

    // ---------- Section3: BH2 seg 4~9 ----------
    if(passBH2_4_9){
      S3.N_BH2sel++;
      if(htMult>=1) S3.N_ge1++;
      if(htMult>=2) S3.N_ge2++;
      if(htMult==2){
        S3.N_eq2++;
        S3.pairSummary.Fill(htofSeg);
      }
      for(int t : htofSeg) S3.hHTOF->Fill(t);
    }
  }

  // 결과 출력
  S1.Print("Section1: BH2(any seg)");
  S2.Print("Section2: BH2 seg 3-10");
  S3.Print("Section3: BH2 seg 4-9");

  // 히스토그램 그리기 & 저장
  TCanvas* c1 = new TCanvas("cHTOF_BH2_all","HTOF | BH2(any)",800,600);
  S1.hHTOF->Draw("hist");
  TCanvas* c2 = new TCanvas("cHTOF_BH2_3_10","HTOF | BH2(3-10)",800,600);
  S2.hHTOF->Draw("hist");
  TCanvas* c3 = new TCanvas("cHTOF_BH2_4_9","HTOF | BH2(4-9)",800,600);
  S3.hHTOF->Draw("hist");

  if(save){
    c1->SaveAs("HTOF_pattern_BH2_any.png");
    c2->SaveAs("HTOF_pattern_BH2_3_10.png");
    c3->SaveAs("HTOF_pattern_BH2_4_9.png");
  }
}

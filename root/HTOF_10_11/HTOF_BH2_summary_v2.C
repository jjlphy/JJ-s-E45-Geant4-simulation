// HTOF_BH2_summary_v2.C  (namespace 격리판 + S2/S3 페어 히스토 + 퍼센트 출력)
// 사용법:
//   root -l
//   .L HTOF_BH2_summary_v2.C+
//   HTOF_BH2_summary_v2("../rootfile/E45_Beam1.root","g4hyptpc",true);

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

namespace HBS { // ======== 모든 전역/헬퍼 격리 ========

// ---------------------- geometry assumptions (E72-like) ----------------------
static const int    kNBH2Seg   = 15;     // 0..14
static const double kBH2_x0    = -10.0;  // [mm] world center
static const double kBH2_y0    =   0.0;
static const double kBH2_z0    = -560.0;
static const double kBH2_pitch = 14.0;   // [mm]
static const double kBH2_total = kNBH2Seg * kBH2_pitch; // 210 mm

static const int    kNHTOF     = 34;     // tiles 0..33

// vector<TParticle> 사전 준비 (고유 이름)
struct DictGuard_HTOFSUM {
  DictGuard_HTOFSUM(){
    gSystem->Load("libPhysics");
    gInterpreter->GenerateDictionary("vector<TParticle>","TParticle.h;vector");
  }
};
static DictGuard_HTOFSUM _dict_guard_;

static bool HasBranch(TTree* tr, const char* bname){
  return tr && tr->GetBranch(bname);
}

// 월드좌표 → BH2 세그(0..14), 실패 시 -1
static int MapBH2_WorldToSeg(double x, double y, double z)
{
  // if (std::fabs(z - kBH2_z0) > 10.0) return -1; // z근접성 쓰고싶으면 해제
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

// 비율 계산
static inline double pct(double num, double den){ return (den>0)? (100.0*num/den):0.0; }

// HTOF==2 페어 집계 (특수 7페어 + others)
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
  void PrintSummary(const char* tag, long long N_ge2){
    std::cout << "  [" << tag << "] HTOF==2 pair summary:\n";
    long long total = 0;
    for(auto &kv : counts) total += kv.second;
    total += others;
    for(auto &kv : counts){
      double p = pct(kv.second, N_ge2);
      std::cout << "    pair(" << kv.first.first << "," << kv.first.second << ") : "
                << kv.second << "  (" << std::fixed << std::setprecision(3) << p << " % of [BH2 & HTOF>=2])\n";
    }
    std::cout << "    others : " << others << "\n";
    std::cout << "    total(HTOF==2) counted in summary : " << total << "\n";
  }
};

// 집계 컨테이너
struct StatsBox {
  long long N_total = 0;   // 전체 이벤트
  long long N_BH2sel = 0;  // 조건 충족 BH2-선택 이벤트
  long long N_ge1 = 0;     // HTOF>=1
  long long N_ge2 = 0;     // HTOF>=2
  long long N_eq2 = 0;     // HTOF==2

  PairCount pairSummary;

  TH1I* hHTOF = nullptr;        // HTOF 1D 패턴(유니크 타일)
  TH1I* hPairAdj = nullptr;     // (옵션) 인접 페어 34개 (i,i+1), i=0..33(33→0)

  void InitHist(const char* hname, const char* htitle){
    hHTOF = new TH1I(hname, htitle, kNHTOF, -0.5, kNHTOF-0.5);
    hHTOF->SetDirectory(nullptr);
  }
  // 섹션2/3 전용: 인접 페어 히스토그램(34bin) 초기화 + 라벨
  void InitPairHist(const char* hname, const char* htitle){
    hPairAdj = new TH1I(hname, htitle, kNHTOF, -0.5, kNHTOF-0.5);
    hPairAdj->SetDirectory(nullptr);
    // bin label: "i-(i+1)" with wrap
    for(int i=0;i<kNHTOF;++i){
      int j = (i+1) % kNHTOF;
      TString lab; lab.Form("%d-%d", i, j);
      hPairAdj->GetXaxis()->SetBinLabel(i+1, lab);
    }
  }
  // HTOF==2일 때 인접페어면 해당 bin에 +1 (wrap 포함), 아니면 others만 증가 (히스토엔 미반영)
  void FillPairHistIfAdjacent(const std::set<int>& tiles){
    if(!hPairAdj) return;
    if((int)tiles.size()!=2) return;
    auto it = tiles.begin();
    int a = *it; ++it; int b = *it;
    if(a>b) std::swap(a,b);
    // 인접성 체크 (wrap 포함)
    if(b==a+1){
      hPairAdj->Fill(a);
    }else if(a==0 && b==(kNHTOF-1)){ // (0,33) → bin index 33
      hPairAdj->Fill(kNHTOF-1);
    }else{
      // 비인접 조합은 히스토그램에는 넣지 않음 (others는 pairSummary에서 집계됨)
    }
  }

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

    pairSummary.PrintSummary(title, N_ge2);
  }
};

static bool BH2HitsInRange(const std::set<int>& bh2Seg, int lo, int hi){
  for(int s : bh2Seg){
    if(lo <= s && s <= hi) return true;
  }
  return false;
}

} // namespace HBS

// ======== 공개 함수 (이름도 v2로 유지) ========
void HTOF_BH2_summary_v2(const char* filename="E45.root",
                         const char* treename="g4hyptpc",
                         bool save=true)
{
  using namespace HBS;

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
  // 섹션2/3: 인접 페어 히스토도 생성
  S2.InitPairHist("hPairs_BH2_3_10", "Adjacent HTOF pairs (i-(i+1)) | BH2 seg 3-10;pair;Events");
  S3.InitPairHist("hPairs_BH2_4_9",  "Adjacent HTOF pairs (i-(i+1)) | BH2 seg 4-9;pair;Events");

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
        S2.FillPairHistIfAdjacent(htofSeg); // 인접 페어 히스토 누적
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
        S3.FillPairHistIfAdjacent(htofSeg); // 인접 페어 히스토 누적
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

  // Section2/3 페어 히스토 (인접 페어 34bin)
  TCanvas* c4 = new TCanvas("cPairs_BH2_3_10","Adjacent pairs | BH2(3-10)",900,600);
  if(S2.hPairAdj){ S2.hPairAdj->LabelsOption("v","X"); S2.hPairAdj->Draw("hist"); }

  TCanvas* c5 = new TCanvas("cPairs_BH2_4_9","Adjacent pairs | BH2(4-9)",900,600);
  if(S3.hPairAdj){ S3.hPairAdj->LabelsOption("v","X"); S3.hPairAdj->Draw("hist"); }

  if(save){
    c1->SaveAs("HTOF_pattern_BH2_any.png");
    c2->SaveAs("HTOF_pattern_BH2_3_10.png");
    c3->SaveAs("HTOF_pattern_BH2_4_9.png");
    if(S2.hPairAdj) c4->SaveAs("HTOF_pairs_BH2_3_10.png");
    if(S3.hPairAdj) c5->SaveAs("HTOF_pairs_BH2_4_9.png");
  }
}

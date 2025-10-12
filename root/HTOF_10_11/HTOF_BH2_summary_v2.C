// HTOF_BH2_summary_v2.C  (구간 가변판: "4-9,4-11,3-8" 같은 문자열로 지정)
// 사용법:
//   root -l
//   .L HTOF_BH2_summary_v2.C+
//   HTOF_BH2_summary_v2("../rootfile/E45_Beam1.root","g4hyptpc","4-9,4-11,3-8",true);
//
// 설명:
//   - Section1은 항상 "BH2(any)"로 집계
//   - ranges 문자열에 적은 각 구간(예: 4-9)을 Section2,3,... 로 자동 생성
//   - 각 섹션에 대해 HTOF 1D 패턴, HTOF==2 인접 페어 히스토그램(i-(i+1)) 생성
//   - 콘솔에는 각 섹션의 카운트/퍼센트와 (18..24) 특수 페어 + others 요약 출력

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1I.h"
#include "TSystem.h"
#include "TInterpreter.h"
#include "TParticle.h"
#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"
#include <vector>
#include <set>
#include <map>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cctype>

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
  TH1I* hPairAdj = nullptr;     // 인접 페어 34개 (i,i+1), i=0..33(33→0)

  void InitHist(const TString& key, const TString& title){
    hHTOF = new TH1I(key, title, kNHTOF, -0.5, kNHTOF-0.5);
    hHTOF->SetDirectory(nullptr);
  }
  void InitPairHist(const TString& key, const TString& title){
    hPairAdj = new TH1I(key, title, kNHTOF, -0.5, kNHTOF-0.5);
    hPairAdj->SetDirectory(nullptr);
    for(int i=0;i<kNHTOF;++i){
      int j = (i+1) % kNHTOF;
      TString lab; lab.Form("%d-%d", i, j);
      hPairAdj->GetXaxis()->SetBinLabel(i+1, lab);
    }
  }
  void FillPairHistIfAdjacent(const std::set<int>& tiles){
    if(!hPairAdj) return;
    if((int)tiles.size()!=2) return;
    auto it = tiles.begin();
    int a = *it; ++it; int b = *it;
    if(a>b) std::swap(a,b);
    if(b==a+1){
      hPairAdj->Fill(a);
    }else if(a==0 && b==(kNHTOF-1)){
      hPairAdj->Fill(kNHTOF-1);
    }else{
      // 비인접 조합은 스킵 (others는 PairCount에서 집계)
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

// ===== 섹션 스펙과 파서 =====
struct SectionSpec {
  bool isAny = false; // true면 BH2(any)
  int  lo = 0;
  int  hi = 0;
  TString tag;   // "BH2(any)" 또는 "BH2 seg lo-hi"
  TString key;   // 캔버스/히스토 키 접두어
};

static inline TString TrimWS(const TString& s){
  TString t=s; t.ReplaceAll(" ",""); t.ReplaceAll("\t",""); return t;
}

// "4-9,4-11,3-8" → { (4,9), (4,11), (3,8) }
static std::vector<SectionSpec> ParseRanges(const char* rangesCSV){
  std::vector<SectionSpec> out;
  // 항상 Section1: any
  SectionSpec any; any.isAny=true; any.tag="BH2(any seg)"; any.key="BH2_any";
  out.push_back(any);

  if(!rangesCSV) return out;
  TString csv = rangesCSV;
  csv = csv.Strip(TString::kBoth);
  if(csv.Length()==0) return out;

  TObjArray* toks = csv.Tokenize(",");
  if(!toks) return out;
  for(int i=0;i<toks->GetEntriesFast(); ++i){
    auto* obj = toks->At(i);
    if(!obj) continue;
    TString tok = ((TObjString*)obj)->GetString();
    tok = TrimWS(tok);
    if(tok.Length()==0) continue;

    // 허용 형식: "a-b"
    Ssiz_t dash = tok.Index("-");
    if(dash==kNPOS){ std::cerr<<"[WARN] bad range token: "<<tok<<"\n"; continue; }
    TString A = tok(0, dash);
    TString B = tok(dash+1, tok.Length()-dash-1);
    A = TrimWS(A); B = TrimWS(B);
    if(A.Length()==0 || B.Length()==0){ std::cerr<<"[WARN] bad range token: "<<tok<<"\n"; continue; }

    int lo = A.Atoi();
    int hi = B.Atoi();
    if(lo>hi) std::swap(lo,hi);
    if(lo<0) lo=0;
    if(hi>=kNBH2Seg) hi=kNBH2Seg-1;

    SectionSpec sp;
    sp.isAny=false; sp.lo=lo; sp.hi=hi;
    sp.tag.Form("BH2 seg %d-%d", lo, hi);
    sp.key.Form("BH2_%d_%d", lo, hi);
    out.push_back(sp);
  }
  delete toks;
  return out;
}

} // namespace HBS

// ======== 공개 함수 ========
void HTOF_BH2_summary_v2(const char* filename="E45.root",
                         const char* treename="g4hyptpc",
                         const char* ranges="4-9,4-11,3-8",
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

  // ===== 섹션 스펙 구성 (Section1은 always any, 이후는 ranges 문자열 기반) =====
  std::vector<SectionSpec> specs = ParseRanges(ranges);

  // 섹션 별 StatsBox 및 히스토 준비
  struct SectionRun {
    SectionSpec spec;
    StatsBox box;
    TCanvas* cHTOF = nullptr;
    TCanvas* cPair = nullptr;
  };
  std::vector<SectionRun> runs;
  runs.reserve(specs.size());

  for(size_t i=0;i<specs.size(); ++i){
    const auto& sp = specs[i];
    SectionRun r; r.spec = sp;

    // HTOF 1D 히스토
    TString hname, htitle;
    hname.Form("hHTOF_%s", sp.key.Data());
    if(sp.isAny){
      htitle = "HTOF pattern | BH2(any seg);HTOF tile;Events";
    }else{
      htitle.Form("HTOF pattern | %s;HTOF tile;Events", sp.tag.Data());
    }
    r.box.InitHist(hname, htitle);

    // 인접 페어 히스토 (모든 섹션에 대해 생성)
    TString hp, tp;
    hp.Form("hPairs_%s", sp.key.Data());
    if(sp.isAny){
      tp = "Adjacent HTOF pairs (i-(i+1)) | BH2(any seg);pair;Events";
    }else{
      tp.Form("Adjacent HTOF pairs (i-(i+1)) | %s;pair;Events", sp.tag.Data());
    }
    r.box.InitPairHist(hp, tp);

    runs.push_back(r);
  }

  // ===== 이벤트 루프 =====
  const Long64_t N = T->GetEntries();
  for(Long64_t ie=0; ie<N; ++ie){
    T->GetEntry(ie);
    // 모든 섹션의 N_total 동시 증가
    for(auto& r : runs) r.box.N_total++;

    // (1) BH2 세그 집합
    std::set<int> bh2Seg;
    if(BH2){
      for(const auto& p : *BH2){
        int s = MapBH2_WorldToSeg(p.Vx(), p.Vy(), p.Vz());
        if(0 <= s && s < kNBH2Seg) bh2Seg.insert(s);
      }
    }
    const bool passBH2_any = !bh2Seg.empty();

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

    // 섹션별 집계
    for(auto& r : runs){
      bool passBH2 = false;
      if(r.spec.isAny){
        passBH2 = passBH2_any;
      }else{
        passBH2 = BH2HitsInRange(bh2Seg, r.spec.lo, r.spec.hi);
      }

      if(passBH2){
        r.box.N_BH2sel++;
        if(htMult>=1) r.box.N_ge1++;
        if(htMult>=2) r.box.N_ge2++;
        if(htMult==2){
          r.box.N_eq2++;
          r.box.pairSummary.Fill(htofSeg);
          r.box.FillPairHistIfAdjacent(htofSeg);
        }
        for(int t : htofSeg) r.box.hHTOF->Fill(t);
      }
    }
  }

  // ===== 출력 + 그림 =====
  for(size_t i=0;i<runs.size(); ++i){
    auto& r = runs[i];
    TString title = r.spec.tag.Length()? r.spec.tag : "BH2(any seg)";
    r.box.Print(title);

    // 1D 패턴
    TString cname1; cname1.Form("cHTOF_%s", r.spec.key.Data());
    r.cHTOF = new TCanvas(cname1, title, 800, 600);
    r.box.hHTOF->Draw("hist");

    // 인접 페어
    TString cname2; cname2.Form("cPairs_%s", r.spec.key.Data());
    r.cPair = new TCanvas(cname2, TString("Adjacent pairs | ")+title, 900, 600);
    r.box.hPairAdj->LabelsOption("v","X");
    r.box.hPairAdj->Draw("hist");

    if(save){
      TString fn1; fn1.Form("HTOF_pattern_%s.png", r.spec.key.Data());
      TString fn2; fn2.Form("HTOF_pairs_%s.png",   r.spec.key.Data());
      r.cHTOF->SaveAs(fn1);
      r.cPair->SaveAs(fn2);
    }
  }
}

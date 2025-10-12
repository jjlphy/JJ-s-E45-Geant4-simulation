// HTOF_BH2_summary_v2.C
// (가변 BH2 구간 + 0.1 MIP 컷 + 임베디드 edep/ID 자동 인식 버전)
// 사용법:
//   root -l
//   .L HTOF_BH2_summary_v2.C+
//   HTOF_BH2_summary_v2("../rootfile/E45_fix_Beam_098.root","g4hyptpc","3-11,4-9",true);
//
// 설명:
//   - Section1은 항상 BH2(any)로 집계, ranges 문자열("a-b,a-b,...")로 추가 섹션 자동 생성
//   - 0.1 MIP 컷 기본 적용 (폴리스티렌 dE/dx=2 MeV/cm, BH2=5 mm, HTOF=10 mm)
//   - edep 소스 자동 탐지 우선순위:
//        (1) *_edep 브랜치 (vector<double>) + ID는 Mother(1) 우선/StatusCode 보조
//        (2) TParticle 임베디드: Mother(1)=ID, Weight=edep  ← BVH 코드와 동일 규약
//        (3) (BH2만) 좌표→세그 매핑(에너지는 0; 최후 폴백)
//
//   - HTOF==2 페어 집계(특수 18..24 페어 + others), 인접(i-(i+1)) 히스토그램 포함

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

// ================================================================
namespace HBS {

// -------- geometry assumptions (E72-like) --------
static const int    kNBH2Seg   = 15;     // BH2 seg: 0..14
static const double kBH2_x0    = -10.0;  // [mm] world center
static const double kBH2_y0    =   0.0;
static const double kBH2_z0    = -560.0;
static const double kBH2_pitch = 14.0;   // [mm]
static const double kBH2_total = kNBH2Seg * kBH2_pitch; // 210 mm

static const int    kNHTOF     = 34;     // HTOF tiles: 0..33

// ROOT dictionary for vector<TParticle>
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
static int MapBH2_WorldToSeg(double x, double, double)
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

static inline double pct(double num, double den){ return (den>0)? (100.0*num/den):0.0; }

// -------- Pair counting for HTOF==2 (special 18..24 pairs + others) --------
struct PairCount {
  std::map<std::pair<int,int>, long long> counts;
  long long others = 0;
  PairCount(){ for(int i=18;i<=24;++i) counts[{i,i+1}] = 0; }

  void Fill(const std::set<int>& tiles){
    if((int)tiles.size()!=2) return;
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

// -------- Stats container --------
struct StatsBox {
  long long N_total = 0;
  long long N_BH2sel = 0;
  long long N_ge1 = 0;
  long long N_ge2 = 0;
  long long N_eq2 = 0;

  PairCount pairSummary;

  TH1I* hHTOF = nullptr;
  TH1I* hPairAdj = nullptr;

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

// -------- Section spec & parser --------
struct SectionSpec {
  bool isAny = false;
  int  lo = 0;
  int  hi = 0;
  TString tag;   // "BH2(any)" or "BH2 seg lo-hi"
  TString key;   // canvas/hist key
};

static inline TString TrimWS(const TString& s){
  TString t=s; t.ReplaceAll(" ",""); t.ReplaceAll("\t",""); return t;
}

// "4-9,4-11,3-8" → [(4,9),(4,11),(3,8)], plus "any"
static std::vector<SectionSpec> ParseRanges(const char* rangesCSV){
  std::vector<SectionSpec> out;
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

// -------- thresholds --------
static inline double MIPThresholdMeV(double mipFrac, double mipMeVperCm, double thickness_mm){
  return mipFrac * mipMeVperCm * (thickness_mm*0.1); // mm→cm
}

// -------- Extract ID(edep) from TParticle (embedded path) --------
// preferEmbedded=true: Mother(1)=copyNo, Weight=edep 우선 사용
static inline bool ExtractIdEdep(const TParticle& p,
                                 int nmax,
                                 int& id_out,
                                 double& edep_out,
                                 bool preferEmbedded=true)
{
  id_out = -1;
  edep_out = 0.0;

  if(preferEmbedded){
    int id = p.GetMother(1);
    double w = p.GetWeight();
    if(0 <= id && id < nmax && w > 0){
      id_out = id;
      edep_out = w;   // [MeV]
      return true;
    }
  }

  // 보조: StatusCode를 ID처럼, Weight를 edep처럼
  {
    int id = p.GetStatusCode();
    double w = p.GetWeight();
    if(0 <= id && id < nmax && w > 0){
      id_out = id;
      edep_out = w;
      return true;
    }
  }

  return false;
}

} // namespace HBS
// ================================================================


// ========================== 공개 함수 ===========================

void HTOF_BH2_summary_v2(const char* filename="E45.root",
                         const char* treename="g4hyptpc",
                         const char* ranges="4-9,4-11,3-8",
                         bool save=true,
                         double mipFrac=0.1,
                         double mipMeVperCm=2.0,
                         double BH2_thickness_mm=5.0,
                         double HTOF_thickness_mm=10.0)
{
  using namespace HBS;

  // 파일/트리
  TFile* f = TFile::Open(filename,"READ");
  if(!f || f->IsZombie()){ std::cerr<<"[ERR] open "<<filename<<" failed\n"; return; }
  TTree* T = (TTree*)f->Get(treename);
  if(!T){ std::cerr<<"[ERR] tree "<<treename<<" not found\n"; return; }

  // 브랜치 존재 확인
  if(!HasBranch(T,"BH2")){ std::cerr<<"[ERR] need branch 'BH2' (vector<TParticle>)\n"; return; }
  const bool hasHTOFvec  = HasBranch(T,"HTOF");
  const bool hasHTOFcopy = HasBranch(T,"HTOF_copyNo"); // 호환용(없어도 됨)
  const bool hasBH2edep  = HasBranch(T,"BH2_edep");
  const bool hasHTOFedep = HasBranch(T,"HTOF_edep");

  // 브랜치 포인터
  std::vector<TParticle>* BH2 = nullptr;
  std::vector<TParticle>* HTOF = nullptr;
  std::vector<int>* HTOF_copyNo = nullptr;
  std::vector<double>* BH2_edep = nullptr;
  std::vector<double>* HTOF_edep = nullptr;

  T->SetBranchAddress("BH2",&BH2);
  if(hasHTOFvec)  T->SetBranchAddress("HTOF",&HTOF);
  if(hasHTOFcopy) T->SetBranchAddress("HTOF_copyNo",&HTOF_copyNo);
  if(hasBH2edep)  T->SetBranchAddress("BH2_edep",&BH2_edep);
  if(hasHTOFedep) T->SetBranchAddress("HTOF_edep",&HTOF_edep);

  // 컷 값(MeV)
  const double thrBH2  = MIPThresholdMeV(mipFrac, mipMeVperCm, BH2_thickness_mm);   // 기본 0.10 MeV
  const double thrHTOF = MIPThresholdMeV(mipFrac, mipMeVperCm, HTOF_thickness_mm);  // 기본 0.20 MeV

  std::cout<<std::fixed<<std::setprecision(3);
  std::cout<<"[INFO] Energy cuts: mipFrac="<<mipFrac
           <<", dE/dx="<<mipMeVperCm<<" MeV/cm"
           <<", BH2_th="<<BH2_thickness_mm<<" mm (thr="<<thrBH2<<" MeV)"
           <<", HTOF_th="<<HTOF_thickness_mm<<" mm (thr="<<thrHTOF<<" MeV)\n";

  if(!(hasBH2edep || hasHTOFedep)){
    std::cout<<"[INFO] No *_edep branches. Will use embedded TParticle fields: "
                "Mother(1)=copyNo, Weight=edep (BVH 방식).\n";
  }

  // 섹션 스펙
  std::vector<SectionSpec> specs = ParseRanges(ranges);

  // 섹션 런 박스
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

    TString hname, htitle;
    hname.Form("hHTOF_%s", sp.key.Data());
    if(sp.isAny) htitle = "HTOF pattern | BH2(any seg);HTOF tile;Events";
    else         htitle.Form("HTOF pattern | %s;HTOF tile;Events", sp.tag.Data());
    r.box.InitHist(hname, htitle);

    TString hp, tp;
    hp.Form("hPairs_%s", sp.key.Data());
    if(sp.isAny) tp = "Adjacent HTOF pairs (i-(i+1)) | BH2(any seg);pair;Events";
    else         tp.Form("Adjacent HTOF pairs (i-(i+1)) | %s;pair;Events", sp.tag.Data());
    r.box.InitPairHist(hp, tp);

    runs.push_back(r);
  }

  // ===== 이벤트 루프 =====
  const Long64_t N = T->GetEntries();
  for(Long64_t ie=0; ie<N; ++ie){
    T->GetEntry(ie);
    for(auto& r : runs) r.box.N_total++;

    // -------- BH2: 세그별 edep 누적 --------
    std::map<int,double> bh2E;
    std::set<int> bh2_anySeg; // 폴백용
    if(BH2){
      for(size_t i=0; i<BH2->size(); ++i){
        const TParticle& p = BH2->at(i);
        int sid = -1;
        double edep = 0.0;
        bool ok=false;

        // (a) *_edep 브랜치가 있다면: edep=브랜치, ID=Mother(1) or StatusCode
        if(hasBH2edep && BH2_edep && i < (size_t)BH2_edep->size()){
          int idM = p.GetMother(1);
          int idS = p.GetStatusCode();
          int id  = (0<=idM && idM<kNBH2Seg) ? idM :
                    (0<=idS && idS<kNBH2Seg) ? idS : -1;
          if(id>=0){
            sid  = id;
            edep = BH2_edep->at(i);
            ok   = true;
          }
        }

        // (b) 임베디드 Mother(1)/Weight 사용 (BVH 규약)
        if(!ok){
          if(ExtractIdEdep(p, kNBH2Seg, sid, edep, /*preferEmbedded=*/true)){
            ok = true;
          }
        }

        // (c) 최후 폴백: 좌표→세그 (에너지는 모름→0)
        if(!ok){
          int s = MapBH2_WorldToSeg(p.Vx(), p.Vy(), p.Vz());
          if(0 <= s && s < kNBH2Seg){
            sid  = s;
            edep = 0.0;
            ok   = true;
          }
        }

        if(ok && 0 <= sid && sid < kNBH2Seg){
          bh2E[sid] += edep;
          bh2_anySeg.insert(sid);
        }
      }
    }

    std::set<int> bh2Seg_valid;
    bool anyEdepBH2 = false;
    for(const auto& kv : bh2E){ if(kv.second>0){ anyEdepBH2=true; break; } }

    if(anyEdepBH2){
      for(const auto& kv : bh2E){
        if(kv.second >= thrBH2) bh2Seg_valid.insert(kv.first);
      }
    }else{
      // edep가 전혀 누적되지 않은 경우(브랜치/Weight가 없거나 0) → 과거동작과 동일 폴백
      bh2Seg_valid = bh2_anySeg;
    }
    const bool passBH2_any = !bh2Seg_valid.empty();

    // -------- HTOF: 타일별 edep 누적 --------
    std::map<int,double> htofE;
    std::set<int> htof_anyId; // 폴백용

    // (1) copyNo/edep 브랜치가 있으면 사용 (호환)
    if(hasHTOFcopy && HTOF_copyNo){
      const size_t n = HTOF_copyNo->size();
      for(size_t i=0; i<n; ++i){
        int cn = HTOF_copyNo->at(i);
        if(0 <= cn && cn < kNHTOF){
          double edep = 0.0;
          if(hasHTOFedep && HTOF_edep && i<(size_t)HTOF_edep->size()) edep = HTOF_edep->at(i);
          htofE[cn] += edep;
          htof_anyId.insert(cn);
        }
      }
    }
    // (2) 벡터<TParticle>에서 Mother(1)/Weight 사용 (BVH 규약)
    else if(hasHTOFvec && HTOF){
      for(size_t i=0; i<HTOF->size(); ++i){
        const TParticle& p = HTOF->at(i);
        int tid = -1;
        double edep = 0.0;
        bool ok=false;

        // (a) *_edep가 있다면 edep=브랜치, ID=Mother(1)/StatusCode
        if(hasHTOFedep && HTOF_edep && i < (size_t)HTOF_edep->size()){
          int idM = p.GetMother(1);
          int idS = p.GetStatusCode();
          int id  = (0<=idM && idM<kNHTOF) ? idM :
                    (0<=idS && idS<kNHTOF) ? idS : -1;
          if(id>=0){
            tid  = id;
            edep = HTOF_edep->at(i);
            ok   = true;
          }
        }

        // (b) 임베디드 Mother(1)/Weight
        if(!ok){
          if(ExtractIdEdep(p, kNHTOF, tid, edep, /*preferEmbedded=*/true)){
            ok = true;
          }
        }

        // (c) 실패 → StatusCode만 ID로 사용, edep=0
        if(!ok){
          int idS = p.GetStatusCode();
          if(0<=idS && idS<kNHTOF){
            tid  = idS;
            edep = 0.0;
            ok   = true;
          }
        }

        if(ok && 0<=tid && tid<kNHTOF){
          htofE[tid] += edep;
          htof_anyId.insert(tid);
        }
      }
    }

    std::set<int> htofSeg_valid;
    bool anyEdepHTOF = false;
    for(const auto& kv : htofE){ if(kv.second>0){ anyEdepHTOF=true; break; } }

    if(anyEdepHTOF){
      for(const auto& kv : htofE){
        if(kv.second >= thrHTOF) htofSeg_valid.insert(kv.first);
      }
    }else{
      // edep가 전혀 없으면 과거처럼 "히트 존재"로 대체
      htofSeg_valid = htof_anyId;
    }
    const int htMult = (int)htofSeg_valid.size();

    // -------- 섹션별 집계 --------
    for(auto& r : runs){
      bool passBH2 = r.spec.isAny ? passBH2_any
                                  : BH2HitsInRange(bh2Seg_valid, r.spec.lo, r.spec.hi);
      if(passBH2){
        r.box.N_BH2sel++;
        if(htMult>=1) r.box.N_ge1++;
        if(htMult>=2) r.box.N_ge2++;
        if(htMult==2){
          r.box.N_eq2++;
          r.box.pairSummary.Fill(htofSeg_valid);
          r.box.FillPairHistIfAdjacent(htofSeg_valid);
        }
        for(int t : htofSeg_valid) r.box.hHTOF->Fill(t);
      }
    }
  }

  // ===== 출력 + 그림 =====
  for(size_t i=0;i<runs.size(); ++i){
    auto& r = runs[i];
    TString title = r.spec.tag.Length()? r.spec.tag : "BH2(any seg)";
    r.box.Print(title);

    TString cname1; cname1.Form("cHTOF_%s", r.spec.key.Data());
    r.cHTOF = new TCanvas(cname1, title, 800, 600);
    r.box.hHTOF->Draw("hist");

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

// Veto_sum.C
// ------------------------------------------------------------------
// 사용법 예시:
//
// root [0] .L gen.dict.C+
// root [1] .L Veto_sum.C+
//
// // 1) 인덱스 구축
// root [2] VetoSum_build("E45_BVH4.root", 0.10, 0.04, 0.04, true);
//
// // 2) 질문에서 주신 (BH2, U, D) 묶음 일괄 출력
// //    분모 = 각 BH2 세그먼트 母수 (BH2=h 기준)
// root [3] VetoSum_run_batch_default(false);
//
// //    분모 = BH2(any) 母수 (BH2에 hit한 전체 이벤트 母수)
// root [4] VetoSum_run_batch_default(true);
//
// // 3) CSV 저장 (선택)
// root [5] VetoSum_export_batch_default_csv("veto_batch.csv", false);
// ------------------------------------------------------------------

#include "Rtypes.h"
#include "TFile.h"
#include "TTree.h"
#include "TParticle.h"
#include <vector>
#include <unordered_set>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>

// -------- segmentation --------
static const int N_BH2  = 15; // 0..14
static const int N_BVHU = 22; // 0..21
static const int N_BVHD = 32; // 0..31

// -------- storage (global) --------
static bool g_vs_built = false;
static std::vector<Long64_t> g_vs_denBH2;  // per-BH2 segment denominator
static Long64_t g_vs_denAnyBH2 = 0;        // any-BH2 denominator
static std::vector<std::vector<std::vector<Long64_t>>> g_vs_cnt; // [BH2][U][D]

// -------- helper: unique segs per event above edep cut --------
static inline void vs_collect_unique_above_cut(const std::vector<TParticle>* v,
                                               int nmax, double cutMeV,
                                               std::vector<int>& out)
{
  out.clear();
  if(!v) return;
  std::unordered_set<int> s;
  s.reserve(8);
  for(const auto& p : *v){
    if(p.GetWeight() <= cutMeV) continue;      // energy deposition cut (MeV)
    int id = p.GetMother(1);
    if(0<=id && id<nmax) s.insert(id);
  }
  out.assign(s.begin(), s.end());
  std::sort(out.begin(), out.end());
}

// ================== Builder ==================
void VetoSum_build(const char* fname="E45_BVH4.root",
                   double ecutBH2MeV = 0.10,  // BH2(5mm) ~ 0.1 MIP
                   double ecutUMeV   = 0.04,  // U/D(2mm) ~ 0.1 MIP
                   double ecutDMeV   = 0.04,
                   bool   dedupPerEvent = true)
{
  g_vs_built = false;
  g_vs_denBH2.assign(N_BH2, 0);
  g_vs_denAnyBH2 = 0;
  g_vs_cnt.assign(N_BH2,
    std::vector<std::vector<Long64_t>>(N_BVHU, std::vector<Long64_t>(N_BVHD, 0)));

  TFile* f = TFile::Open(fname, "READ");
  if(!f || f->IsZombie()){ std::cerr << "[ERR] cannot open " << fname << "\n"; return; }
  TTree* tr = (TTree*)f->Get("g4hyptpc");
  if(!tr){ std::cerr << "[ERR] no TTree g4hyptpc\n"; return; }

  std::vector<TParticle> *BH2=nullptr, *BVHU=nullptr, *BVHD=nullptr;
  tr->SetBranchAddress("BH2",   &BH2);
  tr->SetBranchAddress("BVH_U", &BVHU);
  tr->SetBranchAddress("BVH_D", &BVHD);

  const Long64_t N = tr->GetEntries();
  std::vector<int> H,U,D; H.reserve(8); U.reserve(8); D.reserve(8);

  for(Long64_t i=0;i<N;++i){
    tr->GetEntry(i);

    vs_collect_unique_above_cut(BH2,  N_BH2,  ecutBH2MeV, H);
    vs_collect_unique_above_cut(BVHU, N_BVHU, ecutUMeV,   U);
    vs_collect_unique_above_cut(BVHD, N_BVHD, ecutDMeV,   D);

    if(H.empty()) continue; // BH2 미히트 이벤트 제외(母수 정의에 맞춤)

    // any-BH2 母수(이벤트 단위)
    ++g_vs_denAnyBH2;

    // per-BH2 母수(세그먼트 단위)
    if(dedupPerEvent){
      for(int h : H) g_vs_denBH2[h] += 1;
    }else{
      for(int h : H) g_vs_denBH2[h] += 1; // per-hit로 바꾸려면 hit 루프 필요
    }

    // numerator: (BH2=h) AND (U=u) AND (D=d)
    for(int h: H){
      for(int u: U){
        for(int d: D){
          g_vs_cnt[h][u][d] += 1;
        }
      }
    }
  }

  g_vs_built = true;

  std::cout << "\n==================== VetoSum built ====================\n";
  std::cout << "File                : " << fname << "\n";
  std::cout << "Total events        : " << N << "\n";
  std::cout << "Any-BH2 denominator : " << g_vs_denAnyBH2 << " events\n";
  std::cout << "Per-BH2 events (denominator):\n";
  for(int h=0; h<N_BH2; ++h){
    std::cout << "  BH2 " << std::setw(2) << h << " : " << g_vs_denBH2[h] << "\n";
  }
  std::cout << "Cuts (MeV)          : BH2="<<ecutBH2MeV
            << ", U="<<ecutUMeV << ", D="<<ecutDMeV << "\n";
  std::cout << "=======================================================\n";
}

// ================== Ratio / Printing ==================
struct VS_SegTriple { int h, u, d; }; // (BH2, BVH_U, BVH_D)

// denomAllBH2=false : 분모 = BH2=h 母수 (세그먼트별)
// denomAllBH2=true  : 분모 = BH2(any) 母수 (전체)
double VetoSum_ratio(int bh2, int u, int d, bool denomAllBH2=false)
{
  if(!g_vs_built){ std::cerr << "[ERR] VetoSum_build() first.\n"; return 0.0; }
  if(bh2<0 || bh2>=N_BH2 || u<0 || u>=N_BVHU || d<0 || d>=N_BVHD){
    std::cerr << "[ERR] index out of range.\n"; return 0.0;
  }
  const double num = (double) g_vs_cnt[bh2][u][d];
  const double den = denomAllBH2
                   ? (double) std::max<Long64_t>(1, g_vs_denAnyBH2)
                   : (double) std::max<Long64_t>(1, g_vs_denBH2[bh2]);
  return 100.0 * num / den;
}

void VetoSum_print_one(int bh2, int u, int d, bool denomAllBH2=false)
{
  if(!g_vs_built){ std::cerr << "[ERR] VetoSum_build() first.\n"; return; }
  if(bh2<0 || bh2>=N_BH2 || u<0 || u>=N_BVHU || d<0 || d>=N_BVHD){
    std::cerr << "[ERR] index out of range.\n"; return;
  }
  Long64_t num = g_vs_cnt[bh2][u][d];
  Long64_t den = denomAllBH2 ? g_vs_denAnyBH2 : g_vs_denBH2[bh2];
  double r = VetoSum_ratio(bh2,u,d,denomAllBH2);

  std::cout << "\n--- VETO ratio ---\n"
            << " BH2=" << bh2 << ", U=" << u << ", D=" << d << "\n"
            << " Numerator  : " << num << " events\n"
            << " Denominator: " << den << " events ("
            << (denomAllBH2 ? "any-BH2" : ("BH2="+std::to_string(bh2))) << ")\n"
            << " Ratio      : " << std::fixed << std::setprecision(3) << r << " %\n";
}

void VetoSum_print_batch(const std::vector<VS_SegTriple>& triples, bool denomAllBH2=false)
{
  if(!g_vs_built){ std::cerr << "[ERR] VetoSum_build() first.\n"; return; }
  std::cout << "\n== VETO ratios (denominator: " << (denomAllBH2?"any-BH2":"BH2=h") << ") ==\n";
  std::cout << "   BH2   U    D     Num        Den        Ratio(%)\n";
  std::cout << "---------------------------------------------------\n";

  for(const auto& t : triples){
    if(t.h<0 || t.h>=N_BH2 || t.u<0 || t.u>=N_BVHU || t.d<0 || t.d>=N_BVHD){
      std::cout << " [SKIP] out of range: ("<<t.h<<","<<t.u<<","<<t.d<<")\n";
      continue;
    }
    Long64_t num = g_vs_cnt[t.h][t.u][t.d];
    Long64_t den = denomAllBH2 ? g_vs_denAnyBH2 : g_vs_denBH2[t.h];
    double r = den ? (100.0 * (double)num / (double)den) : 0.0;

    std::cout.setf(std::ios::right);
    std::cout << std::setw(5) << t.h
              << std::setw(5) << t.u
              << std::setw(5) << t.d
              << std::setw(10) << num
              << std::setw(12) << den
              << std::setw(12) << std::fixed << std::setprecision(3) << r
              << "\n";
  }
}

// ================== Default list (질문 목록 그대로) ==================
static std::vector<VS_SegTriple> VetoSum_default_list()
{
  return {
    {7, 9,16},{7, 9,17},{7, 9,18},{7, 9,19},{7, 9,20},{7, 9,21},
    {7,10,21},{7,10,20},{7,10,19},{7,10,18},{7,10,17},
    {8,10,18},{8,10,19},{8,10,20},{8,10,21},{8,10,22},
    {8,11,22},{8,11,21},{8,11,20},{8,11,19},{8,11,18},{8,11,17},
    {8,12,20},{8,12,19},{8,12,18},
    {6, 7,16},{6, 7,17},{6, 7,18},{6, 7,19},{6, 7,20},
    {6, 8,21},{6, 8,20},{6, 8,19},{6, 8,18},{6, 8,17},{6, 8,16},{6, 8,15},
    {6, 9,16},{6, 9,17},{6, 9,18},{6, 9,19},{5, 7,14}, {5, 7,15}, {5, 7,16}, {5, 7,17}, {5, 7,18}, {5, 7,19}
    , {5, 7,20}, {5, 6,19}, {5, 6,18}, {5, 6,17}, {5, 6,16}, {5, 6,15}, {5, 6,14}, {4, 4,15}, {4, 4,16}
    , {4, 4,17}, {4, 5,13}, {4, 5,14}, {4, 5,15}, {4, 5,16}, {4, 5,17}, {4, 5,18}, {4, 5,19}, {4, 6,13}
    , {4, 6,14}, {4, 6,15}, {4, 6,16}, {4, 6,17}, {4, 6,18}, {3, 3,13}, {3, 3,14}, {3, 3,15}, {3, 3,16}, {3, 3,17}
    , {3, 4,12}, {3, 4,13}, {3, 4,14}, {3, 4,15}, {3, 4,16}, {3, 4,17}, {3, 4,18}, {9, 12, 18}, {9, 12, 19}, {9, 12, 20}
    , {9, 12, 21}, {9, 12, 22}, {9, 12, 23}, {9, 12, 24}, {9, 13, 19}, {9, 13, 20}, {9, 13, 21}, {9, 13, 22}, {9, 13, 23}
    , {10, 13, 21}, {10, 13, 22}, {10, 13, 23}, {10, 13, 24}, {10, 14, 20}, {10, 14, 21}, {10, 14, 22}, {10, 14, 23}, {10, 14, 24}
    , {11, 15, 21}, {11, 15, 22}, {11, 15, 23}, {11, 15, 24}, {11, 16, 21}, {11, 16, 22}, {11, 16, 23}, {11, 16, 24}, {12, 16, 22}
    , {12, 16, 23}, {12, 17, 22}, {12, 17, 23}, {12, 17, 24}, {13, 18, 22}, {13, 18, 23}, {13, 18, 24}, {13, 18, 25}, {13, 18, 26}
    , {13, 19, 25}, {13, 19, 24}, {13, 19, 24}, {13, 19, 22}, {2, 2, 12},{2, 2, 13}, {2, 2, 14}, {2, 2, 15}
    , {2, 3, 15}, {2, 3, 14}, {2, 3, 13}
  };
}

void VetoSum_run_batch_default(bool denomAllBH2=false)
{
  VetoSum_print_batch(VetoSum_default_list(), denomAllBH2);
}

// CSV 저장(선택)
void VetoSum_export_batch_csv(const std::vector<VS_SegTriple>& triples,
                              const char* out="veto_batch.csv",
                              bool denomAllBH2=false)
{
  if(!g_vs_built){ std::cerr << "[ERR] VetoSum_build() first.\n"; return; }
  std::ofstream ofs(out);
  if(!ofs){ std::cerr << "[ERR] cannot open " << out << " for write\n"; return; }
  ofs << "BH2,U,D,Num,Den,Ratio(%),DenMode\n";
  for(const auto& t : triples){
    if(t.h<0 || t.h>=N_BH2 || t.u<0 || t.u>=N_BVHU || t.d<0 || t.d>=N_BVHD) continue;
    Long64_t num = g_vs_cnt[t.h][t.u][t.d];
    Long64_t den = denomAllBH2 ? g_vs_denAnyBH2 : g_vs_denBH2[t.h];
    double r = den ? (100.0 * (double)num / (double)den) : 0.0;
    ofs << t.h << "," << t.u << "," << t.d << ","
        << num << "," << den << ","
        << std::fixed << std::setprecision(3) << r << ","
        << (denomAllBH2?"any-BH2":"BH2=h") << "\n";
  }
  ofs.close();
  std::cout << "[OK] CSV written: " << out << "\n";
}

// 기본 목록을 바로 CSV로 저장하는 편의 함수
void VetoSum_export_batch_default_csv(const char* out="veto_batch.csv",
                                      bool denomAllBH2=false)
{
  VetoSum_export_batch_csv(VetoSum_default_list(), out, denomAllBH2);
}

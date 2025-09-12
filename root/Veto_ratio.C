// Veto_ratio.C  — BH2 전체 母수/세그먼트 母수 선택 가능
// ------------------------------------------------------------------
// 예시:
//   .L Veto_ratio.C+
//
//   // 1) 바로 프롬프트 모드(전체 BH2 母수 사용)
//   prompt_veto_ratio("E45_BVH4.root", 0.10, 0.04, 0.04, true, /*denomAllBH2=*/true);
//
//   // 2) 인덱스 구축 후 수동 질의
//   build_veto_index("E45_BVH4.root", 0.10, 0.04, 0.04, true);
//   print_veto_ratio(6,15,20, /*denomAllBH2=*/true); // 전체 BH2 母수
//   print_veto_ratio(6,15,20);                      // BH2=6 母수 (기존 방식)
//   prompt_veto_ratio_loop(/*denomAllBH2=*/true);   // U H D 반복 질의
// ------------------------------------------------------------------

#include "Rtypes.h"
#include "TFile.h"
#include "TTree.h"
#include "TParticle.h"
#include <vector>
#include <unordered_set>
#include <tuple>
#include <algorithm>
#include <iostream>
#include <iomanip>

// segmentation
static const int N_BH2  = 15;
static const int N_BVHU = 22;
static const int N_BVHD = 32;

// storage
static bool g_built = false;
static std::vector<Long64_t> g_denBH2;   // per-BH2 segment denominator
static Long64_t g_denBH2_any = 0;        // any-BH2 denominator (이벤트 단위)
static std::vector<std::vector<std::vector<Long64_t>>> g_cnt;  // [BH2][U][D]

// helper: 이벤트당 중복 제거 + 에너지컷
static inline void collect_unique_above_cut(const std::vector<TParticle>* v,
                                            int nmax, double cutMeV,
                                            std::vector<int>& out)
{
  out.clear();
  if(!v) return;
  std::unordered_set<int> s;
  s.reserve(8);
  for(const auto& p : *v){
    if(p.GetWeight() <= cutMeV) continue; // Edep cut (MeV)
    int id = p.GetMother(1);
    if(0<=id && id<nmax) s.insert(id);
  }
  out.assign(s.begin(), s.end());
  std::sort(out.begin(), out.end());
}

// -------------------- 메인 빌더 --------------------
void Veto_ratio(const char* fname="E45_BVH4.root",
                double ecutBH2MeV = 0.10,   // BH2 ~5 mm → ~0.1 MIP
                double ecutUMeV   = 0.04,   // U/D ~2 mm → ~0.1 MIP
                double ecutDMeV   = 0.04,
                bool   dedupPerEvent = true)
{
  g_built = false;
  g_denBH2.assign(N_BH2, 0);
  g_denBH2_any = 0;
  g_cnt.assign(N_BH2, std::vector<std::vector<Long64_t>>(N_BVHU, std::vector<Long64_t>(N_BVHD, 0)));

  TFile* f = TFile::Open(fname, "READ");
  if(!f || f->IsZombie()){ std::cerr << "[ERR] cannot open " << fname << "\n"; return; }
  TTree* tr = (TTree*)f->Get("g4hyptpc");
  if(!tr){ std::cerr << "[ERR] no TTree g4hyptpc\n"; return; }

  std::vector<TParticle> *BH2=nullptr, *BVHU=nullptr, *BVHD=nullptr;
  tr->SetBranchAddress("BH2",   &BH2);
  tr->SetBranchAddress("BVH_U", &BVHU);
  tr->SetBranchAddress("BVH_D", &BVHD);

  const Long64_t N = tr->GetEntries();
  std::vector<int> H,U,D;
  H.reserve(8); U.reserve(8); D.reserve(8);

  for(Long64_t i=0;i<N;++i){
    tr->GetEntry(i);

    collect_unique_above_cut(BH2,  N_BH2,  ecutBH2MeV, H);
    collect_unique_above_cut(BVHU, N_BVHU, ecutUMeV,   U);
    collect_unique_above_cut(BVHD, N_BVHD, ecutDMeV,   D);

    if(H.empty()) continue;         // BH2가 안 맞은 이벤트는 제외(母수 정의에 맞게)

    // any-BH2 母수(이벤트 단위)
    ++g_denBH2_any;

    // per-BH2 母수(세그먼트 조건)
    if(dedupPerEvent){
      for(int h : H) g_denBH2[h] += 1;
    }else{
      for(int h : H) g_denBH2[h] += 1; // per-hit 로 바꾸려면 원시 hit 루프 필요
    }

    // numerator: (BH2=h) AND (U=u) AND (D=d)
    for(int h: H){
      for(int u: U){
        for(int d: D){
          g_cnt[h][u][d] += 1;
        }
      }
    }
  }

  g_built = true;

  std::cout << "\n==================== Index built ====================\n";
  std::cout << "File                : " << fname << "\n";
  std::cout << "Total events        : " << N << "\n";
  std::cout << "Any-BH2 denominator : " << g_denBH2_any << " events\n";
  std::cout << "Per-BH2 events (denominator):\n";
  for(int h=0; h<N_BH2; ++h){
    std::cout << "  BH2 " << std::setw(2) << h << " : " << g_denBH2[h] << "\n";
  }
  std::cout << "Cuts (MeV)          : BH2="<<ecutBH2MeV
            << ", U="<<ecutUMeV << ", D="<<ecutDMeV << "\n";
  std::cout << "=====================================================\n";
}

void build_veto_index(const char* fname="E45_BVH4.root",
                      double ecutBH2MeV = 0.10,
                      double ecutUMeV   = 0.04,
                      double ecutDMeV   = 0.04,
                      bool   dedupPerEvent = true)
{
  Veto_ratio(fname, ecutBH2MeV, ecutUMeV, ecutDMeV, dedupPerEvent);
}

// denomAllBH2=false : BH2=h 母수 (기존 방식)
// denomAllBH2=true  : BH2(any) 母수 (요청하신 방식)
double veto_ratio(int bh2, int u, int d, bool denomAllBH2 /*=true*/)
{
  if(!g_built){ std::cerr << "[ERR] build_veto_index() first.\n"; return 0.0; }
  if(bh2<0 || bh2>=N_BH2 || u<0 || u>=N_BVHU || d<0 || d>=N_BVHD){
    std::cerr << "[ERR] index out of range.\n"; return 0.0;
  }
  const double num = (double) g_cnt[bh2][u][d];
  const double den = denomAllBH2
                   ? (double) std::max<Long64_t>(1, g_denBH2_any)
                   : (double) std::max<Long64_t>(1, g_denBH2[bh2]);
  return 100.0 * num / den;
}

void print_veto_ratio(int bh2, int u, int d, bool denomAllBH2 /*=false*/)
{
  if(!g_built){ std::cerr << "[ERR] build_veto_index() first.\n"; return; }
  if(bh2<0 || bh2>=N_BH2 || u<0 || u>=N_BVHU || d<0 || d>=N_BVHD){
    std::cerr << "[ERR] index out of range.\n"; return;
  }
  Long64_t num = g_cnt[bh2][u][d];
  Long64_t den = denomAllBH2 ? g_denBH2_any : g_denBH2[bh2];
  double r = veto_ratio(bh2,u,d,denomAllBH2);
  std::cout << "\n--- VETO ratio ---\n"
            << " BH2=" << bh2 << ", U=" << u << ", D=" << d << "\n"
            << " Numerator  : " << num << " events\n"
            << " Denominator: " << den << " events ("
            << (denomAllBH2 ? "any-BH2" : ("BH2="+std::to_string(bh2))) << ")\n"
            << " Ratio      : " << std::fixed << std::setprecision(3) << r << " %\n";
}

// BH2 하나에서 상위 U-D 페어 나열 (표시는 항상 BH2=h 母수 기준)
void list_top_pairs(int bh2, int topN=10, bool byRatio=false)
{
  if(!g_built){ std::cerr << "[ERR] build_veto_index() first.\n"; return; }
  if(bh2<0 || bh2>=N_BH2){ std::cerr << "[ERR] BH2 out of range.\n"; return; }

  struct Item{ int u,d; Long64_t n; double r; };
  std::vector<Item> v;
  v.reserve(N_BVHU*N_BVHD);
  const double den = (double) std::max<Long64_t>(1, g_denBH2[bh2]);
  for(int u=0; u<N_BVHU; ++u){
    for(int d=0; d<N_BVHD; ++d){
      Long64_t n = g_cnt[bh2][u][d];
      if(n>0){
        v.push_back({u,d,n, 100.0 * (double)n / den});
      }
    }
  }
  if(byRatio){
    std::sort(v.begin(), v.end(), [](const Item& a, const Item& b){
      if(a.r==b.r) return a.n>b.n;
      return a.r>b.r;
    });
  }else{
    std::sort(v.begin(), v.end(), [](const Item& a, const Item& b){
      if(a.n==b.n) return a.r>b.r;
      return a.n>b.n;
    });
  }

  const int M = std::min<int>(topN, (int)v.size());
  std::cout << "\nTop " << M << " pairs in BH2=" << bh2
            << " by " << (byRatio?"ratio":"count") << "\n";
  std::cout << "  rank : (U,D)  count   ratio(%)\n";
  for(int i=0;i<M;++i){
    std::cout << "  " << std::setw(4) << (i+1) << " : ("
              << std::setw(2) << v[i].u << ","
              << std::setw(2) << v[i].d << ")  "
              << std::setw(8) << v[i].n << "  "
              << std::fixed << std::setprecision(3) << std::setw(8) << v[i].r << "\n";
  }
}

// -------------------- 인터랙티브 모드 --------------------
void prompt_veto_ratio_loop(bool denomAllBH2 /*=false*/)
{
  if(!g_built){
    std::cerr << "[ERR] build_veto_index() first (or call prompt_veto_ratio(...)).\n";
    return;
  }
  std::cout << "\nEnter segments as:  U  H  D   (negative to quit)"
            << "\nDenominator mode : " << (denomAllBH2 ? "any-BH2" : "BH2=h") << "\n";
  while(true){
    std::cout << "U H D ? ";
    int u,h,d;
    if(!(std::cin >> u)) break;
    if(u < 0) break;
    if(!(std::cin >> h >> d)) break;
    if(h < 0 || d < 0) break;

    if(u>=0 && u<N_BVHU && h>=0 && h<N_BH2 && d>=0 && d<N_BVHD){
      print_veto_ratio(h,u,d,denomAllBH2);
    }else{
      std::cout << "[WARN] out of range. Valid: U[0.." << (N_BVHU-1)
                << "] H[0.." << (N_BH2-1) << "] D[0.." << (N_BVHD-1) << "]\n";
    }
  }
  std::cout << "Bye.\n";
}

void prompt_veto_ratio(const char* fname="E45_BVH4.root",
                       double ecutBH2MeV = 0.10,
                       double ecutUMeV   = 0.04,
                       double ecutDMeV   = 0.04,
                       bool   dedupPerEvent = true,
                       bool   denomAllBH2 = false)
{
  build_veto_index(fname, ecutBH2MeV, ecutUMeV, ecutDMeV, dedupPerEvent);
  prompt_veto_ratio_loop(denomAllBH2);
}

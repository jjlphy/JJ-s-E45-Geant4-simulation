// ComputeOverkillFromOverlay.C
// (Use overlay triplets {BH2, BVH_U, BVH_D} as veto mask to compute Overkill rate on 2pi sample)
//
// 사용법 예:
//   root -l
//   .L ComputeOverkillFromOverlay.C+
//   // 기본값(파일 경로/컷/세그먼트 수)은 아래 run_overkill()에 정의됨
//   run_overkill();
//
// 입력 파일 형식 (혼합 허용):
//   {3, 4, 17}
//   3 4 17
//   3,4,17
//
// 정의:
//   - 분모(denominator): 이벤트에 BH2 && BVH_U 히트가 존재 (에너지컷 적용)
//   - 분자(numerator): 위 분모 조건을 만족하면서, Overlay Triplets 중 하나와
//                     (BH2, BVH_U, SCH) 조합이 일치 → Veto에 의해 잘림
//   - Overkill rate = numerator / denominator
//
// 주의:
//   - EX(제외) 트리플릿은 입력 리스트에 이미 빠졌다고 가정(시각화와 동일).
//     안전을 위해 코드에서 BH2/BVHU/SCH의 유효 범위 검사는 수행함.

#include "Rtypes.h"
#include "TFile.h"
#include "TTree.h"
#include "TParticle.h"
#include "TSystem.h"
#include "TInterpreter.h"
#include "TROOT.h"

#include <vector>
#include <tuple>
#include <unordered_set>
#include <unordered_map>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <cctype>

// -------- Detector Segmentation (BVH_3D_ver4와 동일하게 설정할 것) --------
static const int N_BH2  = 15; // 0..14
static const int N_BVHU = 22; // 0..14
static const int N_SCH  = 64; // 0..63

// ===== Helper: branch auto-bind =====
static bool bind_branch(TTree* tr, const char* preferred,
                        const std::vector<const char*>& fallbacks,
                        std::vector<TParticle>*& ptr) {
  if (tr->GetBranch(preferred)) { tr->SetBranchAddress(preferred, &ptr); return true; }
  for (auto nm : fallbacks) {
    if (tr->GetBranch(nm)) { tr->SetBranchAddress(nm, &ptr);
      Warning("Overkill","using fallback branch '%s' for '%s'", nm, preferred);
      return true;
    }
  }
  Error("Overkill","branch '%s' not found (no fallback matched)", preferred);
  return false;
}

// ===== Helper: unique segment indices passing e-dep cut =====
static inline void get_unique_hits(const std::vector<TParticle>* v,
                                   int nmax, double cutMeV,
                                   std::vector<int>& out)
{
  out.clear();
  if(!v) return;
  std::unordered_set<int> s;
  s.reserve(8);
  for(const auto& p : *v){
    const double edepMeV = p.GetWeight();
    if(edepMeV < cutMeV) continue;
    int seg = p.GetMother(1);
    if(0<=seg && seg<nmax) s.insert(seg);
  }
  out.assign(s.begin(), s.end());
  std::sort(out.begin(), out.end());
}

// ===== Helper: (BH2,BVHU,SCH) → 단일 int 키 =====
static inline int key3(int h,int u,int s){ return (h<<16) | (u<<8) | s; }

// ===== Parser: overlay_triplets.txt 읽어서 (h,u,s) 집합 구성 =====
static bool load_overlay_triplets(const char* path,
                                  std::unordered_set<int>& overlay_set,
                                  std::vector<std::tuple<int,int,int>>& overlay_list)
{
  overlay_set.clear();
  overlay_list.clear();

  std::ifstream fin(path);
  if(!fin.is_open()){
    Error("Overkill","Cannot open overlay list file: %s", path);
    return false;
  }

  auto trim = [](std::string& s){
    size_t a = s.find_first_not_of(" \t\r\n");
    size_t b = s.find_last_not_of(" \t\r\n");
    if(a==std::string::npos) { s.clear(); return; }
    s = s.substr(a, b-a+1);
  };

  std::string line;
  int ln = 0;
  while(std::getline(fin, line)){
    ln++;
    trim(line);
    if(line.empty()) continue;

    // 괄호/중괄호/콤마 제거를 위해 구분자 통일
    for(char& c : line){
      if(c=='{' || c=='}' || c=='[' || c==']' || c=='(' || c==')' || c==',') c=' ';
    }
    std::istringstream iss(line);
    int h,u,s;
    if(!(iss >> h >> u >> s)){
      Warning("Overkill","Cannot parse line %d: '%s' (skipped)", ln, line.c_str());
      continue;
    }
    if(!(0<=h && h<N_BH2 && 0<=u && u<N_BVHU && 0<=s && s<N_SCH)){
      Warning("Overkill","Out-of-range triplet at line %d: {%d,%d,%d} (skipped)", ln, h,u,s);
      continue;
    }
    int k = key3(h,u,s);
    if(overlay_set.insert(k).second){
      overlay_list.emplace_back(h,u,s);
    }
  }
  fin.close();

  // 보기 좋게 정렬
  std::sort(overlay_list.begin(), overlay_list.end(),
            [](const std::tuple<int,int,int>& a, const std::tuple<int,int,int>& b){
              if (std::get<0>(a) != std::get<0>(b)) return std::get<0>(a) < std::get<0>(b);
              if (std::get<1>(a) != std::get<1>(b)) return std::get<1>(a) < std::get<1>(b);
              return std::get<2>(a) < std::get<2>(b);
            });

  Info("Overkill","Loaded %zu overlay triplets from %s", overlay_list.size(), path);
  return true;
}

// ================== Main Overkill Computation ==================
void ComputeOverkillFromOverlay(const char* data_root = "../E45_2pi_Ver3.root",
                                const char* overlay_txt = "overlay_triplets_ver3.txt",
                                double ecutBH2MeV   = 0.10,
                                double ecutUMeV     = 0.04,
                                double ecutSCHMeV   = 0.04)
{
  gSystem->Load("libPhysics");
  gInterpreter->GenerateDictionary("vector<TParticle>", "TParticle.h;vector");

  // 1) 오버레이 트리플릿 로드 → veto mask 구성
  std::unordered_set<int> overlay_set;                // 빠른 조회용
  std::vector<std::tuple<int,int,int>> overlay_list; // 보고용
  if(!load_overlay_triplets(overlay_txt, overlay_set, overlay_list)){
    Error("Overkill","Failed to load overlay triplets; abort.");
    return;
  }

  // 2) 데이터 열기
  TFile* f = TFile::Open(data_root, "READ");
  if(!f || f->IsZombie()){ Error("Overkill","Cannot open file: %s", data_root); return; }
  TTree* tr = (TTree*)f->Get("g4hyptpc");
  if(!tr){ Error("Overkill","Cannot find TTree 'g4hyptpc'"); f->Close(); return; }

  std::vector<TParticle> *BH2=nullptr, *BVHU=nullptr, *SCHv=nullptr;
  bool ok = true;
  ok &= bind_branch(tr, "BH2",   {},                           BH2);
  ok &= bind_branch(tr, "BVH_U", {"BVH"},                      BVHU);
  ok &= bind_branch(tr, "SCH",   {"/SCH","SCH1","SCHHits"},    SCHv);
  if(!ok){ f->Close(); return; }

  // 3) 카운팅 (분모/분자)
  Long64_t denom_all = 0, num_all = 0;
  std::vector<Long64_t> denom_h(N_BH2, 0), num_h(N_BH2, 0);

  const Long64_t N = tr->GetEntries();
  std::vector<int> hitsH, hitsU, hitsS;
  hitsH.reserve(8); hitsU.reserve(8); hitsS.reserve(8);

  std::cout << "[Info] Computing overkill on " << N << " events..." << std::endl;
  for(Long64_t i=0;i<N;++i){
    if(i && (i%100000==0)) std::cout << "  " << i << " / " << N << "\r" << std::flush;
    tr->GetEntry(i);

    get_unique_hits(BH2,  N_BH2,  ecutBH2MeV, hitsH);
    get_unique_hits(BVHU, N_BVHU, ecutUMeV,   hitsU);
    get_unique_hits(SCHv, N_SCH,  ecutSCHMeV, hitsS);

    if(hitsH.empty() || hitsU.empty()) continue; // 분모 조건: BH2 && BVH_U
    denom_all++;

    bool global_veto = false;
    for(int h : hitsH){
      bool veto_h = false;
      denom_h[h]++;

      if(!hitsS.empty()){
        for(int u : hitsU){
          for(int s : hitsS){
            const int k = key3(h,u,s);
            if(overlay_set.count(k)){ veto_h = true; global_veto = true; break; }
          }
          if(veto_h) break;
        }
      }

      if(veto_h) num_h[h]++;
    }
    if(global_veto) num_all++;
  }
  std::cout << "\n[Info] Done." << std::endl;

  // 4) 결과 출력
  std::cout << "\n========== Overkill Summary ==========\n";
  std::cout << std::fixed << std::setprecision(3);
  std::cout << "Data file             : " << data_root << "\n";
  std::cout << "Overlay list          : " << overlay_txt << " (" << overlay_list.size() << " cells)\n";
  std::cout << "Energy cuts (MeV)     : BH2=" << ecutBH2MeV
            << ", BVH_U=" << ecutUMeV << ", SCH=" << ecutSCHMeV << "\n";
  std::cout << "Denominator (BH2&&U)  : " << denom_all << "\n";
  if(denom_all>0){
    double R = (double)num_all / (double)denom_all * 100.0;
    std::cout << "Overkill rate (global): " << num_all << " / " << denom_all
              << "  =  " << R << " %\n";
  }else{
    std::cout << "Overkill rate (global): n/a (denominator=0)\n";
  }
  std::cout << "\nPer-BH2 segment overkill:\n";
  std::cout << "  h :  numerator / denominator  =  rate(%)\n";
  for(int h=0; h<N_BH2; ++h){
    if(denom_h[h]>0){
      double r = (double)num_h[h] / (double)denom_h[h] * 100.0;
      std::cout << " " << std::setw(2) << h << " :  "
                << std::setw(9) << num_h[h] << " / " << std::setw(9) << denom_h[h]
                << "  =  " << std::setw(7) << r << "\n";
    }else{
      std::cout << " " << std::setw(2) << h << " :  "
                << std::setw(9) << 0 << " / " << std::setw(9) << 0
                << "  =      n/a\n";
    }
  }
  std::cout << "=====================================\n";

  tr->ResetBranchAddresses();
  f->Close();
}

// ================== 편의 실행 래퍼 ==================
void run_overkill() {
  // BVH_3D_ver4.C와 동일한 2π 데이터 파일, 컷 값을 기본으로 둠
  const char* data_root   = "../E45_2pi_Ver3.root";
  const char* overlay_txt = "overlay_triplets_ver3.txt"; // BVH_3D_ver4 출력 복사해서 저장

  double ecutBH2 = 0.10;
  double ecutU   = 0.04;
  double ecutSCH = 0.04;

  ComputeOverkillFromOverlay(data_root, overlay_txt, ecutBH2, ecutU, ecutSCH);
}

// -*- C++ -*-
// E45_Summary_FromScreens.C
//   - 스샷에서 읽은 숫자(Beam, Trig2 / BeamVeto)를 넣으면
//     Beam(π−, π+) 효율·BG@1M(+95%CI), Reaction(π−π0p, π−π+n, π+π0p, π+π+n) Overkill(+95%CI) 요약 표 출력
//   - 콘솔 출력에도 BEAM의 95% CI를 함께 표시하도록 수정
//
// 사용법:
//   root -l
//   .L E45_Summary_FromScreens.C+
//   E45_Summary_FromScreens();     // 콘솔 표 출력 + CSV 저장
//
// CSV:
//   - E45_beam_summary.csv
//   - E45_reaction_summary.csv

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <map>
#include <fstream>
#include <algorithm>
#include <sstream>
#include <cmath>

#include "TEfficiency.h"  // Clopper–Pearson CI

// ---------------- CI helper (95%) ----------------
struct CI { double lo, hi; };  // fraction [0..1]

static inline CI CP95(long long k_success, long long n_total){
  // Binomial Clopper–Pearson 95% CI using ROOT's TEfficiency
  const double alpha = 0.05;
  CI c;
  c.lo = TEfficiency::ClopperPearson((int)n_total, (int)k_success, alpha, /*upper=*/false);
  c.hi = TEfficiency::ClopperPearson((int)n_total, (int)k_success, alpha, /*upper=*/true);
  return c;
}

// ---------------- 데이터 모델 ----------------
enum class Kind { kBeam, kReaction };
enum class Mask { tight, fit, wide, ultra };

static inline const char* MaskName(Mask m){
  switch(m){
    case Mask::tight: return "tight";
    case Mask::fit:   return "fit";
    case Mask::wide:  return "wide";
    case Mask::ultra: return "ultra";
  }
  return "?";
}

// 파일명 prefix → 표기 채널명
static inline std::string PrettyChannel(const std::string& prefix){
  if(prefix=="beam")             return "π⁻ beam";
  if(prefix=="beamplus")         return "π⁺ beam";
  if(prefix=="pi0p")             return "π⁻p → π⁻π⁰p";
  if(prefix=="piplusn")          return "π⁻p → π⁻π⁺n";
  if(prefix=="piplus_pi0p")      return "π⁺p → π⁺π⁰p";
  if(prefix=="piplusplus")       return "π⁺p → π⁺π⁺n";
  return prefix;
}

struct Row {
  Kind kind;                 // Beam or Reaction
  std::string prefix;        // beam / beamplus / pi0p / piplusn / piplus_pi0p / piplusplus
  double pGeV;               // 0.98, 1.02, ...
  long long N_beam;          // Beam(BH2 in-range)
  long long num[4];          // Beam: Trig2_* (tight..ultra), Reaction: BeamVeto_* (tight..ultra)
};

struct Tables {
  std::vector<Row> rows;
  void AddBeam(const std::string& prefix, double pGeV,
               long long Nbeam, long long trig2_tight, long long trig2_fit,
               long long trig2_wide, long long trig2_ultra)
  {
    Row r; r.kind=Kind::kBeam; r.prefix=prefix; r.pGeV=pGeV; r.N_beam=Nbeam;
    r.num[0]=trig2_tight; r.num[1]=trig2_fit; r.num[2]=trig2_wide; r.num[3]=trig2_ultra;
    rows.push_back(r);
  }
  void AddReaction(const std::string& prefix, double pGeV,
                   long long Nbeam, long long veto_tight, long long veto_fit,
                   long long veto_wide, long long veto_ultra)
  {
    Row r; r.kind=Kind::kReaction; r.prefix=prefix; r.pGeV=pGeV; r.N_beam=Nbeam;
    r.num[0]=veto_tight; r.num[1]=veto_fit; r.num[2]=veto_wide; r.num[3]=veto_ultra;
    rows.push_back(r);
  }
};

// ---------------- 스샷 값 입력 구역 ----------------
// (스샷에서 읽은 값으로 채워 넣음)
static void FillFromScreens(Tables& T){
  // ===== π− beam : E45_Nov_beam_xxx.root =====
  T.AddBeam("beam", 0.98,  869012,  76,  13,  13,  13);
  T.AddBeam("beam", 1.02,  864891,  75,  57,  57,  57);
  T.AddBeam("beam", 1.05,  863179, 111,  90,  90,  90);
  T.AddBeam("beam", 1.10,  866307, 125,  52,  21,  21);
  T.AddBeam("beam", 1.15,  866191,  78,  78,  78,  78);

  // ===== π+ beam : E45_Nov_beamplus_xxx.root =====
  T.AddBeam("beamplus", 0.98, 868879,  72,  55,  55,  55);
  T.AddBeam("beamplus", 1.02, 864368, 101,  60,  60,  60);
  T.AddBeam("beamplus", 1.05, 866461, 153,  86,  86,  86);
  T.AddBeam("beamplus", 1.10, 866806, 170,  75,  75,  75);
  T.AddBeam("beamplus", 1.15, 864974, 148,  67,  67,  67);

  // ===== π−p → π−π0p : E45_Nov_pi0p_xxx.root =====
  T.AddReaction("pi0p", 0.98, 865372, 13224, 20507, 27117, 31169);
  T.AddReaction("pi0p", 1.02, 865843,  7168, 11015, 15030, 18496);
  T.AddReaction("pi0p", 1.05, 866074, 13588, 20476, 27298, 31442);
  T.AddReaction("pi0p", 1.10, 865839, 14302, 21351, 27955, 32198);
  T.AddReaction("pi0p", 1.15, 864658, 14393, 21251, 28505, 32654);

  // ===== π−p → π−π+ n : E45_Nov_piplusn_xxx.root =====
  T.AddReaction("piplusn", 0.98, 865535,  7150, 11224, 14882, 18346);
  T.AddReaction("piplusn", 1.02, 865396, 13340, 20225, 27080, 31324);
  T.AddReaction("piplusn", 1.05, 866032,  7270, 11454, 15323, 18949);
  T.AddReaction("piplusn", 1.10, 863873,  7603, 12603, 16729, 20305);
  T.AddReaction("piplusn", 1.15, 866719,  8017, 12272, 17132, 20789);

  // ===== π+ p → π+π0 p : E45_Nov_piplus_pi0p_xxx.root =====
  T.AddReaction("piplus_pi0p", 0.98,  867711,  7335, 10516, 14182, 16296);
  T.AddReaction("piplus_pi0p", 1.02,  866499,  7519, 10840, 14593, 16767);
  T.AddReaction("piplus_pi0p", 1.05,  866798,  7927, 11321, 15065, 17357);
  T.AddReaction("piplus_pi0p", 1.10,  867704,  8368, 11771, 15523, 17870);
  T.AddReaction("piplus_pi0p", 1.15,  865450,  8808, 12033, 16018, 18615);

  // ===== π+ p → π+π+ n : E45_Nov_piplusplus_xxx.root =====
  T.AddReaction("piplusplus", 0.98,  867549, 12564, 15876, 22910, 25131);
  T.AddReaction("piplusplus", 1.02,  866103, 11821, 15555, 20387, 25401);
  T.AddReaction("piplusplus", 1.05,  867751, 12436, 16026, 24013, 26878);
  T.AddReaction("piplusplus", 1.10,  867658, 12860, 16707, 24406, 27342);
  T.AddReaction("piplusplus", 1.15,  866676, 13245, 17241, 25024, 27996);
}

// ---------------- 표 출력/CSV ----------------
static void PrintAndSave(const Tables& T){
  if (T.rows.empty()) {
    std::cout << "[WARN] FillFromScreens()에 데이터가 없습니다. AddBeam/AddReaction을 채우세요.\n";
    return;
  }

  // 채널·종류별로 모아서 정렬(모멘텀 오름차순)
  std::map<std::string, std::vector<const Row*>> beam, reac;
  for(const auto& r: T.rows){
    if(r.kind==Kind::kBeam) beam[r.prefix].push_back(&r);
    else                    reac[r.prefix].push_back(&r);
  }
  auto byP = [](const Row* a, const Row* b){ return a->pGeV < b->pGeV; };
  for(auto& kv: beam) std::sort(kv.second.begin(), kv.second.end(), byP);
  for(auto& kv: reac) std::sort(kv.second.begin(), kv.second.end(), byP);

  // 콘솔 헤더 (폭 확장: BEAM은 CI까지 보이도록)
  auto hdr_beam = [](){
    std::cout<<"\n--- BEAM: Beam Veto efficiency & BG@1M (eff [CI95]) ---\n";
    std::cout<<std::left<<std::setw(18)<<"Channel"
             <<std::setw(8)<<"p[GeV]"
             <<std::setw(32)<<"tight BG (eff [CI95])"
             <<std::setw(32)<<"fit BG (eff [CI95])"
             <<std::setw(32)<<"wide BG (eff [CI95])"
             <<std::setw(32)<<"ultra BG (eff [CI95])" <<"\n";
  };
  auto hdr_reac = [](){
    std::cout<<"\n--- REACTION: Overkill [%] (95% CI) ---\n";
    std::cout<<std::left<<std::setw(18)<<"Channel"
             <<std::setw(8)<<"p[GeV]"
             <<std::setw(24)<<"tight"
             <<std::setw(24)<<"fit"
             <<std::setw(24)<<"wide"
             <<std::setw(24)<<"ultra" <<"\n";
  };

  // CSV 파일
  std::ofstream csvB("E45_beam_summary.csv");
  csvB<<"channel,pGeV,N_beam,mask,trig2,eff_pct,eff_CI95_low,eff_CI95_high,bg_at_1M\n";
  std::ofstream csvR("E45_reaction_summary.csv");
  csvR<<"channel,pGeV,N_beam,mask,veto,ovk_pct,ovk_CI95_low,ovk_CI95_high\n";

  // ---- BEAM 출력 ----
  hdr_beam();
  for(const auto& kv: beam){
    for(const Row* r: kv.second){
      std::cout<<std::left<<std::setw(18)<<PrettyChannel(r->prefix)
               <<std::setw(8)<<std::fixed<<std::setprecision(2)<<r->pGeV;
      for(int i=0;i<4;++i){
        const long long trig2 = r->num[i];   // survivors (!veto)
        const long long n     = r->N_beam;
        const long long k_eff = n - trig2;   // success = veto-fired
        const CI ci = CP95(k_eff, n);
        const double eff = 100.0 * (double)k_eff / (double)n;
        const long long bg = (long long) std::llround( 1e6 * (double)trig2 / (double)n );

        std::ostringstream cell;
        // 예: "   87  (99.991% [99.989,99.993])"
        cell<<std::setw(6)<<bg<<"  ("
            <<std::fixed<<std::setprecision(3)<<eff<<"% "
            <<"["<<std::setprecision(3)<<100.0*ci.lo<<","<<100.0*ci.hi<<"])";
        std::cout<<std::setw(32)<<cell.str();

        csvB<<PrettyChannel(r->prefix)<<","<<r->pGeV<<","<<n<<","<<MaskName((Mask)i)
            <<","<<trig2<<","<<eff<<","<<(100.0*ci.lo)<<","<<(100.0*ci.hi)<<","<<bg<<"\n";
      }
      std::cout<<"\n";
    }
  }

  // ---- REACTION 출력 ----
  hdr_reac();
  for(const auto& kv: reac){
    for(const Row* r: kv.second){
      std::cout<<std::left<<std::setw(18)<<PrettyChannel(r->prefix)
               <<std::setw(8)<<std::fixed<<std::setprecision(2)<<r->pGeV;
      for(int i=0;i<4;++i){
        const long long veto = r->num[i]; // BeamVeto fired
        const long long n    = r->N_beam; // Beam
        const CI ci = CP95(veto, n);
        const double ovk = 100.0 * (double)veto / (double)n;

        std::ostringstream cell;
        cell<<std::fixed<<std::setprecision(3)<<ovk<<"% "
            <<"["<<std::setprecision(3)<<100.0*ci.lo<<","<<100.0*ci.hi<<"]";
        std::cout<<std::setw(24)<<cell.str();

        csvR<<PrettyChannel(r->prefix)<<","<<r->pGeV<<","<<n<<","<<MaskName((Mask)i)
            <<","<<veto<<","<<ovk<<","<<(100.0*ci.lo)<<","<<(100.0*ci.hi)<<"\n";
      }
      std::cout<<"\n";
    }
  }

  std::cout<<"\nCSV saved: E45_beam_summary.csv, E45_reaction_summary.csv\n";
}

// ---------------- 엔트리 포인트 ----------------
void E45_Summary_FromScreens(){
  Tables T;
  FillFromScreens(T);   // ← 스샷 수치 이미 채워 넣음
  PrintAndSave(T);
}

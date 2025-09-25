// BVH_veto_from_beammask.C
// 목적:
//   Beam 파일로부터 "Veto 마스크(빨강 ∪ 보라 = Beam TRUE)"를 만들고,
//   그 마스크를 적용해
//     (1) Beam-through 파일의 Veto 효율(%)  = 빔 억제율
//     (2) Overkill 파일의 오버킬 비율(%)     = false veto 비율
//   을 계산/출력.
//
// 분모/분자 정의(이전 BVH_3D_english 방식과 동일):
//   - 분모: (BH2 && BVH_U)인 이벤트 수
//   - 분자: 분모 이벤트 중 (BVH_D가 존재) && (선택된 마스크 셀과 매칭) 인 이벤트 수
//
// 사용법:
//   root -l
//   .L Veto_3D.C+
//   BVH_veto_from_beammask("E45_BVH1_60mm_54.root","E45_overkill_54.root", 0.10,0.04,0.04, 1, /*printPerBH2=*/1);
//
// 주의:
//   - 세그 개수는 BH2=15, BVH_U=15, BVH_D=54로 가정 (네 파일 확인 결과)
//   - threshold는 Beam 마스크 만들 때만 사용(thr_beam). Overkill에는 같은 마스크를 적용.
//   - 에너지 컷은 두 파일 모두 동일하게 적용(eBH2,eU,eD)

#include "Rtypes.h"
#include "TFile.h"
#include "TTree.h"
#include "TParticle.h"
#include "TH1.h"
#include "TH2F.h"
#include "TStyle.h"
#include "TString.h"
#include "TSystem.h"
#include "TInterpreter.h"
#include <vector>
#include <unordered_set>
#include <algorithm>
#include <iostream>
#include <iomanip>

static const int N_BH2  = 15;  // 0..14
static const int N_BVHU = 15;  // 0..14
static const int N_BVHD = 54;  // 0..53

// ---------- 유틸: 유니크 세그먼트 추출 ----------
static inline void get_unique_hits(const std::vector<TParticle>* v,
                                   int nmax, double cutMeV,
                                   std::vector<int>& out)
{
  out.clear();
  if(!v) return;
  std::unordered_set<int> s; s.reserve(8);
  for(const auto& p : *v){
    if(p.GetWeight() <= cutMeV) continue;  // weight=Edep[MeV]
    int id = p.GetMother(1);               // seg ID
    if(0<=id && id<nmax) s.insert(id);
  }
  out.assign(s.begin(), s.end());
  std::sort(out.begin(), out.end());
}

// ---------- 마스크(Beam TRUE) 만들기: per-BH2 U×D에서 thr 이상 셀만 1 ----------
struct MaskPack {
  std::vector<std::vector<char>> mask; // size: N_BH2 × (N_BVHU*N_BVHD)
};

static bool build_beam_mask(const char* beamFile,
                            double eBH2,double eU,double eD,
                            int thr_beam,
                            MaskPack& M)
{
  gSystem->Load("libPhysics");
  gInterpreter->GenerateDictionary("vector<TParticle>", "TParticle.h;vector");
  TH1::AddDirectory(kFALSE);

  if (gSystem->AccessPathName(beamFile)) {
    std::cerr << "[Error] File not found: " << beamFile << "\n";
    return false;
  }
  TFile* f=TFile::Open(beamFile,"READ");
  if(!f || f->IsZombie()){
    std::cerr << "[Error] Cannot open: " << beamFile << "\n";
    return false;
  }
  TTree* tr=(TTree*)f->Get("g4hyptpc");
  if(!tr){
    std::cerr << "[Error] no g4hyptpc in " << beamFile << "\n";
    f->Close(); return false;
  }
  if(!tr->GetBranch("BH2")||!tr->GetBranch("BVH_U")||!tr->GetBranch("BVH_D")){
    std::cerr << "[Error] missing branches in " << beamFile << "\n";
    tr->GetListOfBranches()->Print(); f->Close(); return false;
  }

  std::vector<TParticle> *BH2=nullptr,*BVHU=nullptr,*BVHD=nullptr;
  if(tr->SetBranchAddress("BH2",&BH2)<0 ||
     tr->SetBranchAddress("BVH_U",&BVHU)<0 ||
     tr->SetBranchAddress("BVH_D",&BVHD)<0){
    std::cerr << "[Error] SetBranchAddress failed\n";
    f->Close(); return false;
  }

  // per-BH2 카운트 히스토 준비 (파일 분리)
  std::vector<TH2F*> H(N_BH2,nullptr);
  for(int h=0;h<N_BH2;++h){
    TH2F* hh=new TH2F(Form("hUD_beam_%d",h),Form("beam BH2=%d;BVH_U;BVH_D",h),
                      N_BVHU,-0.5,N_BVHU-0.5, N_BVHD,-0.5,N_BVHD-0.5);
    hh->SetDirectory(nullptr);
    H[h]=hh;
  }

  // 첫 패스: U×D 빈도 채우기 (H,U,D 모두 존재하는 이벤트만)
  const Long64_t N = tr->GetEntries();
  std::vector<int> h,u,d; h.reserve(8); u.reserve(8); d.reserve(8);
  for(Long64_t i=0;i<N;++i){
    if(i && (i%200000==0)) std::cout<<"  [beam pass] "<<i<<"/"<<N<<"\r"<<std::flush;
    tr->GetEntry(i);
    get_unique_hits(BH2 ,N_BH2 ,eBH2,h);
    get_unique_hits(BVHU,N_BVHU,eU  ,u);
    get_unique_hits(BVHD,N_BVHD,eD  ,d);
    if(h.empty()||u.empty()||d.empty()) continue;
    for(int hh:h){
      TH2F* HH=H[hh]; if(!HH) continue;
      for(int uu:u){ if(uu<0||uu>=N_BVHU) continue;
        for(int dd:d){ if(dd<0||dd>=N_BVHD) continue; HH->Fill(uu,dd); }
      }
    }
  }
  std::cout << "\n";

  // 마스크 구성
  const int CELLS = N_BVHU*N_BVHD;
  M.mask.assign(N_BH2,std::vector<char>(CELLS,0));
  for(int hh=0;hh<N_BH2;++hh){
    TH2F* HH=H[hh]; if(!HH) continue;
    for(int bx=1;bx<=HH->GetNbinsX();++bx){
      for(int by=1;by<=HH->GetNbinsY();++by){
        if(HH->GetBinContent(bx,by) >= thr_beam){
          int uu=bx-1, dd=by-1;
          int idx = uu*N_BVHD + dd;
          if(0<=idx && idx<CELLS) M.mask[hh][idx]=1;
        }
      }
    }
  }

  f->Close();
  return true;
}

// ---------- 마스크 적용하여 (denom, num) 계산 ----------
struct VetoStats {
  Long64_t denom_all=0;   // BH2 && U
  Long64_t num_all=0;     // (mask match && D 존재)
  std::vector<Long64_t> denom_h, num_h; // per-BH2
};

static bool compute_veto_rate(const char* fname,
                              double eBH2,double eU,double eD,
                              const MaskPack& M,
                              VetoStats& S,
                              bool printPerBH2=false)
{
  gSystem->Load("libPhysics");
  gInterpreter->GenerateDictionary("vector<TParticle>", "TParticle.h;vector");

  if (gSystem->AccessPathName(fname)) {
    std::cerr << "[Error] File not found: " << fname << "\n";
    return false;
  }
  TFile* f=TFile::Open(fname,"READ");
  if(!f || f->IsZombie()){
    std::cerr << "[Error] Cannot open: " << fname << "\n";
    return false;
  }
  TTree* tr=(TTree*)f->Get("g4hyptpc");
  if(!tr){
    std::cerr << "[Error] no g4hyptpc in " << fname << "\n";
    f->Close(); return false;
  }
  std::vector<TParticle> *BH2=nullptr,*BVHU=nullptr,*BVHD=nullptr;
  tr->SetBranchAddress("BH2",&BH2);
  tr->SetBranchAddress("BVH_U",&BVHU);
  tr->SetBranchAddress("BVH_D",&BVHD);

  S.denom_all=0; S.num_all=0;
  S.denom_h.assign(N_BH2,0);
  S.num_h.assign(N_BH2,0);

  std::vector<int> h,u,d; h.reserve(8); u.reserve(8); d.reserve(8);
  const Long64_t N = tr->GetEntries();

  for(Long64_t i=0;i<N;++i){
    if(i && (i%200000==0)) std::cout<<"  [apply "<<gSystem->BaseName(fname)<<"] "<<i<<"/"<<N<<"\r"<<std::flush;
    tr->GetEntry(i);
    get_unique_hits(BH2 ,N_BH2 ,eBH2,h);
    get_unique_hits(BVHU,N_BVHU,eU  ,u);
    get_unique_hits(BVHD,N_BVHD,eD  ,d);

    const bool hasH=!h.empty(), hasU=!u.empty(), hasD=!d.empty();
    if(!(hasH && hasU)) continue; // 분모 조건
    S.denom_all++;

    bool global_match=false;
    for(int hh:h){
      S.denom_h[hh]++;
      bool matched_h=false;
      if(hasD){
        for(int uu:u){
          const int base = uu * N_BVHD;
          for(int dd:d){
            // 마스크 체크
            if(M.mask[hh][base+dd]){
              matched_h=true;
              global_match=true;
              break;
            }
          }
          if(matched_h) break;
        }
      }
      if(matched_h) S.num_h[hh]++;
    }
    if(global_match) S.num_all++;
  }
  std::cout << "\n";
  f->Close();

  if(printPerBH2){
    std::cout << "\nPer-BH2 veto rates (denom = events with BH2=h && U):\n";
    std::cout << "  h :  numerator / denominator  =  rate(%)\n";
    for(int hidx=0; hidx<N_BH2; ++hidx){
      if(S.denom_h[hidx]>0){
        double r = 100.0 * (double)S.num_h[hidx]/(double)S.denom_h[hidx];
        std::cout << " " << std::setw(2) << hidx << " :  "
                  << std::setw(9) << S.num_h[hidx] << " / " << std::setw(9) << S.denom_h[hidx]
                  << "  =  " << std::setw(7) << std::setprecision(3) << r << "\n";
      }else{
        std::cout << " " << std::setw(2) << hidx << " :  "
                  << std::setw(9) << 0 << " / " << std::setw(9) << 0
                  << "  =      n/a\n";
      }
    }
  }
  return true;
}

// ---------- 드라이버: 전체 실행 ----------
void BVH_veto_from_beammask(const char* beamFile, const char* sigFile,
                            double eBH2=0.10, double eU=0.04, double eD=0.04,
                            int thr_beam=1,
                            int printPerBH2=0)
{
  gStyle->SetOptStat(0);
  std::cout << "[Config] ecuts(MeV): BH2="<<eBH2<<" U="<<eU<<" D="<<eD
            << " | Beam-threshold="<<thr_beam<<"\n";

  // 1) Beam에서 Veto 마스크 생성 (빨강 ∪ 보라)
  MaskPack MB;
  std::cout << "[Step] Build beam mask from " << beamFile << "\n";
  if(!build_beam_mask(beamFile, eBH2, eU, eD, thr_beam, MB)){
    std::cerr << "[Abort] build_beam_mask failed.\n"; return;
  }

  // 2) Beam 파일에 적용 → 빔 억제율
  std::cout << "[Step] Apply mask to BEAM: " << beamFile << "\n";
  VetoStats SB;
  if(!compute_veto_rate(beamFile, eBH2, eU, eD, MB, SB, printPerBH2)){
    std::cerr << "[Abort] compute_veto_rate(beam) failed.\n"; return;
  }

  // 3) Overkill 파일에 적용 → 오버킬 비율(false veto)
  std::cout << "[Step] Apply mask to OVERKILL: " << sigFile << "\n";
  VetoStats SS;
  if(!compute_veto_rate(sigFile, eBH2, eU, eD, MB, SS, printPerBH2)){
    std::cerr << "[Abort] compute_veto_rate(signal) failed.\n"; return;
  }

  // 4) 요약
  auto pct = [](Long64_t num, Long64_t den)->double{
    return (den>0) ? 100.0*(double)num/(double)den : 0.0;
  };

  std::cout << "\n==================== Summary (mask = Beam TRUE) ====================\n";
  std::cout << "Beam file      : " << beamFile << "\n";
  std::cout << "  denom_all(B) : " << SB.denom_all << " (BH2 && U)\n";
  std::cout << "  num_all(B)   : " << SB.num_all   << " (mask match && D)\n";
  std::cout << "  VETO efficiency on Beam-through = "
            << std::fixed << std::setprecision(3) << pct(SB.num_all, SB.denom_all) << " %\n";
  std::cout << "-------------------------------------------------------------------\n";
  std::cout << "Signal file    : " << sigFile << "\n";
  std::cout << "  denom_all(S) : " << SS.denom_all << " (BH2 && U)\n";
  std::cout << "  num_all(S)   : " << SS.num_all   << " (mask match && D)\n";
  std::cout << "  Overkill ratio (false veto on signal) = "
            << std::fixed << std::setprecision(3) << pct(SS.num_all, SS.denom_all) << " %\n";
  std::cout << "===================================================================\n";
}

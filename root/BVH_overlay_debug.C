// BVH_3D_overlay_merged_out_verbose.C
// 콘솔에 (BH2, BVH_U, BVH_D) 좌표를 출력 + CSV 저장
// print_mode: 0=none, 1=sample (per BH2 최대 max_print개), 2=all
// one_based: true면 좌표를 1부터 표기
//
// 사용:
// .L BVH_3D_overlay_merged_out_verbose.C+
// BVH_3D_overlay_out("E45_BVH1_60mm_54.root","E45_overkill_54.root", 0.10,0.04,0.04, 1, 1, 1, 30, false);

#include "Rtypes.h"
#include "TFile.h"
#include "TTree.h"
#include "TParticle.h"
#include "TH1.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TString.h"
#include "TBox.h"
#include "TLine.h"
#include "TPad.h"
#include "TROOT.h"
#include "TInterpreter.h"
#include "TSystem.h"
#include <vector>
#include <unordered_set>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>

static const int N_BH2  = 15;  // 0..14
static const int N_BVHU = 15;  // 0..14
static const int N_BVHD = 54;  // 0..53

static inline void get_unique_hits(const std::vector<TParticle>* v,
                                   int nmax, double cutMeV,
                                   std::vector<int>& out){
  out.clear();
  if(!v) return;
  std::unordered_set<int> s; s.reserve(8);
  for(const auto& p:*v){
    if(p.GetWeight()<=cutMeV) continue;
    int id=p.GetMother(1);
    if(0<=id && id<nmax) s.insert(id);
  }
  out.assign(s.begin(), s.end());
  std::sort(out.begin(), out.end());
}

struct CountPack { std::vector<TH2F*> hUD; Long64_t entries=0; };
struct MaskPack  { std::vector<std::vector<char>> mask; };

static bool build_counts(const char* fname,double eBH2,double eU,double eD,CountPack& C){
  gSystem->Load("libPhysics");
  gInterpreter->GenerateDictionary("vector<TParticle>", "TParticle.h;vector");
  TH1::AddDirectory(kFALSE);

  if (gSystem->AccessPathName(fname)) { std::cerr<<"[Error] File not found: "<<fname<<"\n"; return false; }
  TFile* f=TFile::Open(fname,"READ");
  if(!f || f->IsZombie()){ std::cerr<<"[Error] Cannot open: "<<fname<<"\n"; return false; }
  TTree* tr=(TTree*)f->Get("g4hyptpc");
  if(!tr){ std::cerr<<"[Error] no g4hyptpc in "<<fname<<"\n"; f->Close(); return false; }
  if(!tr->GetBranch("BH2")||!tr->GetBranch("BVH_U")||!tr->GetBranch("BVH_D")){
    std::cerr<<"[Error] missing branches in "<<fname<<"\n"; tr->GetListOfBranches()->Print(); f->Close(); return false;
  }
  std::vector<TParticle> *BH2=nullptr,*BVHU=nullptr,*BVHD=nullptr;
  if(tr->SetBranchAddress("BH2",&BH2)<0||tr->SetBranchAddress("BVH_U",&BVHU)<0||tr->SetBranchAddress("BVH_D",&BVHD)<0){
    std::cerr<<"[Error] SetBranchAddress failed\n"; f->Close(); return false;
  }
  C.hUD.assign(N_BH2,nullptr);
  for(int h=0;h<N_BH2;++h){
    TH2F* H=new TH2F(Form("hUD_%s_bh2_%d",gSystem->BaseName(fname),h),
                     Form("%s BH2=%d;BVH_U;BVH_D",fname,h),
                     N_BVHU,-0.5,N_BVHU-0.5, N_BVHD,-0.5,N_BVHD-0.5);
    H->SetDirectory(nullptr);
    C.hUD[h]=H;
  }
  Long64_t N=tr->GetEntries(); C.entries=N;
  std::vector<int> h,u,d; h.reserve(8); u.reserve(8); d.reserve(8);
  for(Long64_t i=0;i<N;++i){
    tr->GetEntry(i);
    get_unique_hits(BH2 ,N_BH2 ,eBH2,h);
    get_unique_hits(BVHU,N_BVHU,eU  ,u);
    get_unique_hits(BVHD,N_BVHD,eD  ,d);
    if(h.empty()||u.empty()||d.empty()) continue;
    for(int hh:h){ if(hh<0||hh>=N_BH2) continue; TH2F* H=C.hUD[hh]; if(!H) continue;
      for(int uu:u){ if(uu<0||uu>=N_BVHU) continue;
        for(int dd:d){ if(dd<0||dd>=N_BVHD) continue; H->Fill(uu,dd); }
      }
    }
  }
  f->Close(); return true;
}

static void build_mask_from_counts(const CountPack& C,int thr,MaskPack& M){
  const int CELLS=N_BVHU*N_BVHD;
  M.mask.assign(N_BH2,std::vector<char>(CELLS,0));
  for(int h=0;h<N_BH2;++h){
    TH2F* H=(h<(int)C.hUD.size())?C.hUD[h]:nullptr; if(!H) continue;
    for(int bx=1;bx<=H->GetNbinsX();++bx){
      for(int by=1;by<=H->GetNbinsY();++by){
        if(H->GetBinContent(bx,by)>=thr){
          int u=bx-1,d=by-1; int idx=u*N_BVHD+d;
          if(0<=idx && idx<CELLS) M.mask[h][idx]=1;
        }
      }
    }
  }
}

static TH2F* make_sum_hist(TH2F* Hb, TH2F* Hs, const char* name){
  TH2F* Hsum=nullptr;
  if(Hb){ Hsum=(TH2F*)Hb->Clone(name); Hsum->SetDirectory(nullptr); if(Hs) Hsum->Add(Hs); }
  else if(Hs){ Hsum=(TH2F*)Hs->Clone(name); Hsum->SetDirectory(nullptr); }
  if(Hsum){ Hsum->GetXaxis()->SetTitle("BVH_U Seg"); Hsum->GetYaxis()->SetTitle("BVH_D Seg"); }
  return Hsum;
}

void BVH_3D_overlay_out(const char* beamFile,
                        const char* sigFile,
                        double eBH2=0.10,double eU=0.04,double eD=0.04,
                        int thr_beam=1,int thr_sig=1,
                        int print_mode=1, int max_print=50, bool one_based=false)
{
  gSystem->Load("libPhysics");
  gInterpreter->GenerateDictionary("vector<TParticle>", "TParticle.h;vector");
  gStyle->SetOptStat(0);
  TH1::AddDirectory(kFALSE);

  CountPack B,S;
  if(!build_counts(beamFile,eBH2,eU,eD,B)){ std::cerr<<"[Error] build_counts(beam) failed\n"; return; }
  if(!build_counts(sigFile ,eBH2,eU,eD,S)){ std::cerr<<"[Error] build_counts(sig) failed\n";  return; }

  MaskPack MB,MS;
  build_mask_from_counts(B,thr_beam,MB);
  build_mask_from_counts(S,thr_sig ,MS);

  // 파일명
  TString baseB = gSystem->BaseName(beamFile);
  TString baseS = gSystem->BaseName(sigFile);
  TString outBeam  = Form("beam_union_%s_thr%d.csv", baseB.Data(), thr_beam);
  TString outGreen = Form("green_only_%s_thr%d_vs_%s_thr%d.csv",
                          baseS.Data(), thr_sig, baseB.Data(), thr_beam);

  std::ofstream fBeam(outBeam.Data());
  std::ofstream fGreen(outGreen.Data());
  if(!fBeam){  std::cerr<<"[Error] cannot open "<<outBeam<<"\n"; }
  if(!fGreen){ std::cerr<<"[Error] cannot open "<<outGreen<<"\n"; }
  if(fBeam)  fBeam  << "BH2,BVH_U,BVH_D\n";
  if(fGreen) fGreen << "BH2,BVH_U,BVH_D\n";

  const int CELLS=N_BVHU*N_BVHD;
  size_t cntBeam=0, cntGreen=0;

  // 콘솔 헤더
  if(print_mode>0){
    std::cout << "\n=== Beam ∪ Magenta (Beam=TRUE) coordinates ===\n";
  }

  for(int h=0; h<N_BH2; ++h){
    const auto& mB = MB.mask[h];
    const auto& mS = MS.mask[h];
    if((int)mB.size()!=CELLS || (int)mS.size()!=CELLS){
      std::cerr<<"[Warn] mask size mismatch at BH2="<<h<<"\n";
      continue;
    }

    int printedBeam=0, printedGreen=0;
    if(print_mode>0){
      std::cout << "[BH2="<<(one_based? h+1:h)<<"]\n";
    }

    // (1) Beam TRUE (빨강 + 보라의 합집합)
    for(int u=0; u<N_BVHU; ++u){
      for(int d=0; d<N_BVHD; ++d){
        int idx=u*N_BVHD+d;
        bool b=mB[idx];
        bool s=mS[idx];

        if(b){
          int H = one_based? h+1:h;
          int U = one_based? u+1:u;
          int D = one_based? d+1:d;
          if(fBeam) fBeam << H << "," << U << "," << D << "\n";
          ++cntBeam;

          if(print_mode==2 || (print_mode==1 && printedBeam<max_print)){
            std::cout << "  (Beam) (BH2,U,D)=("<<H<<","<<U<<","<<D<<")\n";
            ++printedBeam;
          }
        }
      }
    }

    // (2) Green only (Signal TRUE && !Beam)
    if(print_mode>0){
      std::cout << "  -- Green-only --\n";
    }
    for(int u=0; u<N_BVHU; ++u){
      for(int d=0; d<N_BVHD; ++d){
        int idx=u*N_BVHD+d;
        bool b=mB[idx];
        bool s=mS[idx];

        if(s && !b){
          int H = one_based? h+1:h;
          int U = one_based? u+1:u;
          int D = one_based? d+1:d;
          if(fGreen) fGreen << H << "," << U << "," << D << "\n";
          ++cntGreen;

          if(print_mode==2 || (print_mode==1 && printedGreen<max_print)){
            std::cout << "  (GreenOnly) (BH2,U,D)=("<<H<<","<<U<<","<<D<<")\n";
            ++printedGreen;
          }
        }
      }
    }

    if(print_mode==1){
      if(printedBeam>=max_print)  std::cout << "  ... Beam more omitted ...\n";
      if(printedGreen>=max_print) std::cout << "  ... Green-only more omitted ...\n";
    }
  }

  if(fBeam)  fBeam.close();
  if(fGreen) fGreen.close();

  std::cout << "\n[OUT] Beam ∪ Magenta -> " << outBeam  << "  ("<<cntBeam<<" rows)\n";
  std::cout << "[OUT] Green-only      -> " << outGreen << "  ("<<cntGreen<<" rows)\n";

  // (옵션) 그림도 보고 싶으면 아래 호출(주석 해제):
  // ---- quick_view(B,S,MB,MS);  // 필요하면 시각화 함수 만들어서 붙일 수 있음
}

// BVH_3D_overlay_merged.C
//
// 목적:
//   - 두 ROOT 파일(E45_BVH1_60mm_54.root = Beam-through, E45_overkill_54.root = Overkill)을
//     각각 읽어 BH2별 U×D(Upstream×Downstream) 2D 빈도 히스토를 만든 뒤,
//     같은 캔버스에 빨강(Beam), 초록(Overkill), 보라(겹침) 박스로 오버레이.
//   - 배경 히트맵은 두 파일의 U×D 빈도를 **합산**해서 보여줌(시각적 참고용).
//   - 3개씩 5개 캔버스(= 15개 BH2)로 출력.
//   - threshold를 beam/signal 각각 독립적으로 설정할 수 있음 (기본 1,1).
//
// 사용법:
//   root -l
//   .L BVH_3D_overlay_merged.C+
//   BVH_3D_overlay("E45_BVH1_60mm_54.root","E45_overkill_54.root",
//                  0.10,0.04,0.04, 1, 1);
//
// 필요 브랜치(파일): BH2(vector<TParticle>), BVH_U(vector<TParticle>), BVH_D(vector<TParticle>)
// 세그먼트 수(이 파일 가정): BH2=15, BVH_U=15, BVH_D=54
//
// 메모/안전:
//   - vector<TParticle> 딕셔너리 자동 생성
//   - 파일/트리/브랜치 존재 확인
//   - 히스토그램은 파일과 분리(SetDirectory(nullptr)) → 파일 close 후 안전
//   - 인덱스/마스크 범위 보호

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

static const int N_BH2  = 15;  // BH2 segments (0..14)
static const int N_BVHU = 15;  // BVH_U segments (0..14)
static const int N_BVHD = 54;  // BVH_D segments (0..53)

static inline void get_unique_hits(const std::vector<TParticle>* v,
                                   int nmax, double cutMeV,
                                   std::vector<int>& out)
{
  out.clear();
  if(!v) return;
  std::unordered_set<int> s; s.reserve(8);
  for(const auto& p : *v){
    if(p.GetWeight() <= cutMeV) continue;      // weight = Edep[MeV]
    int id = p.GetMother(1);                   // seg ID
    if(0<=id && id<nmax) s.insert(id);
  }
  out.assign(s.begin(), s.end());
  std::sort(out.begin(), out.end());
}

struct CountPack {
  std::vector<TH2F*> hUD;     // per-BH2 counts
  Long64_t entries = 0;
};

struct MaskPack {
  std::vector<std::vector<char>> mask; // per-BH2, size = N_BVHU*N_BVHD
};

static bool build_counts(const char* fname,
                         double ecutBH2MeV, double ecutUMeV, double ecutDMeV,
                         CountPack& C)
{
  // 딕셔너리 보장
  gSystem->Load("libPhysics");
  gInterpreter->GenerateDictionary("vector<TParticle>", "TParticle.h;vector");
  TH1::AddDirectory(kFALSE); // 생성 히스토 파일과 분리

  if (gSystem->AccessPathName(fname)) {
    std::cerr << "[Error] File not found: " << fname << "\n";
    return false;
  }
  TFile* f = TFile::Open(fname,"READ");
  if(!f || f->IsZombie()){
    std::cerr << "[Error] Cannot open: " << fname << "\n";
    return false;
  }
  TTree* tr = (TTree*)f->Get("g4hyptpc");
  if(!tr){
    std::cerr << "[Error] TTree 'g4hyptpc' not found in " << fname << "\n";
    f->Close(); return false;
  }
  if(!tr->GetBranch("BH2") || !tr->GetBranch("BVH_U") || !tr->GetBranch("BVH_D")){
    std::cerr << "[Error] Missing branches (need BH2, BVH_U, BVH_D) in " << fname << "\n";
    tr->GetListOfBranches()->Print();
    f->Close(); return false;
  }

  std::vector<TParticle> *BH2=nullptr, *BVHU=nullptr, *BVHD=nullptr;
  if(tr->SetBranchAddress("BH2",&BH2)<0 ||
     tr->SetBranchAddress("BVH_U",&BVHU)<0 ||
     tr->SetBranchAddress("BVH_D",&BVHD)<0){
    std::cerr << "[Error] SetBranchAddress failed (dictionary/type)\n";
    f->Close(); return false;
  }

  // per-BH2 히스토 (파일과 분리)
  C.hUD.assign(N_BH2, nullptr);
  for(int h=0; h<N_BH2; ++h){
    TH2F* H = new TH2F(Form("hUD_%s_bh2_%d", gSystem->BaseName(fname), h),
                       Form("%s BH2=%d;BVH_U;BVH_D", fname, h),
                       N_BVHU,-0.5,N_BVHU-0.5,
                       N_BVHD,-0.5,N_BVHD-0.5);
    H->SetDirectory(nullptr);
    C.hUD[h] = H;
  }

  const Long64_t N = tr->GetEntries();
  C.entries = N;
  std::vector<int> hitsH, hitsU, hitsD;
  hitsH.reserve(8); hitsU.reserve(8); hitsD.reserve(8);

  for(Long64_t i=0; i<N; ++i){
    tr->GetEntry(i);
    get_unique_hits(BH2 ,N_BH2 , ecutBH2MeV, hitsH);
    get_unique_hits(BVHU,N_BVHU, ecutUMeV  , hitsU);
    get_unique_hits(BVHD,N_BVHD, ecutDMeV  , hitsD);
    if(hitsH.empty() || hitsU.empty() || hitsD.empty()) continue;

    for(int h : hitsH){
      if(h<0 || h>=N_BH2) continue;
      TH2F* H = C.hUD[h];
      if(!H) continue;
      for(int u : hitsU){
        if(u<0 || u>=N_BVHU) continue;
        for(int d : hitsD){
          if(d<0 || d>=N_BVHD) continue;
          H->Fill(u,d);
        }
      }
    }
  }
  f->Close();
  return true;
}

static void build_mask_from_counts(const CountPack& C, int thr, MaskPack& M){
  const int CELLS = N_BVHU * N_BVHD;
  M.mask.assign(N_BH2, std::vector<char>(CELLS, 0));
  for(int h=0; h<N_BH2; ++h){
    TH2F* H = (h<(int)C.hUD.size()) ? C.hUD[h] : nullptr;
    if(!H) continue;
    for(int bx=1; bx<=H->GetNbinsX(); ++bx){
      for(int by=1; by<=H->GetNbinsY(); ++by){
        if(H->GetBinContent(bx,by) >= thr){
          int u = bx-1, d = by-1;
          int idx = u*N_BVHD + d;
          if(idx>=0 && idx<CELLS) M.mask[h][idx] = 1;
        }
      }
    }
  }
}

// 배경 히트맵 = Beam+Signal 합
static TH2F* make_sum_hist(TH2F* Hb, TH2F* Hs, const char* name){
  TH2F* Hsum = nullptr;
  if(Hb){
    Hsum = (TH2F*)Hb->Clone(name);
    Hsum->SetDirectory(nullptr);
    if(Hs) Hsum->Add(Hs);
  }else if(Hs){
    Hsum = (TH2F*)Hs->Clone(name);
    Hsum->SetDirectory(nullptr);
  }
  if(Hsum){
    Hsum->GetXaxis()->SetTitle("BVH_U Seg");
    Hsum->GetYaxis()->SetTitle("BVH_D Seg");
  }
  return Hsum;
}

void BVH_3D_overlay(const char* beamFile,
                    const char* sigFile,
                    double ecutBH2MeV = 0.10,
                    double ecutUMeV   = 0.04,
                    double ecutDMeV   = 0.04,
                    int    thr_beam   = 1,
                    int    thr_sig    = 1)
{
  gSystem->Load("libPhysics");
  gInterpreter->GenerateDictionary("vector<TParticle>", "TParticle.h;vector");
  gStyle->SetOptStat(0);
  TH1::AddDirectory(kFALSE);

  // 1) 두 파일에서 per-BH2 2D count 만들기
  CountPack B, S;
  if(!build_counts(beamFile, ecutBH2MeV, ecutUMeV, ecutDMeV, B)){ std::cerr<<"[Error] build_counts(beam) failed.\n"; return; }
  if(!build_counts(sigFile , ecutBH2MeV, ecutUMeV, ecutDMeV, S)){ std::cerr<<"[Error] build_counts(sig) failed.\n";  return; }

  // 2) 두 파일 각각의 마스크 만들기
  MaskPack MB, MS;
  build_mask_from_counts(B, thr_beam, MB);
  build_mask_from_counts(S, thr_sig , MS);

  // 3) BH2를 3개씩 묶어 총 5개 캔버스로 그리기
  TLine grid; grid.SetLineStyle(kDotted); grid.SetLineColor(kGray+1);

  for(int h=0; h<N_BH2; ++h){
    if(h % 3 == 0){
      // 새 캔버스 (가로 3개)
      TCanvas* c = new TCanvas(Form("c_bh2_%02d_to_%02d", h, std::min(h+2,N_BH2-1)),
                               Form("BH2 Segments %d to %d", h, std::min(h+2,N_BH2-1)),
                               1500, 500);
      c->Divide(3,1);
    }

    // 현재 캔버스 = ROOT의 마지막 생성 캔버스
    TCanvas* ccur = (TCanvas*)gPad->GetCanvas();
    ccur->cd((h%3)+1);

    // 배경 = 합산 히스토(있으면), 없으면 Signal만/Beam만
    TH2F* Hb = (h<(int)B.hUD.size()) ? B.hUD[h] : nullptr;
    TH2F* Hs = (h<(int)S.hUD.size()) ? S.hUD[h] : nullptr;
    TH2F* Hsum = make_sum_hist(Hb, Hs, Form("hUD_sum_bh2_%d", h));
    if(!Hsum){ std::cerr << "[Warn] No hist to draw at BH2="<<h<<"\n"; continue; }

    Hsum->Draw("COLZ");
    gPad->Update();

    // 격자선
    double xmn = Hsum->GetXaxis()->GetXmin(), xmx = Hsum->GetXaxis()->GetXmax();
    double ymn = Hsum->GetYaxis()->GetXmin(), ymx = Hsum->GetYaxis()->GetXmax();
    for(int i=0;i<=N_BVHU;++i){
      double x = Hsum->GetXaxis()->GetBinLowEdge(i+1);
      grid.DrawLine(x, ymn, x, ymx);
    }
    for(int j=0;j<=N_BVHD;++j){
      double y = Hsum->GetYaxis()->GetBinLowEdge(j+1);
      grid.DrawLine(xmn, y, xmx, y);
    }

    // 박스 오버레이: 빨강=Beam, 초록=Overkill, 보라=둘 다
    const int CELLS = N_BVHU*N_BVHD;
    const auto& mB = MB.mask[h];
    const auto& mS = MS.mask[h];
    if((int)mB.size()!=CELLS || (int)mS.size()!=CELLS){
      std::cerr << "[Warn] mask size mismatch at BH2="<<h<<"\n";
      continue;
    }

    for(int u=0; u<N_BVHU; ++u){
      for(int d=0; d<N_BVHD; ++d){
        bool b = mB[u*N_BVHD + d];
        bool s = mS[u*N_BVHD + d];
        if(!b && !s) continue;

        double x1=u-0.5, x2=u+0.5;
        double y1=d-0.5, y2=d+0.5;
        TBox* box = new TBox(x1,y1,x2,y2);
        box->SetFillStyle(0);
        box->SetLineWidth(2);
        if(b && s)      box->SetLineColor(kMagenta); // 겹침
        else if(b)      box->SetLineColor(kRed);     // Beam only
        else            box->SetLineColor(kGreen+2); // Overkill only
        box->Draw("SAME");
      }
    }
    gPad->RedrawAxis();
  }

  std::cout << "[Info] Done. Beam(thr="<<thr_beam<<")=red, Overkill(thr="<<thr_sig<<")=green, both=magenta.\n";
}

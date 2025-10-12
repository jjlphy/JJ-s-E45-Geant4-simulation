// VP1_BH2_Study_v4.C
//  - PNG 저장 기본 끔
//  - 캔버스 강제 업데이트
//  - VP1 판정: seg==1/1001 OR |z - (-337)| < zTol (2→10→50mm 단계적 폴백)
//  - BH2 세그: X좌표(폭14mm, 15분할, 중심X=-10mm)로 역추정

#include <TFile.h>
#include <TTree.h>
#include <TParticle.h>
#include <TCanvas.h>
#include <TH2F.h>
#include <TBox.h>
#include <TLine.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TInterpreter.h>
#include <vector>
#include <cmath>
#include <iostream>
#include <limits>

namespace V4 {

// ===== 설정 =====
constexpr double BoxCx = 0.0;     // mm
constexpr double BoxCy = 12.0;    // mm
constexpr double BoxW  = 136.0;   // mm
constexpr double BoxH  = 112.0;   // mm
constexpr double BoxHX = BoxW*0.5;
constexpr double BoxHY = BoxH*0.5;

constexpr int    BH2_NSEG = 15;   // 0..14
constexpr double BH2_W    = 14.0; // mm
constexpr double BH2_Xc   = -10.0;
constexpr double BH2_Ymin = -100.0, BH2_Ymax=100.0;

constexpr double VP1_Z    = -337.0; // mm

constexpr int    Nbin = 200;
constexpr double Xmin=-500, Xmax=500, Ymin=-500, Ymax=500;

inline int getSegId(const TParticle& p){
  int sc = p.GetStatusCode();
  if(sc>=0 && sc<100000) return sc;
  int fm = p.GetFirstMother();
  if(fm>=0 && fm<100000) return fm;
  return static_cast<int>(std::lround(p.Px()));
}

inline int BH2_seg_fromX(double x){
  const double totalW = BH2_NSEG*BH2_W;
  const double x0 = BH2_Xc - 0.5*totalW;
  int seg = static_cast<int>(std::floor((x - x0)/BH2_W));
  if(seg<0) seg=0; if(seg>=BH2_NSEG) seg=BH2_NSEG-1;
  return seg;
}

inline void drawBox(){
  auto box=new TBox(BoxCx-BoxHX,BoxCy-BoxHY,BoxCx+BoxHX,BoxCy+BoxHY);
  box->SetFillStyle(0); box->SetLineColor(kRed+1); box->SetLineWidth(3); box->Draw("same");
}

inline void drawBH2Overlay(){
  const double totalW = BH2_NSEG*BH2_W;
  double x = BH2_Xc - 0.5*totalW;
  std::vector<double> centers; centers.reserve(BH2_NSEG);
  for(int i=0;i<=BH2_NSEG;i++){
    auto ln=new TLine(x,BH2_Ymin,x,BH2_Ymax);
    ln->SetLineColor(kGray+2); ln->SetLineStyle(3); ln->Draw("same");
    if(i<BH2_NSEG){ centers.push_back(x+0.5*BH2_W); }
    x += BH2_W;
  }
  TLatex tx; tx.SetTextSize(0.02);
  for(int s=0;s<BH2_NSEG;s++) tx.DrawLatex(centers[s], BH2_Ymax*0.9, Form("%d",s));
}

inline void countInOut(const TH2* h,double xmin,double xmax,double ymin,double ymax,
                       Long64_t& in,Long64_t& out){
  in=out=0;
  for(int ix=1; ix<=h->GetNbinsX(); ++ix){
    double xc=h->GetXaxis()->GetBinCenter(ix);
    for(int iy=1; iy<=h->GetNbinsY(); ++iy){
      double yc=h->GetYaxis()->GetBinCenter(iy);
      double c=h->GetBinContent(ix,iy);
      if(c<=0) continue;
      if(xc>=xmin && xc<=xmax && yc>=ymin && yc<=ymax) in += (Long64_t)c;
      else out += (Long64_t)c;
    }
  }
}

} // namespace V4

void VP1_BH2_Study(const char* fname="../../E45_with_SCH.root",
                   const char* treeName="g4hyptpc",
                   bool savePNG=false)  // 기본: PNG 저장 안 함
{
  using namespace V4;
  gStyle->SetOptStat(0);

  // vector<TParticle> 딕셔너리
  gSystem->Load("libEG");
  gInterpreter->GenerateDictionary("vector<TParticle>","TParticle.h;vector");

  // 파일/트리
  TFile* f=TFile::Open(fname,"READ");
  if(!f||f->IsZombie()){ std::cerr<<"[ERR] open "<<fname<<" failed\n"; return; }
  TTree* tr = dynamic_cast<TTree*>(f->Get(treeName));
  if(!tr){ std::cerr<<"[ERR] tree "<<treeName<<" not found\n"; f->Close(); return; }

  std::vector<TParticle>* BH2=nullptr;
  std::vector<TParticle>* VP=nullptr;

  bool hasBH2 = tr->GetBranch("BH2");
  bool hasVP  = tr->GetBranch("VP");

  if(hasBH2) tr->SetBranchAddress("BH2",&BH2); else std::cerr<<"[WRN] no 'BH2' branch\n";
  if(hasVP ) tr->SetBranchAddress("VP" ,&VP ); else std::cerr<<"[WRN] no 'VP'  branch\n";

  const Long64_t N = tr->GetEntries();
  std::cout<<"[INFO] Entries = "<<N<<"\n";

  // 히스토그램
  TH2F* hVP1_all  = new TH2F("hVP1_all","VP1 XY (all);X [mm];Y [mm]", Nbin,Xmin,Xmax, Nbin,Ymin,Ymax);
  TH2F* hVP1_bh49 = new TH2F("hVP1_bh49","VP1 XY (BH2 seg4-9);X [mm];Y [mm]", Nbin,Xmin,Xmax, Nbin,Ymin,Ymax);
  TH2F* hBH2_xy   = new TH2F("hBH2_xy","BH2 XY;X [mm];Y [mm]", Nbin,-150,150, Nbin,-100,100);

  // 디버그용 Z 분포( VP1 판정 보조 )
  double vpZmin=+std::numeric_limits<double>::infinity();
  double vpZmax=-std::numeric_limits<double>::infinity();

  Long64_t nBH2_any=0, nBH2_3to10=0, nBH2_4to9=0, nBH2_4to9_AND_VP1=0;

  for(Long64_t i=0;i<N;++i){
    tr->GetEntry(i);

    bool bh2_any=false, bh2_3to10=false, bh2_4to9=false;
    bool vp1_any=false;

    // BH2: 좌표→세그
    if(BH2){
      bool fa=false,f310=false,f49=false;
      for(const auto& p: *BH2){
        const double x=p.Vx(), y=p.Vy();
        hBH2_xy->Fill(x,y);
        int seg = BH2_seg_fromX(x);
        if(!fa){ bh2_any=true; fa=true; }
        if(!f310 && (seg>=3 && seg<=10)){ bh2_3to10=true; f310=true; }
        if(!f49  && (seg>=4 && seg<= 9)){ bh2_4to9 =true; f49 =true; }
        if(fa && f310 && f49) break;
      }
    }

    // VP1: ID 또는 Z-근접(순차 확대)로 판정
    if(VP){
      bool filled=false;
      for(const auto& p: *VP){
        const int seg = getSegId(p);
        const double x=p.Vx(), y=p.Vy(), z=p.Vz();

        if(z<vpZmin) vpZmin=z;
        if(z>vpZmax) vpZmax=z;

        bool isVP1 = (seg==1 || seg==1001);
        double dz = std::abs(z - VP1_Z);
        if(!isVP1) {
          if(dz<=2.0)      isVP1=true;
          else if(dz<=10.0) isVP1=true;
          else if(dz<=50.0) isVP1=true;
        }
        if(!isVP1) continue;

        hVP1_all->Fill(x,y);
        vp1_any=true;
        if(bh2_4to9) hVP1_bh49->Fill(x,y);
        if(filled) break; filled=true;
      }
    }

    if(bh2_any)   ++nBH2_any;
    if(bh2_3to10) ++nBH2_3to10;
    if(bh2_4to9)  ++nBH2_4to9;
    if(bh2_4to9 && vp1_any) ++nBH2_4to9_AND_VP1;
  }

  auto pct=[&](Long64_t k){ return N>0 ? 100.0*double(k)/double(N) : 0.0; };

  std::cout<<"\n================ SUMMARY ================\n";
  std::cout<<"Total events                    : "<<N<<"\n";
  std::cout<<"BH2 any-hit (N, %)             : "<<nBH2_any<<"  ("<<pct(nBH2_any) <<" %)\n";
  std::cout<<"BH2 seg 3–10 pass (N, %)       : "<<nBH2_3to10<<"  ("<<pct(nBH2_3to10)<<" %)\n";
  std::cout<<"BH2 seg 4–9  pass (N, %)       : "<<nBH2_4to9 <<"  ("<<pct(nBH2_4to9) <<" %)\n";
  std::cout<<"BH2 seg 4–9  ∧ VP1 pass (N, %) : "<<nBH2_4to9_AND_VP1<<"  ("<<pct(nBH2_4to9_AND_VP1)<<" %)\n";

  // 히스토 엔트리 확인(캔버스 비어보일 때 원인 파악용)
  std::cout<<"[DEBUG] Entries: hBH2_xy="<<hBH2_xy->GetEntries()
           <<", hVP1_all="<<hVP1_all->GetEntries()
           <<", hVP1_bh49="<<hVP1_bh49->GetEntries()<<"\n";
  if(hasVP){
    std::cout<<"[DEBUG] VP z-range seen: ["<<vpZmin<<", "<<vpZmax<<"] mm (VP1_Z="<<VP1_Z<<")\n";
  }

  // 오버레이 in/out
  const double xmin=BoxCx-BoxHX, xmax=BoxCx+BoxHX;
  const double ymin=BoxCy-BoxHY, ymax=BoxCy+BoxHY;
  auto ratio=[](Long64_t a,Long64_t b){ double s=a+b; return s>0? 100.0*double(a)/s : 0.0; };

  Long64_t nin1=0,nout1=0,nin3=0,nout3=0;
  countInOut(hVP1_all ,xmin,xmax,ymin,ymax,nin1,nout1);
  countInOut(hVP1_bh49,xmin,xmax,ymin,ymax,nin3,nout3);

  std::cout<<"\n---- Overlay box: center(0,12), size 136x112 mm ----\n";
  std::cout<<"[1] VP1(all)    : inside="<<nin1<<"  outside="<<nout1
           <<"  (inside "<<ratio(nin1,nout1)<<" %)\n";
  std::cout<<"[3] VP1(BH2 4–9): inside="<<nin3<<"  outside="<<nout3
           <<"  (inside "<<ratio(nin3,nout3)<<" %)\n";
  std::cout<<"===================================================\n\n";

  // 캔버스(강제 업데이트; PNG 저장 X)
  TCanvas* c1=new TCanvas("c1","VP1 (all)",900,800);
  hVP1_all->Draw("COLZ"); drawBox(); c1->Modified(); c1->Update();

  TCanvas* c3=new TCanvas("c3","VP1 (BH2 seg4-9)",900,800);
  hVP1_bh49->Draw("COLZ"); drawBox(); c3->Modified(); c3->Update();

  TCanvas* c5=new TCanvas("c5","BH2 XY",900,700);
  hBH2_xy->Draw("COLZ"); drawBH2Overlay(); c5->Modified(); c5->Update();

  if(savePNG){
    c1->SaveAs("VP1_all_xy.png");
    c3->SaveAs("VP1_BH2seg4to9_xy.png");
    c5->SaveAs("BH2_xy_with_segments.png");
  }

  f->Close();
  std::cout<<"[DONE]\n";
}

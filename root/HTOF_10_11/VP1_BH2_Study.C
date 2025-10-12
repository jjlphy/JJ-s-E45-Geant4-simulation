// VP1_BH2_Study.C  (ROOT 6.32 호환: TTimer::SingleShot 사용 안 함)
// - BH2 세그 구간을 인자로 지정(segLo, segHi) → 4–9, 3–10, 4–10 등 쉽게 변경
// - macOS/화면 비표시 이슈 회피: GL painter 끄고 강제 Update
// - 좌표가 로컬일 수 있어 BH2/KVC/HTOF는 로컬→글로벌 보정 후 사용
// - VP 브랜치가 없거나 맞지 않아도 BH2↔(KVC 우선, 없으면 HTOF) 직선으로
//   Z = VP1_Z 에서의 교점(x,y)을 투영해 VP1 분포를 채움
// - PNG 저장은 기본 끔
/*root -l
.L VP1_BH2_Study.C+
VP1_BH2_Study("../../E45_with_SCH.root","g4hyptpc",  /*segLo=*/3, /*segHi=*/10, /*VP1_Z=*/-337.0,/*savePNG=*/
              /*segLo=*/3, /*segHi=*/10,  // ← 여기만 바꾸면 3–10 등으로 OK
              /*VP1_Z=*/-337.0,
          //    /*savePNG=*/false) 
              

#include <TFile.h>
#include <TTree.h>
#include <TParticle.h>
#include <TCanvas.h>
#include <TH2F.h>
#include <TBox.h>
#include <TLine.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TInterpreter.h>
#include <vector>
#include <cmath>
#include <iostream>

namespace V7F {

// ====== 기하 상수(글로벌, mm) ======
struct Ctr { double x,y,z; };
constexpr Ctr C_BH2  = {-10.0,  0.0, -560.0};  // DCGeomParam_E72_0 기준
constexpr Ctr C_HTOF = {  0.0, 12.0,    0.0};
constexpr Ctr C_KVC  = {160.0,  0.0,  600.0};

// BH2 세그 정의(글로벌 X)
constexpr int    BH2_NSEG = 15;
constexpr double BH2_W    = 14.0;              // seg width
constexpr double BH2_Xc   = -10.0;             // BH2 center X
constexpr double BH2_Ymin = -100.0, BH2_Ymax = 100.0;

// 오버레이 박스
constexpr double BoxCx = 0.0;                  // mm
constexpr double BoxCy = 12.0;
constexpr double BoxW  = 136.0;
constexpr double BoxH  = 112.0;
constexpr double BoxHX = BoxW*0.5, BoxHY = BoxH*0.5;

// 히스토 범위
constexpr int    Nbin = 200;
constexpr double Xmin=-500, Xmax=500, Ymin=-500, Ymax=500;

inline bool looksLocalZ(double z_local, double z_center){
  // 로컬이면 z≈0, 글로벌이면 ≈ 중심 z
  return (std::fabs(z_local) < 80.0 && std::fabs(z_local - z_center) > 120.0);
}
inline void toGlobal(const TParticle& p, const Ctr& C,
                     double& xg,double& yg,double& zg, bool& usedLocal)
{
  if(looksLocalZ(p.Vz(), C.z)){
    xg = p.Vx()+C.x; yg = p.Vy()+C.y; zg = p.Vz()+C.z; usedLocal=true;
  }else{
    xg = p.Vx();     yg = p.Vy();     zg = p.Vz();     usedLocal=false;
  }
}

inline int BH2_seg_fromX(double xg){
  const double tot = BH2_NSEG*BH2_W;
  const double x0  = BH2_Xc - 0.5*tot;    // 왼쪽 끝
  int seg = (int)std::floor((xg - x0)/BH2_W);
  if(seg<0) seg=0; if(seg>=BH2_NSEG) seg=BH2_NSEG-1;
  return seg;
}

inline bool lineIntersectAtZ(double x0,double y0,double z0,
                             double x1,double y1,double z1,
                             double Zgoal, double& xi,double& yi)
{
  const double dz = z1 - z0;
  if(std::fabs(dz) < 1e-6) return false;
  const double t = (Zgoal - z0)/dz;
  if(t < -0.5 || t > 1.5) return false;   // 과도한 외삽 방지
  xi = x0 + t*(x1 - x0);
  yi = y0 + t*(y1 - y0);
  return true;
}

// Draw helpers
inline void drawBox(){
  auto box=new TBox(BoxCx-BoxHX, BoxCy-BoxHY, BoxCx+BoxHX, BoxCy+BoxHY);
  box->SetFillStyle(0); box->SetLineColor(kRed+1); box->SetLineWidth(3); box->Draw("same");
}
inline void drawBH2Overlay(){
  const double tot = BH2_NSEG*BH2_W;
  double x = BH2_Xc - 0.5*tot;
  std::vector<double> centers; centers.reserve(BH2_NSEG);
  for(int i=0;i<=BH2_NSEG;i++){
    auto ln=new TLine(x,BH2_Ymin,x,BH2_Ymax);
    ln->SetLineColor(kGray+2); ln->SetLineStyle(3); ln->Draw("same");
    if(i<BH2_NSEG) centers.push_back(x+0.5*BH2_W);
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

inline void forceRefresh(){
  gPad->Modified(); gPad->Update();
  gSystem->ProcessEvents();
  gSystem->Sleep(80);
}

} // namespace V7F


void VP1_BH2_Study(const char* fname="../../E45_with_SCH.root",
                   const char* treeName="g4hyptpc",
                   int segLo=4, int segHi=9,     // ★ BH2 세그 구간 인자
                   double VP1_Z=-337.0,          // VP1 평면 Z (필요시 변경)
                   bool savePNG=false)           // PNG 저장 여부(기본 끔)
{
  using namespace V7F;

  // 렌더러/배치 설정: 창 비표시 이슈 회피
  gROOT->SetBatch(kFALSE);
  gStyle->SetCanvasPreferGL(kFALSE);

  // 딕셔너리
  gSystem->Load("libEG");
  gInterpreter->GenerateDictionary("vector<TParticle>","TParticle.h;vector");

  // 파일/트리
  TFile* f=TFile::Open(fname,"READ");
  if(!f||f->IsZombie()){ std::cerr<<"[ERR] open "<<fname<<" failed\n"; return; }
  TTree* tr = dynamic_cast<TTree*>(f->Get(treeName));
  if(!tr){ std::cerr<<"[ERR] tree "<<treeName<<" not found\n"; f->Close(); return; }

  // 브랜치
  std::vector<TParticle>* BH2 = nullptr;
  std::vector<TParticle>* KVC = nullptr;
  std::vector<TParticle>* HTOF= nullptr;
  std::vector<TParticle>* VP  = nullptr;

  bool hasBH2 = tr->GetBranch("BH2");
  bool hasKVC = tr->GetBranch("KVC");
  bool hasHTOF= tr->GetBranch("HTOF");
  bool hasVP  = tr->GetBranch("VP");

  if(hasBH2) tr->SetBranchAddress("BH2",&BH2); else std::cerr<<"[WRN] no 'BH2'\n";
  if(hasKVC) tr->SetBranchAddress("KVC",&KVC);
  if(hasHTOF)tr->SetBranchAddress("HTOF",&HTOF);
  if(hasVP ) tr->SetBranchAddress("VP" ,&VP );

  const Long64_t N = tr->GetEntries();
  std::cout<<"[INFO] Entries="<<N<<"  (BH2/KVC/HTOF/VP="
           <<(hasBH2?"Y":"N")<<"/"<<(hasKVC?"Y":"N")<<"/"
           <<(hasHTOF?"Y":"N")<<"/"<<(hasVP?"Y":"N")<<")\n";
  std::cout<<"[INFO] BH2 segment window = ["<<segLo<<","<<segHi<<"] (inclusive)\n";

  // 히스토
  TH2F* hVP1_all  = new TH2F("hVP1_all","VP1 XY (all; projected/globalized);X [mm];Y [mm]",
                              Nbin,Xmin,Xmax, Nbin,Ymin,Ymax);
  TH2F* hVP1_sel  = new TH2F("hVP1_sel","VP1 XY (BH2 selected; projected/globalized);X [mm];Y [mm]",
                              Nbin,Xmin,Xmax, Nbin,Ymin,Ymax);
  TH2F* hBH2_xy   = new TH2F("hBH2_xy","BH2 XY (globalized);X [mm];Y [mm]",
                              Nbin,-150,150, Nbin,-100,100);

  // 카운트
  Long64_t nBH2_any=0, nBH2_win=0, nBH2_3to10=0, nBH2_4to9=0, nBH2_win_AND_VP1=0;
  Long64_t nProjOK=0, nVPbranchUsed=0;

  for(Long64_t i=0;i<N;++i){
    tr->GetEntry(i);

    bool bh2_any=false, bh2_inwin=false;
    bool vp1_any=false;

    // ---- BH2 ----
    const TParticle* bh2_first = nullptr;
    double bxg=0,byg=0,bzg=0; bool ul=false;
    if(BH2 && !BH2->empty()){
      bh2_first = &(*BH2)[0];
      toGlobal(*bh2_first, C_BH2, bxg,byg,bzg, ul);

      bool flagged_any=false, flagged_win=false;
      bool f310=false, f49=false; // 통계용
      for(const auto& p: *BH2){
        double xg,yg,zg; bool tmp=false;
        toGlobal(p, C_BH2, xg,yg,zg, tmp);
        hBH2_xy->Fill(xg,yg);
        int seg = BH2_seg_fromX(xg);

        if(!flagged_any){ bh2_any=true; flagged_any=true; }
        if(!flagged_win && seg>=segLo && seg<=segHi){ bh2_inwin=true; flagged_win=true; }
        if(!f310 && seg>=3 && seg<=10) f310=true;
        if(!f49  && seg>=4 && seg<= 9) f49=true;

        if(flagged_any && flagged_win && f310 && f49) break;
      }
      if(f310) ++nBH2_3to10;     // 전체 이벤트 대비 통계
      if(f49)  ++nBH2_4to9;
    }

    // ---- VP1 채우기: (1) VP 브랜치 사용 시도 ----
    bool filled_all=false, filled_sel=false;
    if(VP && !VP->empty()){
      for(const auto& p: *VP){
        // 보통 VP는 글로벌 저장 → 그대로 사용, z만 확인
        const double x=p.Vx(), y=p.Vy(), z=p.Vz();
        if(std::fabs(z - VP1_Z) > 5.0) continue;
        hVP1_all->Fill(x,y); vp1_any=true; filled_all=true;
        if(bh2_inwin){ hVP1_sel->Fill(x,y); filled_sel=true; }
        ++nVPbranchUsed;
        break;
      }
    }

    // ---- (2) 투영: BH2 ↔ (KVC 우선, 없으면 HTOF) ----
    if((!filled_all) || (bh2_inwin && !filled_sel)){
      const TParticle* downRaw = nullptr; Ctr C_down = C_KVC;
      if(KVC && !KVC->empty()){ downRaw=&(*KVC)[0]; C_down=C_KVC; }
      else if(HTOF && !HTOF->empty()){ downRaw=&(*HTOF)[0]; C_down=C_HTOF; }

      if(bh2_first && downRaw){
        double dxg,dyg,dzg; bool ul2=false;
        toGlobal(*downRaw, C_down, dxg,dyg,dzg, ul2);
        double xi,yi;
        if(lineIntersectAtZ(bxg,byg,bzg, dxg,dyg,dzg, VP1_Z, xi,yi)){
          hVP1_all->Fill(xi,yi); vp1_any=true;
          if(bh2_inwin) hVP1_sel->Fill(xi,yi);
          ++nProjOK;
        }
      }
    }

    if(bh2_any)   ++nBH2_any;
    if(bh2_inwin) ++nBH2_win;
    if(bh2_inwin && vp1_any) ++nBH2_win_AND_VP1;
  }

  auto pct=[&](Long64_t k){ return (N>0 ? 100.0*double(k)/double(N) : 0.0); };

  // ===== 출력 요약 =====
  std::cout<<"\n================ SUMMARY ================\n";
  std::cout<<"Total events                        : "<<N<<"\n";
  std::cout<<"BH2 any-hit (N, %)                 : "<<nBH2_any<<"  ("<<pct(nBH2_any) <<" %)\n";
  std::cout<<"BH2 seg ["<<segLo<<","<<segHi<<"] (N, %)  : "<<nBH2_win<<"  ("<<pct(nBH2_win) <<" %)\n";
  std::cout<<"BH2 seg 3–10 (N, %)                : "<<nBH2_3to10<<"  ("<<pct(nBH2_3to10)<<" %)\n";
  std::cout<<"BH2 seg 4–9  (N, %)                : "<<nBH2_4to9 <<"  ("<<pct(nBH2_4to9) <<" %)\n";
  std::cout<<"BH2 seg ["<<segLo<<","<<segHi<<"] ∧ VP1 (N, %): "<<nBH2_win_AND_VP1
           <<"  ("<<pct(nBH2_win_AND_VP1)<<" %)\n";
  std::cout<<"VP-branch used="<<nVPbranchUsed<<", Projected="<<nProjOK<<"\n";

  // 오버레이 상자 in/out
  auto ratio=[](Long64_t a,Long64_t b){ double s=a+b; return s>0? 100.0*double(a)/s : 0.0; };
  const double xmin=BoxCx-BoxHX, xmax=BoxCx+BoxHX;
  const double ymin=BoxCy-BoxHY, ymax=BoxCy+BoxHY;
  Long64_t nin_all=0, nout_all=0, nin_sel=0, nout_sel=0;
  countInOut(hVP1_all , xmin,xmax,ymin,ymax, nin_all,nout_all);
  countInOut(hVP1_sel , xmin,xmax,ymin,ymax, nin_sel,nout_sel);

  std::cout<<"\n---- Overlay (center (0,12), 136x112 mm) ----\n";
  std::cout<<"[VP1 all] inside="<<nin_all<<"  outside="<<nout_all
           <<"  (inside "<<ratio(nin_all,nout_all)<<" %)\n";
  std::cout<<"[VP1 BH2 ["<<segLo<<","<<segHi<<"]] inside="<<nin_sel<<"  outside="<<nout_sel
           <<"  (inside "<<ratio(nin_sel,nout_sel)<<" %)\n";
  std::cout<<"=============================================\n\n";

  // ===== 그림 (PNG 저장 기본 끔 + 강제 갱신) =====
  TCanvas* c1=new TCanvas("c1","VP1 XY (all)",900,800);
  hVP1_all->Draw("COLZ"); drawBox(); forceRefresh();
  if(savePNG) c1->SaveAs("VP1_all_xy.png");

  TCanvas* c2=new TCanvas("c2","VP1 XY (BH2 selected)",900,800);
  hVP1_sel->Draw("COLZ"); drawBox(); forceRefresh();
  if(savePNG) c2->SaveAs("VP1_selected_xy.png");

  TCanvas* c3=new TCanvas("c3","BH2 XY",900,700);
  hBH2_xy->Draw("COLZ"); drawBH2Overlay(); forceRefresh();
  if(savePNG) c3->SaveAs("BH2_xy_with_segments.png");

  f->Close();
  std::cout<<"[DONE]\n";
}

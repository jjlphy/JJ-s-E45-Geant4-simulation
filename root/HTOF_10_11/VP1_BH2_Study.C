// VP1_BH2_Study.C
// Build: root -l; .L VP1_BH2_Study.C+; VP1_BH2_Study("E45_VPHTOF_098.root","g4hyptpc",true)
//
// 본 매크로는 다음 가정을 사용합니다:
//  - TTree 이름: treeName (기본 "g4hyptpc")
//  - 가지(브랜치):
//      * "BH2" : std::vector<TParticle>  (세그먼트 히트들)
//      * "VP"  : std::vector<TParticle>  (VP1..VPn 히트들, TParticle의 fMother[0] 또는 fStatusCode 또는 Px()에 seg/copyNo 저장하는 관례 중 하나)
//  - 위치는 TParticle::Vx(), Vy(), Vz() [mm]에 저장
//  - 세그먼트 ID는 다음 우선순위로 추출: StatusCode -> FirstMother -> Px()
//    (환경에 맞춰 하나만 써도 되며, 필요 시 아래 getSegId()를 수정하세요)

#include <TFile.h>
#include <TTree.h>
#include <TParticle.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TBox.h>
#include <TLine.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TMath.h>
#include <vector>
#include <iostream>
#include <cassert>

namespace {
  // --------- 설정값 ----------
  // 오버레이 박스: 중심 (0, 12 mm), 크기 136 x 112 mm
  constexpr double kBoxCx = 0.0;     // mm
  constexpr double kBoxCy = 12.0;    // mm
  constexpr double kBoxW  = 136.0;   // mm (X 크기)
  constexpr double kBoxH  = 112.0;   // mm (Y 크기)
  // 박스 반폭/반높이
  constexpr double kBoxHX = kBoxW*0.5;
  constexpr double kBoxHY = kBoxH*0.5;

  // BH2 세그먼트 폭 (mm) — DetSize_E72_0: Bh2SegWidth 15개 전부 14 mm
  const int    kNBH2Seg = 15;
  const double kBH2SegWidth[kNBH2Seg] =
    {14,14,14,14,14,14,14,14,14,14,14,14,14,14,14};
  // BH2 세그먼트 Y 범위 (대충 넉넉히)
  constexpr double kBH2Ymin = -100.0;
  constexpr double kBH2Ymax =  100.0;

  // 히스토그램 범위 (넉넉히)
  // VP는 HistParam에서 VPX/VPY가 ±500으로 정의되어 있어 그에 맞춤
  constexpr int    kNbin = 200;
  constexpr double kXmin = -500.0, kXmax = 500.0;
  constexpr double kYmin = -500.0, kYmax = 500.0;

  // 세그먼트 판별 범위
  inline bool inRange(int s, int lo, int hi){ return (s>=lo && s<=hi); }

  // 이벤트-중복(다중 히트) 대비: 이벤트당 1카운트만 세고 싶을 때 true
  constexpr bool kDistinctPerEventForCounts = true;

  // --------- 유틸 ----------
  // TParticle로부터 "세그먼트 ID"를 추출한다.
  // 실 코드 구현과 다르면 이 함수를 한 줄만 수정해서 맞추면 됩니다.
  int getSegId(const TParticle& p) {
    // 우선순위 1: StatusCode에 세그 ID 저장
    int sc = p.GetStatusCode();
    if (sc>=0 && sc<10000) return sc;

    // 우선순위 2: FirstMother를 세그 ID로 사용 (환경에 따라)
    int fm = p.GetFirstMother();
    if (fm>=0 && fm<10000) return fm;

    // 우선순위 3: Px()를 세그 ID로 사용 (정수로 캐스팅)
    int pxAsSeg = static_cast<int>(std::lround(p.Px()));
    if (pxAsSeg>=0 && pxAsSeg<10000) return pxAsSeg;

    // 마지막: 실패시 -1
    return -1;
  }

  // 박스 내부 여부
  bool insideBox(double x, double y){
    return (x >= (kBoxCx - kBoxHX) && x <= (kBoxCx + kBoxHX) &&
            y >= (kBoxCy - kBoxHY) && y <= (kBoxCy + kBoxHY));
  }

  // BH2 세그먼트 경계 X 좌표들 (중심 0에 대칭 배치 가정)
  void buildBH2Edges(std::vector<double>& edges, std::vector<double>& centers){
    edges.clear(); centers.clear();
    double total = 0.0;
    for(int i=0;i<kNBH2Seg;i++) total += kBH2SegWidth[i];
    double x0 = -0.5*total; // 왼쪽 끝
    edges.push_back(x0);
    double x = x0;
    for(int i=0;i<kNBH2Seg;i++){
      double next = x + kBH2SegWidth[i];
      centers.push_back(0.5*(x+next));
      edges.push_back(next);
      x = next;
    }
  }

  void drawBH2Overlay(){
    std::vector<double> edges, centers;
    buildBH2Edges(edges, centers);
    // 세그먼트 경계선
    for(size_t i=0;i<edges.size();++i){
      TLine* ln = new TLine(edges[i], kBH2Ymin, edges[i], kBH2Ymax);
      ln->SetLineColor(kGray+2);
      ln->SetLineStyle(3);
      ln->Draw("same");
    }
    // seg index 라벨
    TLatex tx;
    tx.SetTextSize(0.02);
    for(int seg=0; seg<kNBH2Seg; ++seg){
      tx.DrawLatex(centers[seg], kBH2Ymax*0.9, Form("%d",seg));
    }
  }

  void drawOverlayBox(){
    TBox* box = new TBox(kBoxCx-kBoxHX, kBoxCy-kBoxHY, kBoxCx+kBoxHX, kBoxCy+kBoxHY);
    box->SetFillStyle(0);
    box->SetLineColor(kRed+1);
    box->SetLineWidth(3);
    box->Draw("same");
  }

  // 히스토그램 안/밖 카운트
  void countInOut(const TH2* h, double xmin, double xmax, double ymin, double ymax,
                  Long64_t& nin, Long64_t& nout)
  {
    nin = nout = 0;
    int nx = h->GetNbinsX();
    int ny = h->GetNbinsY();
    for(int ix=1; ix<=nx; ++ix){
      double xc = h->GetXaxis()->GetBinCenter(ix);
      for(int iy=1; iy<=ny; ++iy){
        double yc = h->GetYaxis()->GetBinCenter(iy);
        double c  = h->GetBinContent(ix,iy);
        if(c<=0) continue;
        if(xc>=xmin && xc<=xmax && yc>=ymin && yc<=ymax) nin += static_cast<Long64_t>(c);
        else                                           nout+= static_cast<Long64_t>(c);
      }
    }
  }
}

void VP1_BH2_Study(const char* fname="../E45_VPHTOF_098.root",
                   const char* treeName="g4hyptpc",
                   bool savePNG=true)
{
  gStyle->SetOptStat(0);

  // 파일/트리 열기
  TFile* f = TFile::Open(fname,"READ");
  if(!f || f->IsZombie()){ std::cerr<<"[ERR] cannot open "<<fname<<"\n"; return; }

  TTree* tr = dynamic_cast<TTree*>(f->Get(treeName));
  if(!tr){ std::cerr<<"[ERR] TTree '"<<treeName<<"' not found.\n"; f->Close(); return; }

  // 브랜치 포인터
  std::vector<TParticle>* BH2 = nullptr;
  std::vector<TParticle>* VP  = nullptr;

  // 브랜치 설정 (없으면 경고만)
  if(tr->GetBranch("BH2")) tr->SetBranchAddress("BH2",&BH2);
  else std::cerr<<"[WRN] branch 'BH2' not found. BH2 관련 카운트는 0이 됩니다.\n";

  if(tr->GetBranch("VP"))  tr->SetBranchAddress("VP",&VP);
  else std::cerr<<"[WRN] branch 'VP' not found. VP 관련 히스토그램은 비게 됩니다.\n";

  const Long64_t nent = tr->GetEntries();
  std::cout<<"[INFO] Entries = "<<nent<<"\n";

  // 히스토그램
  TH2F* hVP1_all   = new TH2F("hVP1_all","VP1 XY (all events);X [mm];Y [mm]", kNbin,kXmin,kXmax, kNbin,kYmin,kYmax);
  TH2F* hVP1_bh49  = new TH2F("hVP1_bh49","VP1 XY (BH2 seg4-9 pass);X [mm];Y [mm]", kNbin,kXmin,kXmax, kNbin,kYmin,kYmax);
  TH2F* hBH2_xy    = new TH2F("hBH2_xy","BH2 XY;X [mm];Y [mm]", kNbin,-150,150, kNbin,-100,100);

  // 카운트
  Long64_t nBH2_any=0, nBH2_3to10=0, nBH2_4to9=0;
  Long64_t nBH2_4to9_AND_VP1=0;

  // 이벤트 루프
  for(Long64_t i=0;i<nent;++i){
    tr->GetEntry(i);

    bool bh2_any=false, bh2_3to10=false, bh2_4to9=false;
    bool vp1_any=false;

    // --- BH2 판정 ---
    if(BH2){
      // 이벤트당 중복제거 위한 플래그
      bool flagged_any=false, flagged_3to10=false, flagged_4to9=false;

      for(const auto& p : *BH2){
        const int seg = getSegId(p);
        const double x = p.Vx();
        const double y = p.Vy();

        // 분포용 XY 채우기 (원하면 조건 추가)
        hBH2_xy->Fill(x,y);

        // 세그 판정
        if(!flagged_any){ bh2_any=true; flagged_any=true; }
        if(!flagged_3to10 && inRange(seg,3,10)){ bh2_3to10=true; flagged_3to10=true; }
        if(!flagged_4to9  && inRange(seg,4, 9)){ bh2_4to9 =true; flagged_4to9 =true; }

        if(kDistinctPerEventForCounts && flagged_any && flagged_3to10 && flagged_4to9) break;
      }
    }

    // --- VP1 히트 및 좌표 채우기 ---
    if(VP){
      // VP 세그는 copyNo(=seg)로 들어온다고 가정 → seg==1이 VP1
      bool filled_vp1_this_event=false;
      for(const auto& p : *VP){
        const int seg = getSegId(p);
        if(seg!=1) continue; // VP1만
        const double x = p.Vx();
        const double y = p.Vy();
        hVP1_all->Fill(x,y);
        vp1_any = true;
        filled_vp1_this_event = true;
        if(bh2_4to9){
          hVP1_bh49->Fill(x,y);
        }
        if(kDistinctPerEventForCounts && filled_vp1_this_event) break;
      }
    }

    // --- 카운트 집계 ---
    if(bh2_any)   ++nBH2_any;
    if(bh2_3to10) ++nBH2_3to10;
    if(bh2_4to9)  ++nBH2_4to9;
    if(bh2_4to9 && vp1_any) ++nBH2_4to9_AND_VP1;
  }

  // 비율 계산 (전체 이벤트 대비)
  auto pct = [&](Long64_t k){ return (nent>0 ? 100.0*double(k)/double(nent) : 0.0); };

  std::cout<<"\n================ SUMMARY ================\n";
  std::cout<<"Total events                    : "<<nent<<"\n";
  std::cout<<"BH2 any-hit (N, %)             : "<<nBH2_any<<"  ("<<pct(nBH2_any) <<" %)\n";
  std::cout<<"BH2 seg 3–10 pass (N, %)       : "<<nBH2_3to10<<"  ("<<pct(nBH2_3to10)<<" %)\n";
  std::cout<<"BH2 seg 4–9  pass (N, %)       : "<<nBH2_4to9 <<"  ("<<pct(nBH2_4to9) <<" %)\n";
  std::cout<<"BH2 seg 4–9  ∧ VP1 pass (N, %) : "<<nBH2_4to9_AND_VP1<<"  ("<<pct(nBH2_4to9_AND_VP1)<<" %)\n";

  // --- 오버레이 박스 내부/외부 카운트 (히스토그램 bin 합계 기준) ---
  const double xmin = kBoxCx - kBoxHX;
  const double xmax = kBoxCx + kBoxHX;
  const double ymin = kBoxCy - kBoxHY;
  const double ymax = kBoxCy + kBoxHY;

  Long64_t nin1=0, nout1=0, nin3=0, nout3=0;
  countInOut(hVP1_all , xmin,xmax,ymin,ymax, nin1,nout1);
  countInOut(hVP1_bh49, xmin,xmax,ymin,ymax, nin3,nout3);

  auto ratio = [](Long64_t a, Long64_t b){ double s=a+b; return (s>0? 100.0*double(a)/s : 0.0); };

  std::cout<<"\n---- Overlay box: center(0,12 mm), size 136x112 mm ----\n";
  std::cout<<"[1] VP1(all)    : inside="<<nin1<<"  outside="<<nout1
           <<"  (inside "<<ratio(nin1,nout1)<<" %)\n";
  std::cout<<"[3] VP1(BH2 4–9): inside="<<nin3<<"  outside="<<nout3
           <<"  (inside "<<ratio(nin3,nout3)<<" %)\n";
  std::cout<<"=======================================================\n\n";

  // ---------- 그림 그리기 ----------
  // (1) VP1 all
  TCanvas* c1 = new TCanvas("c1","VP1 (all)",900,800);
  hVP1_all->Draw("COLZ");
  drawOverlayBox();
  TLatex t1; t1.SetNDC(); t1.SetTextSize(0.035);
  t1.DrawLatex(0.15,0.93,"VP1 XY (all events)");
  if(savePNG) c1->SaveAs("VP1_all_xy.png");

  // (3) VP1 with BH2 seg4–9
  TCanvas* c3 = new TCanvas("c3","VP1 (BH2 seg4-9)",900,800);
  hVP1_bh49->Draw("COLZ");
  drawOverlayBox();
  TLatex t3; t3.SetNDC(); t3.SetTextSize(0.035);
  t3.DrawLatex(0.12,0.93,"VP1 XY (events with BH2 seg 4-9)");
  if(savePNG) c3->SaveAs("VP1_BH2seg4to9_xy.png");

  // (5) BH2 xy + 세그먼트 오버레이
  TCanvas* c5 = new TCanvas("c5","BH2 XY",900,700);
  hBH2_xy->Draw("COLZ");
  drawBH2Overlay();
  TLatex t5; t5.SetNDC(); t5.SetTextSize(0.035);
  t5.DrawLatex(0.18,0.93,"BH2 XY with segment boundaries (0..14)");
  if(savePNG) c5->SaveAs("BH2_xy_with_segments.png");

  // 끝
  f->Close();
  std::cout<<"[DONE]\n";
}

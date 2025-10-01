// HTOF_Summary_BH2_Cuts.C
// - 전체 이벤트 통계
// - BH2 통과 비율, BH2&HTOF>=1 비율, BH2&(HTOF>=1)&&(!{16,17,18,19}) 비율
// - 위 "깨끗한" 분모에 대해: Multiplicity, Hit pattern, Plane×Local 히스토그램
// - 분모에서  M>=2 개수/비율
// - 분모에서  M==2 인 이벤트 중 (2,3), (3,4), (4,5), (5,6) 개수/비율

#include <TFile.h>
#include <TTree.h>
#include <TParticle.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <vector>
#include <set>
#include <iostream>
#include <iomanip>
#include <TSystem.h>
#include <TInterpreter.h>

static const int    kNPlaneHTOF  = 8;
static const int    kSegPerPlane = 4;
static const int    kNSegHTOF    = kNPlaneHTOF * kSegPerPlane; // 32
static const double kPitch       = 68.0;   // [mm]
static const double kL           = 337.0;  // [mm] (radial)
static const double kHTOFx       = 0.0;
static const double kHTOFy       = 12.0;
static const double kHTOFz       = 0.0;

struct _DICT_ {
  _DICT_(){
    gSystem->Load("libPhysics");
    gInterpreter->GenerateDictionary("vector<TParticle>","TParticle.h;vector");
  }
} _dict_;

// ---- 좌표를 HTOF 중심기준으로 이동
inline void toLocalCenter(double x,double y,double z,double& lx,double& ly,double& lz){
  lx = x - kHTOFx;  ly = y - kHTOFy;  lz = z - kHTOFz;
}

// ---- plane i의 로컬 축 (ex: 가로, ez: 반경)
inline void planeAxes(int i, double& exx, double& exz, double& ezx, double& ezz){
  const double th = i * 2.0 * TMath::Pi() / kNPlaneHTOF; // i*45°
  exx =  TMath::Cos(th);  exz = -TMath::Sin(th);
  ezx =  TMath::Sin(th);  ezz =  TMath::Cos(th);
}

// ---- xloc → local(0..3)
inline int inferLocalFromX(double xloc, bool flip_local=false){
  double idx_f = xloc / kPitch + 1.5; // 중앙 오프셋(1.5)
  int j = TMath::Nint(idx_f);
  if(j<0) j=0; if(j>=kSegPerPlane) j=kSegPerPlane-1;
  if(flip_local) j = (kSegPerPlane-1)-j; // 0↔3, 1↔2
  return j;
}

// ---- (x,y,z) → seg(0..31)
inline bool mapHitToSeg(double x,double y,double z,
                        int L_sign,int plane_offset,bool flip_local,
                        int& plane,int& local,int& seg,double R_tol_mm)
{
  double lx,ly,lz; toLocalCenter(x,y,z,lx,ly,lz);
  const double R = std::hypot(lx,lz);
  if(R_tol_mm>0 && std::abs(R - kL) > R_tol_mm) return false;

  int best_i=-1; double best_abs=1e99, best_xloc=0.0;
  for(int i=0;i<kNPlaneHTOF;++i){
    double exx,exz,ezx,ezz; planeAxes(i,exx,exz,ezx,ezz);
    const double zloc = lx*ezx + lz*ezz;           // 반경
    const double diff = std::abs(zloc - (L_sign*kL));
    if(diff < best_abs){
      best_abs = diff; best_i = i;
      best_xloc = lx*exx + lz*exz;                 // 가로
    }
  }
  if(best_i<0) return false;

  plane = (best_i + plane_offset) % kNPlaneHTOF;
  if(plane<0) plane += kNPlaneHTOF;

  local = inferLocalFromX(best_xloc, flip_local);
  seg   = plane*kSegPerPlane + local;
  return true;
}

void HTOF_Summary_BH2_Cuts(const char* filename="E45_VP5.root",
                           const char* treename="g4hyptpc",
                           double R_tol_mm = 200.0, // 337±R_tol_mm 내 hit만 사용(0이면 비활성)
                           int    L_sign    = +1,   // plane 선택시 목표 zloc 부호(+1: +L, -1: -L)
                           int    plane_offset = 0, // GUI/정의가 한 칸 어긋나면 보정
                           bool   flip_local  = false,
                           int    mult_max    = 7   // Multiplicity 히스토그램 상한(0..mult_max)
){
  // 파일/트리
  TFile* f=TFile::Open(filename);
  if(!f||f->IsZombie()){ printf("Cannot open %s\n",filename); return; }
  TTree* T=(TTree*)f->Get(treename);
  if(!T){ printf("Tree %s not found\n",treename); return; }

  // 브랜치
  std::vector<TParticle>* BH2=nullptr;
  std::vector<TParticle>* HTOF=nullptr;
  if(!T->GetBranch("BH2") || !T->GetBranch("HTOF")){
    printf("Branch 'BH2' or 'HTOF' not found.\n"); return;
  }
  T->SetBranchAddress("BH2",&BH2);
  T->SetBranchAddress("HTOF",&HTOF);

  // 히스토그램 (분모: “BH2 && HTOF>=1 && not{16..19}” 에 대해서만)
  TH1D* hSeg_unique = new TH1D("hHTOF_seg_unique",
    "HTOF Hit Pattern (UNIQUE/event) — BH2 && HTOF>=1 && !{16..19};Seg ID;Counts",
    32,-0.5,31.5);
  TH2D* hPlaneLocal = new TH2D("hHTOF_plane_local",
    "Plane vs Local — BH2 && HTOF>=1 && !{16..19};Plane;Local",
    8,-0.5,7.5, 4,-0.5,3.5);
  TH1D* hMult = new TH1D("hHTOF_mult",
    "Multiplicity per event — BH2 && HTOF>=1 && !{16..19};# unique hit seg;Events",
    mult_max+1,-0.5,mult_max+0.5);

  // 카운터
  Long64_t N_total         = T->GetEntries();
  Long64_t N_bh2           = 0; // BH2 통과
  Long64_t N_bh2_htofAny   = 0; // BH2 && HTOF>=1
  Long64_t N_cleanDenom    = 0; // BH2 && HTOF>=1 && !{16..19}
  Long64_t N_clean_Mge2    = 0; // 위 분모에서 M>=2
  Long64_t N_clean_Meq2    = 0; // 위 분모에서 M==2

  // M==2인 이벤트에서의 페어 카운트
  Long64_t C23=0, C34=0, C45=0, C56=0;

  // 루프
  for(Long64_t ie=0; ie<N_total; ++ie){
    T->GetEntry(ie);
    if(!BH2 || !HTOF) continue;

    // (1) BH2 통과?
    if(BH2->empty()) continue;
    N_bh2++;

    // (2) HTOF 세그먼트 매핑
    std::set<int> segset;
    for(const auto& p : *HTOF){
      int plane,local,seg;
      if(!mapHitToSeg(p.Vx(),p.Vy(),p.Vz(),L_sign,plane_offset,flip_local,
                      plane,local,seg,R_tol_mm)) continue;
      if(seg>=0 && seg<kNSegHTOF) segset.insert(seg);
    }

    if(segset.empty()) continue;          // BH2지만 HTOF가 전혀 없음
    N_bh2_htofAny++;

    // (3) {16,17,18,19} 포함하면 제외
    bool has16to19 = segset.count(16) || segset.count(17) ||
                     segset.count(18) || segset.count(19);
    if(has16to19) continue;

    // ---- 여기부터가 "분모" ----
    N_cleanDenom++;

    // 히스토그램 채우기 (UNIQUE/event)
    int M = 0;
    for(int s : segset){
      hSeg_unique->Fill(s);
      hPlaneLocal->Fill(s/4, s%4);
      M++;
    }
    if(M>mult_max) M = mult_max;
    hMult->Fill(M);

    // M>=2 / M==2
    if(M>=2) N_clean_Mge2++;
    if(M==2){
      N_clean_Meq2++;
      // 페어 판정
      bool has23 = segset.count(2) && segset.count(3);
      bool has34 = segset.count(3) && segset.count(4);
      bool has45 = segset.count(4) && segset.count(5);
      bool has56 = segset.count(5) && segset.count(6);
      if(has23) C23++;
      if(has34) C34++;
      if(has45) C45++;
      if(has56) C56++;
    }
  }

  // 비율 계산 도우미
  auto pct = [](Long64_t num, Long64_t den)->double{
    return (den>0) ? (100.0 * double(num) / double(den)) : 0.0;
  };

  // ---- (1)~(4) 요약 출력 ----
  std::cout << std::fixed << std::setprecision(2);
  std::cout << "========== SUMMARY ==========\n";
  std::cout << "Total events                                   : " << N_total << "\n";
  std::cout << "BH2 pass                                       : " << N_bh2
            << "  (" << pct(N_bh2, N_total) << " % of total)\n";
  std::cout << "BH2 pass & HTOF >=1                            : " << N_bh2_htofAny
            << "  (" << pct(N_bh2_htofAny, N_bh2) << " % of BH2)\n";
  std::cout << "BH2 pass & HTOF >=1 & !{16,17,18,19}           : " << N_cleanDenom
            << "  (" << pct(N_cleanDenom, N_bh2) << " % of BH2)\n";

  // (3) 분모에서 Multiplicity>=2
  std::cout << "\n[Denominator = BH2 & HTOF>=1 & !{16..19}]\n";
  std::cout << "Multiplicity >= 2                              : " << N_clean_Mge2
            << "  (" << pct(N_clean_Mge2, N_cleanDenom) << " % of denominator)\n";

  // (4) 분모에서 M==2 중 각 페어
  std::cout << "Multiplicity == 2                              : " << N_clean_Meq2
            << "  (" << pct(N_clean_Meq2, N_cleanDenom) << " % of denominator)\n";
  std::cout << "  - pair (2,3)                                 : " << C23
            << "  (" << pct(C23, N_clean_Meq2) << " % of M==2)\n";
  std::cout << "  - pair (3,4)                                 : " << C34
            << "  (" << pct(C34, N_clean_Meq2) << " % of M==2)\n";
  std::cout << "  - pair (4,5)                                 : " << C45
            << "  (" << pct(C45, N_clean_Meq2) << " % of M==2)\n";
  std::cout << "  - pair (5,6)                                 : " << C56
            << "  (" << pct(C56, N_clean_Meq2) << " % of M==2)\n";
  std::cout << "=============================\n";

  // ---- 그림: 분모에 대해서만 히스토그램 시각화 ----
  TCanvas* c = new TCanvas("cHTOF_SUM","BH2-gated, no(16-19) — Distributions",1400,800);
  c->Divide(2,2);
  c->cd(1); hSeg_unique->Draw("hist");
  c->cd(2); hPlaneLocal->Draw("colz");
  c->cd(3); hMult->Draw("hist");
}

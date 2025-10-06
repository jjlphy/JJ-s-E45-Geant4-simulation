// HTOF_Summary_BH2_Cuts_Exclude15to20.C
// - BH2 게이트 & HTOF 세그먼트 매핑 유지
// - 배제: seg 15..20
// - 2D 히트패턴은 "이벤트 내 서로 다른 세그먼트 집합" 간 교차곱으로 채움
//   (같은 세그먼트 중복 히트는 1회만 기록)
// - Multiplicity(M) = 유일 세그 수, M>=2 통계와 % 출력 포함
//
// 사용 예:
//   root -l
//   .L HTOF_Summary_BH2_Cuts_Exclude15to20.C+
//   HTOF_Summary_BH2_Cuts_Exclude15to20("E45_VP5.root","g4hyptpc",200.0,+1,0,false,7);

#include <TFile.h>
#include <TTree.h>
#include <TParticle.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TSystem.h>
#include <TInterpreter.h>
#include <vector>
#include <set>
#include <iostream>
#include <iomanip>

static const int    kNPlaneHTOF  = 8;
static const int    kSegPerPlane = 4;
static const int    kNSegHTOF    = kNPlaneHTOF * kSegPerPlane; // 32
static const double kPitch       = 68.0;   // [mm]
static const double kL           = 337.0;  // [mm]
static const double kHTOFx       = 0.0;
static const double kHTOFy       = 12.0;
static const double kHTOFz       = 0.0;

enum : int {
  PDG_PiPlus  = +211,
  PDG_PiMinus = -211,
  PDG_Proton  = 2212
};

struct _DICT_ {
  _DICT_(){
    gSystem->Load("libPhysics");
    gInterpreter->GenerateDictionary("vector<TParticle>","TParticle.h;vector");
  }
} _dict_;

inline void toLocalCenter(double x,double y,double z,double& lx,double& ly,double& lz){
  lx = x - kHTOFx;  ly = y - kHTOFy;  lz = z - kHTOFz;
}
inline void planeAxes(int i, double& exx, double& exz, double& ezx, double& ezz){
  const double th = i * 2.0 * TMath::Pi() / kNPlaneHTOF; // i*45°
  exx =  TMath::Cos(th);  exz = -TMath::Sin(th);
  ezx =  TMath::Sin(th);  ezz =  TMath::Cos(th);
}
inline int inferLocalFromX(double xloc, bool flip_local=false){
  double idx_f = xloc / kPitch + 1.5;
  int j = TMath::Nint(idx_f);
  if(j<0) j=0; if(j>=kSegPerPlane) j=kSegPerPlane-1;
  if(flip_local) j = (kSegPerPlane-1)-j;
  return j;
}
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
    const double zloc = lx*ezx + lz*ezz;
    const double diff = std::abs(zloc - (L_sign*kL));
    if(diff < best_abs){
      best_abs = diff; best_i = i;
      best_xloc = lx*exx + lz*exz;
    }
  }
  if(best_i<0) return false;

  plane = (best_i + plane_offset) % kNPlaneHTOF;
  if(plane<0) plane += kNPlaneHTOF;

  local = inferLocalFromX(best_xloc, flip_local);
  seg   = plane*kSegPerPlane + local;
  return true;
}

void HTOF_Summary_BH2_Cuts_Exclude15to20(const char* filename="E45_VP5.root",
                           const char* treename="g4hyptpc",
                           double R_tol_mm = 200.0,
                           int    L_sign    = +1,
                           int    plane_offset = 0,
                           bool   flip_local  = false,
                           int    mult_max    = 7)
{
  // --- 파일/트리 ---
  TFile* f=TFile::Open(filename);
  if(!f||f->IsZombie()){ printf("Cannot open %s\n",filename); return; }
  TTree* T=(TTree*)f->Get(treename);
  if(!T){ printf("Tree %s not found\n",treename); return; }

  // --- 브랜치 ---
  std::vector<TParticle>* BH2 = nullptr;
  std::vector<TParticle>* HTOF = nullptr;
  T->SetBranchAddress("BH2",&BH2);
  T->SetBranchAddress("HTOF",&HTOF);

  // --- 히스토그램 ---
  TH1D* hSeg_unique = new TH1D("hHTOF_seg_unique",
    "HTOF Hit Pattern (UNIQUE/event, exclude 15-20);Seg ID;Counts",32,-0.5,31.5);
  TH2D* hPlaneLocal = new TH2D("hHTOF_plane_local",
    "Plane vs Local (exclude 15-20);Plane;Local",8,-0.5,7.5,4,-0.5,3.5);
  TH1D* hMult = new TH1D("hHTOF_mult",
    "Multiplicity per event (exclude 15-20);# unique hit seg;Events",
    mult_max+1,-0.5,mult_max+0.5);

  // --- 2D 히트패턴(32×32) ---
  TH2D* hPiPlus_vs_PiMinus = new TH2D("hPiPlus_vs_PiMinus",
    "HTOF (pi^{+} seg) vs (pi^{-} seg);#pi^{+} seg ID;#pi^{-} seg ID",
    32,-0.5,31.5, 32,-0.5,31.5);

  TH2D* hProton_vs_PiMinus = new TH2D("hProton_vs_PiMinus",
    "HTOF (p seg) vs (pi^{-} seg);#p seg ID;#pi^{-} seg ID",
    32,-0.5,31.5, 32,-0.5,31.5);

  // --- 카운터 ---
  Long64_t N_total=T->GetEntries();
  Long64_t N_bh2=0, N_bh2_htofAny=0, N_cleanDenom=0;
  Long64_t N_clean_Mge2=0, N_clean_Meq2=0;
  Long64_t C23=0,C34=0,C45=0,C56=0;

  for(Long64_t ie=0; ie<N_total; ++ie){
    T->GetEntry(ie);
    if(!BH2||!HTOF) continue;

    // BH2 게이트
    if(BH2->empty()) continue;
    N_bh2++;

    // --- 이벤트 내 "유일 세그먼트" 집합들(같은 세그 중복 제거) ---
    std::set<int> seg_all, seg_piPlus, seg_piMinus, seg_proton;

    // --- HTOF 입자 분류 & 세그먼트 매핑 ---
    for(const auto& p : *HTOF){
      int plane,local,seg;
      if(!mapHitToSeg(p.Vx(),p.Vy(),p.Vz(),
                      L_sign,plane_offset,flip_local,
                      plane,local,seg,R_tol_mm)) continue;
      if(seg<0 || seg>=kNSegHTOF) continue;

      seg_all.insert(seg);
      const int pdg = p.GetPdgCode();
      if(pdg==PDG_PiPlus)  seg_piPlus.insert(seg);
      if(pdg==PDG_PiMinus) seg_piMinus.insert(seg);
      if(pdg==PDG_Proton)  seg_proton.insert(seg);
    }

    if(seg_all.empty()) continue;
    N_bh2_htofAny++;

    // --- 배제 세그먼트(15..20) 포함 시 스킵 ---
    bool hasExcl=false;
    for(int s : {15,16,17,18,19,20})
      if(seg_all.count(s)) { hasExcl=true; break; }
    if(hasExcl) continue;

    N_cleanDenom++;

    // --- Multiplicity (유일 세그 수) & 기본 분포 채움 ---
    int M=0;
    for(int s: seg_all){
      hSeg_unique->Fill(s);
      hPlaneLocal->Fill(s/4, s%4);
      M++;
    }
    if(M>mult_max) M=mult_max;
    hMult->Fill(M);
    if(M>=2) N_clean_Mge2++;
    if(M==2){
      N_clean_Meq2++;
      if(seg_all.count(2)&&seg_all.count(3)) C23++;
      if(seg_all.count(3)&&seg_all.count(4)) C34++;
      if(seg_all.count(4)&&seg_all.count(5)) C45++;
      if(seg_all.count(5)&&seg_all.count(6)) C56++;
    }

    // === 2D 히트패턴 채우기(세그먼트 집합 교차곱; 같은 쌍은 1회만) ===
    // π+ × π−  (X=π+, Y=π−)
    for(int sx : seg_piPlus){
      for(int sy : seg_piMinus){
        hPiPlus_vs_PiMinus->Fill(sx, sy);
      }
    }
    // p × π−   (요청 축: X=p, Y=π−)
    for(int sx : seg_proton){
      for(int sy : seg_piMinus){
        hProton_vs_PiMinus->Fill(sx, sy);
      }
    }
  }

  auto pct=[](Long64_t a,Long64_t b){ return (b>0)?100.0*double(a)/double(b):0.0; };

  std::cout<<std::fixed<<std::setprecision(2);
  std::cout<<"========== SUMMARY (exclude 15-20) ==========\n";
  std::cout<<"Total events                          : "<<N_total<<"\n";
  std::cout<<"BH2 pass                              : "<<N_bh2<<" ("<<pct(N_bh2,N_total)<<" % of total)\n";
  std::cout<<"BH2 pass & HTOF>=1                    : "<<N_bh2_htofAny<<" ("<<pct(N_bh2_htofAny,N_bh2)<<" % of BH2)\n";
  std::cout<<"BH2 pass & HTOF>=1 & !{15..20}        : "<<N_cleanDenom<<" ("<<pct(N_cleanDenom,N_bh2)<<" % of BH2)\n";
  std::cout<<"\nMultiplicity>=2                       : "<<N_clean_Mge2<<" ("<<pct(N_clean_Mge2,N_cleanDenom)<<" %)\n";
  std::cout<<"Multiplicity==2                       : "<<N_clean_Meq2<<" ("<<pct(N_clean_Meq2,N_cleanDenom)<<" %)\n";
  std::cout<<"  - pair(2,3)                         : "<<C23<<" ("<<pct(C23,N_clean_Meq2)<<" % of M==2)\n";
  std::cout<<"  - pair(3,4)                         : "<<C34<<" ("<<pct(C34,N_clean_Meq2)<<" % of M==2)\n";
  std::cout<<"  - pair(4,5)                         : "<<C45<<" ("<<pct(C45,N_clean_Meq2)<<" % of M==2)\n";
  std::cout<<"  - pair(5,6)                         : "<<C56<<" ("<<pct(C56,N_clean_Meq2)<<" % of M==2)\n";
  std::cout<<"============================================\n";

  TCanvas* c=new TCanvas("cHTOF_SUM_excl15to20","BH2-gated, exclude 15-20",1700,950);
  c->Divide(3,2);
  c->cd(1); hSeg_unique->Draw("hist");
  c->cd(2); hPlaneLocal->Draw("colz");
  c->cd(3); hMult->Draw("hist");
  c->cd(4); hPiPlus_vs_PiMinus->Draw("colz");
  c->cd(5); hProton_vs_PiMinus->Draw("colz");
  c->cd(6); gPad->Clear(); TLatex tl; tl.SetNDC(); tl.SetTextSize(0.04);
           tl.DrawLatex(0.05,0.85,"HTOF 2D hit-patterns are 32#times32 (1#times1 cells).");
           tl.DrawLatex(0.05,0.75,"Counts per unique-segment pair per event.");
}

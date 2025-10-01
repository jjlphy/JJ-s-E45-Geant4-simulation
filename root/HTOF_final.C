// HTOF_Summary_BH2_Cuts_Exclude15to20.C
// - Exclude seg 15,16,17,18,19,20

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
static const double kL           = 337.0;  // [mm]
static const double kHTOFx       = 0.0;
static const double kHTOFy       = 12.0;
static const double kHTOFz       = 0.0;

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
  TFile* f=TFile::Open(filename);
  if(!f||f->IsZombie()){ printf("Cannot open %s\n",filename); return; }
  TTree* T=(TTree*)f->Get(treename);
  if(!T){ printf("Tree %s not found\n",treename); return; }

  std::vector<TParticle>* BH2=nullptr;
  std::vector<TParticle>* HTOF=nullptr;
  T->SetBranchAddress("BH2",&BH2);
  T->SetBranchAddress("HTOF",&HTOF);

  TH1D* hSeg_unique = new TH1D("hHTOF_seg_unique",
    "HTOF Hit Pattern (UNIQUE/event, exclude 15-20);Seg ID;Counts",32,-0.5,31.5);
  TH2D* hPlaneLocal = new TH2D("hHTOF_plane_local",
    "Plane vs Local (exclude 15-20);Plane;Local",8,-0.5,7.5,4,-0.5,3.5);
  TH1D* hMult = new TH1D("hHTOF_mult",
    "Multiplicity per event (exclude 15-20);# unique hit seg;Events",
    mult_max+1,-0.5,mult_max+0.5);

  Long64_t N_total=T->GetEntries();
  Long64_t N_bh2=0, N_bh2_htofAny=0, N_cleanDenom=0;
  Long64_t N_clean_Mge2=0, N_clean_Meq2=0;
  Long64_t C23=0,C34=0,C45=0,C56=0;

  for(Long64_t ie=0; ie<N_total; ++ie){
    T->GetEntry(ie);
    if(!BH2||!HTOF) continue;

    if(BH2->empty()) continue;
    N_bh2++;

    std::set<int> segset;
    for(const auto& p : *HTOF){
      int plane,local,seg;
      if(!mapHitToSeg(p.Vx(),p.Vy(),p.Vz(),
                      L_sign,plane_offset,flip_local,
                      plane,local,seg,R_tol_mm)) continue;
      if(seg>=0 && seg<kNSegHTOF) segset.insert(seg);
    }

    if(segset.empty()) continue;
    N_bh2_htofAny++;

    // --- 확장된 배제 조건 ---
    bool hasExcl=false;
    for(int s : {15,16,17,18,19,20})
      if(segset.count(s)) { hasExcl=true; break; }
    if(hasExcl) continue;

    N_cleanDenom++;

    int M=0;
    for(int s: segset){
      hSeg_unique->Fill(s);
      hPlaneLocal->Fill(s/4, s%4);
      M++;
    }
    if(M>mult_max) M=mult_max;
    hMult->Fill(M);

    if(M>=2) N_clean_Mge2++;
    if(M==2){
      N_clean_Meq2++;
      if(segset.count(2)&&segset.count(3)) C23++;
      if(segset.count(3)&&segset.count(4)) C34++;
      if(segset.count(4)&&segset.count(5)) C45++;
      if(segset.count(5)&&segset.count(6)) C56++;
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

  TCanvas* c=new TCanvas("cHTOF_SUM_excl15to20","BH2-gated, exclude 15-20",1400,800);
  c->Divide(2,2);
  c->cd(1); hSeg_unique->Draw("hist");
  c->cd(2); hPlaneLocal->Draw("colz");
  c->cd(3); hMult->Draw("hist");
}

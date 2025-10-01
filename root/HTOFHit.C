// HTOF_HitPattern_ALIGNED.C
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
#include <TSystem.h>
#include <TInterpreter.h>

static const int    kNPlaneHTOF  = 8;
static const int    kSegPerPlane = 4;
static const int    kNSegHTOF    = kNPlaneHTOF * kSegPerPlane;
static const double kPitch       = 68.0;   // mm
static const double kL           = 337.0;  // mm
static const double kHTOFx       = 0.0;
static const double kHTOFy       = 12.0;
static const double kHTOFz       = 0.0;

struct _DICT_ {
  _DICT_(){ gSystem->Load("libPhysics");
            gInterpreter->GenerateDictionary("vector<TParticle>","TParticle.h;vector"); }
} _dict_;

inline void toLocalCenter(double x,double y,double z,double& lx,double& ly,double& lz){
  lx=x-kHTOFx; ly=y-kHTOFy; lz=z-kHTOFz;
}

// 각 plane의 로컬 축 (ex: 가로, ez: 반경)
inline void planeAxes(int i, double& exx, double& exz, double& ezx, double& ezz){
  const double th = i * 2.0 * TMath::Pi() / kNPlaneHTOF; // i*45°
  exx =  TMath::Cos(th);     exz = -TMath::Sin(th);  // ex_i = (cos, -sin)
  ezx =  TMath::Sin(th);     ezz =  TMath::Cos(th);  // ez_i = (sin,  cos)
}

// xloc→local(0..3)
inline int inferLocalFromX(double xloc, bool flip_local=false){
  // 중앙 오프셋 1.5: (-1.5,-0.5,+0.5,+1.5)*Pitch
  double idx_f = xloc/kPitch + 1.5;
  int j = TMath::Nint(idx_f);
  if(j<0) j=0; if(j>=kSegPerPlane) j=kSegPerPlane-1;
  if(flip_local) j = (kSegPerPlane-1)-j; // 0<->3, 1<->2
  return j;
}

/**
 * filename, treename: 입력 ROOT
 * R_tol_mm: 반경 필터 (L±R_tol_mm만 카운트)
 * L_sign:   plane 판정시 목표 반경 부호 (+1이면 zloc≈+L, -1이면 zloc≈-L)
 * plane_offset: plane 번호를 최종적으로 +offset (GUI 번호와 정렬용)
 * flip_local:   local(0↔3,1↔2) 뒤집기
 * mult_max:     Multiplicity 히스토그램 최대 bin (기본 7 → 0..7)
 */
void HTOF_HitPattern_ALIGNED(const char* filename="E45_with_HTOF.root",
                             const char* treename="g4hyptpc",
                             double R_tol_mm=200.0,
                             int    L_sign=-1,          // GUI가 plane0을 +L로 본다면 +1로!
                             int    plane_offset=0,     // GUI와 1칸/4칸 어긋나면 여기로 보정
                             bool   flip_local=false,
                             int    mult_max=7)
{
  TFile* f=TFile::Open(filename);
  if(!f||f->IsZombie()){ printf("Cannot open %s\n",filename); return; }
  TTree* T=(TTree*)f->Get(treename);
  if(!T){ printf("Tree %s not found\n",treename); return; }

  std::vector<TParticle>* hits=nullptr;
  if(!T->GetBranch("HTOF")){ printf("Branch 'HTOF' not found.\n"); return; }
  T->SetBranchAddress("HTOF",&hits);

  // 히스토그램: ALL_HITS / UNIQUE / Plane×Local / Multiplicity
  TH1D* hSeg_all    = new TH1D("hHTOF_seg_all","HTOF Hit Pattern (ALL HITS);Seg ID (0-31);Counts",32,-0.5,31.5);
  TH1D* hSeg_unique = new TH1D("hHTOF_seg_unique","HTOF Hit Pattern (UNIQUE / event);Seg ID (0-31);Counts",32,-0.5,31.5);
  TH2D* hPL         = new TH2D("hHTOF_plane_local","Plane vs Local;Plane;Local",8,-0.5,7.5,4,-0.5,3.5);
  TH1D* hMult       = new TH1D("hHTOF_mult","HTOF Multiplicity per event;# unique hit seg;Events",
                                mult_max+1, -0.5, mult_max+0.5);

  Long64_t nent=T->GetEntries();
  std::vector<char> seen(kNSegHTOF, 0);

  for(Long64_t ie=0; ie<nent; ++ie){
    T->GetEntry(ie);
    if(!hits) continue;
    std::fill(seen.begin(), seen.end(), 0);
    int nuniq=0;

    for(const auto& p : *hits){
      double lx,ly,lz; toLocalCenter(p.Vx(),p.Vy(),p.Vz(),lx,ly,lz);
      const double R = std::hypot(lx,lz);
      if(R_tol_mm>0 && std::abs(R - kL) > R_tol_mm) continue;

      // --- plane 고르기: zloc_i가 L_sign*kL 에 가장 가까운 i ---
      int best_i = -1; double best_abs = 1e99;
      double best_xloc = 0.0, best_zloc = 0.0;
      for(int i=0;i<kNPlaneHTOF;++i){
        double exx,exz,ezx,ezz; planeAxes(i,exx,exz,ezx,ezz);
        const double zloc = lx*ezx + lz*ezz;         // 반경 성분
        const double diff = std::abs(zloc - (L_sign * kL));
        if(diff < best_abs){
          best_abs = diff;
          best_i   = i;
          best_zloc= zloc;
          best_xloc= lx*exx + lz*exz;               // 가로 성분
        }
      }
      if(best_i<0) continue;

      int plane = (best_i + plane_offset) % kNPlaneHTOF;
      if(plane<0) plane += kNPlaneHTOF;

      int local = inferLocalFromX(best_xloc, flip_local);
      int seg   = plane*kSegPerPlane + local;

      // === ALL_HITS ===
      hSeg_all->Fill(seg);

      // === UNIQUE_PER_EVENT ===
      if(!seen[seg]){
        seen[seg]=1; nuniq++;
        hSeg_unique->Fill(seg);
        hPL->Fill(plane, local);
      }
    }
    if(nuniq>mult_max) nuniq = mult_max;
    hMult->Fill(nuniq);
  }

  // 그리기
  TCanvas* c=new TCanvas("cHTOF","HTOF Hits (Aligned)",1400,700);
  c->Divide(2,2);
  c->cd(1); hSeg_all->Draw("hist");
  c->cd(2); hSeg_unique->SetLineColor(kBlue+1); hSeg_unique->Draw("hist");
  c->cd(3); hPL->Draw("colz");
  c->cd(4); hMult->Draw("hist");

  std::cout<<"[ALIGN] L_sign="<<L_sign<<" (target zloc="<<L_sign*kL<<")"
           <<", plane_offset="<<plane_offset
           <<", flip_local="<<(flip_local?"true":"false")<<"\n";
}

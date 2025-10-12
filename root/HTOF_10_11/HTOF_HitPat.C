// HTOF_HitPat34.C  — Draw 34 physical-tile hit pattern (upper/lower split on plane=0, local=1,2)
//
// Usage:
//   root -l
//   .L HTOF_HitPat34.C+
//   HTOF_HitPat34("E45_VP5.root","g4hyptpc", /*Rtol_mm=*/25.0, /*L_sign=*/+1,/*plane_offset=*/0, /*flip_local=*/false, /*UNIQUE_PER_EVENT=*/true);

//                 /*plane_offset=*/0, /*flip_local=*/false, /*UNIQUE_PER_EVENT=*/true);

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TSystem.h"
#include "TInterpreter.h"
#include "TParticle.h"
#include "TString.h"
#include <vector>
#include <set>
#include <cmath>
#include <iostream>

// ======================== geometry constants (HTOF) =========================
static const int    kNPlaneHTOF  = 8;
static const int    kSegPerPlane = 4;             // logical bars per plane
static const double kPitch       = 68.0;          // [mm] tangential pitch
static const double kL           = 337.0;         // [mm] ring radius
static const double kHTOFx       = 0.0;           // HTOF center (global)
static const double kHTOFy       = 12.0;
static const double kHTOFz       = 0.0;

// ---------------- dictionary for vector<TParticle> --------------------------
struct _DICT_ {
  _DICT_(){
    gSystem->Load("libPhysics");
    gInterpreter->GenerateDictionary("vector<TParticle>","TParticle.h;vector");
  }
} _dict_;

// ---------------- small helpers (local axes & mapping) ----------------------
inline void toLocalCenter(double x,double y,double z,double& lx,double& ly,double& lz){
  lx = x - kHTOFx;  ly = y - kHTOFy;  lz = z - kHTOFz;
}
inline void planeAxes(int i, double& exx, double& exz, double& ezx, double& ezz){
  const double th = i * 2.0 * M_PI / kNPlaneHTOF; // i*45°
  exx =  std::cos(th);  exz = -std::sin(th);       // tangential axis
  ezx =  std::sin(th);  ezz =  std::cos(th);       // radial axis
}
inline int inferLocalFromX(double xloc, bool flip_local=false){
  double idx_f = xloc / kPitch + 1.5;     // centers near -1.5,-0.5,+0.5,+1.5
  int j = (int)std::lround(idx_f);        // 0..3
  if(j<0) j=0; if(j>=kSegPerPlane) j=kSegPerPlane-1;
  if(flip_local) j = (kSegPerPlane-1) - j;
  return j;
}

// Map a hit position to (plane, local, seg) using nearest-radial-plane method.
// Returns false if outside radial tolerance.
inline bool mapHitToSeg(double x,double y,double z,
                        int L_sign,int plane_offset,bool flip_local,
                        int& plane,int& local,int& seg,double R_tol_mm)
{
  double lx,ly,lz; toLocalCenter(x,y,z,lx,ly,lz);
  const double R = std::hypot(lx,lz);
  if(R_tol_mm>0 && std::abs(R - kL) > R_tol_mm) return false;

  int best_i=-1; double best_abs=1e99; double best_xloc=0.0;
  for(int i=0;i<kNPlaneHTOF;++i){
    double exx,exz,ezx,ezz; planeAxes(i,exx,exz,ezx,ezz);
    const double zloc = lx*ezx + lz*ezz;          // radial coordinate along plane i
    const double diff = std::abs(zloc - (L_sign*kL));
    if(diff < best_abs){
      best_abs = diff; best_i = i;
      best_xloc = lx*exx + lz*exz;                // tangential coord for local index
    }
  }
  if(best_i<0) return false;

  plane = (best_i + plane_offset) % kNPlaneHTOF;
  if(plane<0) plane += kNPlaneHTOF;

  local = inferLocalFromX(best_xloc, flip_local);
  seg   = plane*kSegPerPlane + local;             // 0..31 (logical)
  return true;
}

// --------- map (plane,local,ly) to physical-tile ID (0..33) = 34 tiles -------
// Only plane=0, local in {1,2} are split into U/L using ly>0 or <0.
inline int mapToTile34(int plane, int local, double ly){
  int seg = plane*4 + local; // logical seg (0..31)
  if(plane==0){
    if(local==0) return 0;                       // plane0-local0  -> tile 0
    if(local==3) return 5;                       // plane0-local3  -> tile 5
    if(local==1) return (ly>0.0 ? 1 : 2);        // split upper/lower
    if(local==2) return (ly>0.0 ? 3 : 4);        // split upper/lower
  }
  // planes 1..7 shift by +2 (since 0..5 already used)
  return seg + 2;                                // seg=4..31 -> tile=6..33
}

// =============================== main =======================================
void HTOF_HitPat34(const char* filename="E45_VP5.root",
                   const char* treename="g4hyptpc",
                   double R_tol_mm = 25.0,
                   int    L_sign    = +1,
                   int    plane_offset = 0,
                   bool   flip_local  = false,
                   bool   UNIQUE_PER_EVENT = true)
{
  // open & branches
  TFile* f=TFile::Open(filename);
  if(!f||f->IsZombie()){ std::cerr<<"Cannot open "<<filename<<"\n"; return; }
  TTree* T=(TTree*)f->Get(treename);
  if(!T){ std::cerr<<"Tree "<<treename<<" not found\n"; return; }

  std::vector<TParticle>* HTOF = nullptr;
  T->SetBranchAddress("HTOF",&HTOF);

  // histo: 34-bin physical tiles
  TH1D* h = new TH1D("hHTOF_tile34",
     "HTOF Physical Tile Hit Pattern;Tile ID (0..33);Counts", 34, -0.5, 33.5);

  const Long64_t N = T->GetEntries();
  for(Long64_t ie=0; ie<N; ++ie){
    T->GetEntry(ie);
    if(!HTOF) continue;

    std::set<int> filled_this_event; // for UNIQUE_PER_EVENT

    for(const auto& p : *HTOF){
      int plane, local, seg;
      if(!mapHitToSeg(p.Vx(), p.Vy(), p.Vz(),
                      L_sign, plane_offset, flip_local,
                      plane, local, seg, R_tol_mm)) continue;

      double lx,ly,lz; toLocalCenter(p.Vx(), p.Vy(), p.Vz(), lx, ly, lz);
      int tile = mapToTile34(plane, local, ly);

      if(UNIQUE_PER_EVENT){
        if(filled_this_event.insert(tile).second){
          h->Fill(tile);
        }
      }else{
        h->Fill(tile);
      }
    }
  }

  // draw
  TCanvas* c = new TCanvas("cHTOF34","HTOF 34-tile Hit Pattern",900,500);
  h->SetLineWidth(2);
  h->Draw("hist");
  c->Update();

  std::cout << "[INFO] Filled " << (UNIQUE_PER_EVENT?"unique ":"all ")
            << "hits into 34-tile pattern. Entries=" << h->GetEntries() << "\n";
}

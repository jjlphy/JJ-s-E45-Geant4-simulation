// HTOF_seg.C
// - 목적: "HTOF의 특정 세그먼트(기본 15-20)에 히트한 이벤트"만 모아,
//         그 이벤트들에서 BH2가 어느 세그먼트를 맞았는지 확인
// - 결과:
//    (1) BH2 HitPattern (1D)
//    (2) (HTOF seg in selection) vs (BH2 seg) 2D
//    (3) 터미널에 BH2 세그먼트별 카운트 출력 (이벤트당 유일 세그 기준)
//
// 사용 예:
//   root -l
//   .L HTOF_seg.C+
//   // 기본: HTOF 15-20 선택, BH2 세그 수 15(0..14)
//   HTOF_seg("E45_piplusn.root","g4hyptpc");
//   // 커스텀: HTOF 2,4,7-9 선택 / BH2 세그 수 16(0..15)
//   HTOF_seg("E45_piplusn.root","g4hyptpc", "2,4,7-9", 200.0, +1, 0, false, 16);

#include <TFile.h>
#include <TTree.h>
#include <TParticle.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TMath.h>
#include <TSystem.h>
#include <TInterpreter.h>

#include <vector>
#include <set>
#include <unordered_set>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <algorithm>

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
    if(diff < best_abs){ best_abs = diff; best_i = i; best_xloc = lx*exx + lz*exz; }
  }
  if(best_i<0) return false;

  int plane_raw = best_i;
  plane = (plane_raw + plane_offset) % kNPlaneHTOF;
  if(plane<0) plane += kNPlaneHTOF;

  local = inferLocalFromX(best_xloc, flip_local);
  seg   = plane*kSegPerPlane + local;
  return true;
}

// spec 파서: "15-20,22" / "all" / "none" / ""
static std::unordered_set<int> parse_spec(const char* spec, int lo, int hi){
  std::unordered_set<int> out;
  if(!spec) return out;
  std::string s(spec);
  s.erase(std::remove_if(s.begin(), s.end(), ::isspace), s.end());
  if(s.empty() || s=="none" || s=="NONE") return out;      // 빈 집합
  if(s=="all"  || s=="ALL" ) { for(int v=lo; v<=hi; ++v) out.insert(v); return out; }
  std::stringstream ss(s); std::string tok;
  while(std::getline(ss, tok, ',')){
    if(tok.empty()) continue;
    auto pos = tok.find('-');
    if(pos==std::string::npos){
      int v = std::stoi(tok); if(v>=lo && v<=hi) out.insert(v);
    }else{
      int a = std::stoi(tok.substr(0,pos));
      int b = std::stoi(tok.substr(pos+1));
      if(a>b) std::swap(a,b);
      a = std::max(a, lo); b = std::min(b, hi);
      for(int v=a; v<=b; ++v) out.insert(v);
    }
  }
  return out;
}

// 메인
void HTOF_seg(const char* filename="E45_piplusn.root",
              const char* treename="g4hyptpc",
              const char* htof_select_spec="15-20", // 이 세그먼트에 히트한 이벤트만 선별
              double R_tol_mm = 200.0,
              int    L_sign    = +1,
              int    plane_offset = 0,
              bool   flip_local  = false,
              int    BH2_NSEG    = 15)               // BH2 세그 수(0..BH2_NSEG-1)
{
  // 파일/트리
  TFile* f=TFile::Open(filename);
  if(!f||f->IsZombie()){ printf("Cannot open %s\n",filename); return; }
  TTree* T=(TTree*)f->Get(treename);
  if(!T){ printf("Tree %s not found\n",treename); return; }

  // 브랜치
  std::vector<TParticle>* BH2 = nullptr;
  std::vector<TParticle>* HTOF = nullptr;
  T->SetBranchAddress("BH2",&BH2);
  T->SetBranchAddress("HTOF",&HTOF);

  // 선택 세그먼트 집합(HTOF)
  const auto HSEL = parse_spec(htof_select_spec, 0, kNSegHTOF-1);

  // 히스토그램
  TH1D* hBH2 = new TH1D("hBH2_for_HTOFsel",
    Form("BH2 HitPattern (events with HTOF in {%s});BH2 SegID;Counts", htof_select_spec),
    BH2_NSEG, -0.5, BH2_NSEG-0.5);

  TH2D* hBH2_vs_HTOF = new TH2D("hBH2_vs_HTOFsel",
    Form("BH2 seg vs HTOF seg (only HTOF in {%s});HTOF seg;BH2 seg", htof_select_spec),
    kNSegHTOF, -0.5, kNSegHTOF-0.5, BH2_NSEG, -0.5, BH2_NSEG-0.5);

  std::vector<Long64_t> bh2_counts(BH2_NSEG, 0);

  Long64_t N_total = T->GetEntries();
  Long64_t N_sel   = 0;

  for(Long64_t ie=0; ie<N_total; ++ie){
    T->GetEntry(ie);
    if(!BH2||!HTOF) continue;

    // HTOF 유일 세그먼트 집합
    std::set<int> htof_segs_unique;
    for(const auto& p : *HTOF){
      int pl,loc,seg;
      if(!mapHitToSeg(p.Vx(),p.Vy(),p.Vz(), L_sign,plane_offset,flip_local, pl,loc,seg, R_tol_mm)) continue;
      if(seg<0 || seg>=kNSegHTOF) continue;
      htof_segs_unique.insert(seg);
    }
    if(htof_segs_unique.empty()) continue;

    // 선택 세그먼트와 교집합 검사
    bool pass=false;
    if(HSEL.empty()){
      // "none"이면 선택 없음 → 아무 이벤트도 통과시키지 않음 (혼동 방지용)
      pass=false;
    }else{
      for(int s: htof_segs_unique){ if(HSEL.count(s)){ pass=true; break; } }
    }
    if(!pass) continue;
    N_sel++;

    // BH2 유일 세그먼트
    std::set<int> bh2_segs_unique;
    for(const auto& p:*BH2){
      int sid = p.GetMother(1);
      if(0<=sid && sid<BH2_NSEG) bh2_segs_unique.insert(sid);
    }

    // 채우기
    for(int sid : bh2_segs_unique){
      hBH2->Fill(sid);
      bh2_counts[sid]++;
      // 어떤 HTOF 세그(선택집합 내)와 함께였는지 2D도 채움
      for(int hs : htof_segs_unique){
        if(HSEL.count(hs)) hBH2_vs_HTOF->Fill(hs, sid);
      }
    }
  }

  // 출력
  std::cout<<std::fixed<<std::setprecision(2);
  std::cout<<"========== HTOF_seg SUMMARY (HTOF select: {"<<htof_select_spec<<"}) ==========\n";
  std::cout<<"Total events                 : "<<N_total<<"\n";
  std::cout<<"Selected (has HTOF in set)   : "<<N_sel<<" ("<< (N_total?100.0*double(N_sel)/double(N_total):0.0) <<" %)\n";
  std::cout<<"BH2 hit counts (unique seg / event):\n";
  for(int i=0;i<BH2_NSEG;++i){
    std::cout<<"  BH2 seg "<<std::setw(2)<<i<<" : "<<bh2_counts[i]<<"\n";
  }
  std::cout<<"=================================================================\n";

  // 캔버스
  TCanvas* c = new TCanvas("cHTOF_seg","HTOF-selected → BH2 mapping",1400,600);
  c->Divide(2,1);
  c->cd(1); hBH2->SetFillColorAlpha(kAzure+1,0.35); hBH2->Draw("hist");
  c->cd(2); hBH2_vs_HTOF->Draw("colz");
}

// BH2_HTOF_SCH_Multi2.C  (v4)
// - EVENT-LEVEL EXCLUDE (코드 1 정책)
// - (1) BH2×SCH 2D에서 히트>=1인 (BH2,SCH) 좌표를 터미널/파일(BH2&&SCH_tight_cut.txt)로 출력/저장
// - (2) 2D를 복제해, 히트 있는 셀을 기준으로 양축 ±1 확장(dilation) 칸에 빨간 박스 오버레이
//       확장으로 선택된 모든 (BH2,SCH) 좌표를 파일(BH2&&SCH_fit_cut.txt)로 저장
/*root -l
.L BH2_HTOF_SCH_Multi2.C+
BH2_HTOF_SCH_Multi2("E45_beamthrough_10_6.root","g4hyptpc","15-20");
// 

BH2_HTOF_SCH_Multi2("E45_piplusn.root","g4hyptpc","15-20");
// 
*/

#include <TFile.h>
#include <TTree.h>
#include <TParticle.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TSystem.h>
#include <TInterpreter.h>
#include <TLatex.h>
#include <TBox.h>
#include <TStyle.h>
#include <TPaletteAxis.h>

#include <vector>
#include <set>
#include <unordered_set>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <cmath>

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
  exx =  std::cos(th);  exz = -std::sin(th);
  ezx =  std::sin(th);  ezz =  std::cos(th);
}
inline int inferLocalFromX(double xloc, bool flip_local=false){
  double idx_f = xloc / kPitch + 1.5;
  int j = TMath::Nint(idx_f);
  if(j<0) j=0; if(j>=kSegPerPlane) j=kSegPerPlane-1;
  if(flip_local) j = (kSegPerPlane-1)-j;
  return j;
}
inline bool mapHitToHTOFseg(double x,double y,double z,
                            int L_sign,int plane_offset,bool flip_local,
                            int& seg,double R_tol_mm)
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

  int plane = (best_i + plane_offset) % kNPlaneHTOF;
  if(plane<0) plane += kNPlaneHTOF;

  int local = inferLocalFromX(best_xloc, flip_local);
  seg   = plane*kSegPerPlane + local;
  return (seg>=0 && seg<kNSegHTOF);
}

static std::unordered_set<int> parse_spec(const char* spec, int lo, int hi){
  std::unordered_set<int> out;
  if(!spec) return out;
  std::string s(spec);
  s.erase(std::remove_if(s.begin(), s.end(), ::isspace), s.end());
  if(s.empty() || s=="none" || s=="NONE") return out; // 빈 집합 = 배제 없음
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

static inline void unique_seg_from_branch(const std::vector<TParticle>* v,
                                          int nmax, double ecutMeV,
                                          std::set<int>& out)
{
  out.clear();
  if(!v) return;
  for(const auto& p : *v){
    if(ecutMeV>0.0 && p.GetWeight()<ecutMeV) continue;
    int seg = p.GetMother(1);
    if(0<=seg && seg<nmax) out.insert(seg);
  }
}

void BH2_HTOF_SCH_Multi2(const char* filename="E45_piplusn.root",
                         const char* treename="g4hyptpc",
                         const char* htof_exclude_spec="15-20",
                         double ecutBH2MeV=0.0,
                         double ecutHTOFMeV=0.0,
                         double ecutSCHMeV=0.0,
                         int    BH2_NSEG=15,
                         int    SCH_NSEG=64,
                         // HTOF 매핑 옵션
                         double R_tol_mm = 200.0,
                         int    L_sign    = +1,
                         int    plane_offset = 0,
                         bool   flip_local  = false)
{
  gStyle->SetOptStat(0);

  // 파일/트리
  TFile* f=TFile::Open(filename);
  if(!f||f->IsZombie()){ printf("Cannot open %s\n",filename); return; }
  TTree* T=(TTree*)f->Get(treename);
  if(!T){ printf("Tree %s not found\n",treename); return; }

  // 브랜치
  std::vector<TParticle> *BH2=nullptr, *HTOF=nullptr, *SCH=nullptr;
  T->SetBranchAddress("BH2",&BH2);
  T->SetBranchAddress("HTOF",&HTOF);
  if(T->GetBranch("SCH")) T->SetBranchAddress("SCH",&SCH);
  else { printf("Branch 'SCH' not found.\n"); return; }

  const auto HTOF_EXCL = parse_spec(htof_exclude_spec, 0, kNSegHTOF-1);

  // 히스토그램
  TH1D* hSCH = new TH1D("hSCH_Multi2",
    Form("SCH HitPattern (BH2 hit & HTOF>=1 & !{excl}, Multi#geq2);SCH SegID;Counts"),
    SCH_NSEG, -0.5, SCH_NSEG-0.5);

  TH2D* hBH2_vs_SCH = new TH2D("hBH2_vs_SCH",
    "BH2 seg vs SCH seg (denominator events);BH2 seg;SCH seg",
    BH2_NSEG, -0.5, BH2_NSEG-0.5,
    SCH_NSEG, -0.5, SCH_NSEG-0.5);

  // (2) 표시용: 동일 binning의 복제
  TH2D* hBH2_vs_SCH_copy = (TH2D*)hBH2_vs_SCH->Clone("hBH2_vs_SCH_copy");
  hBH2_vs_SCH_copy->SetTitle("BH2 seg vs SCH seg (with \\pm1 dilation overlay);BH2 seg;SCH seg");

  std::vector<Long64_t> sch_counts(SCH_NSEG, 0);

  // 카운터
  const Long64_t N_total = T->GetEntries();
  Long64_t N_bh2_any = 0;
  Long64_t N_pass1   = 0; // BH2 hit & HTOF>=1 & !{excl}
  Long64_t N_denom   = 0; // + multiplicity>=2 (htof_all 기준)
  Long64_t N_sch_any = 0; // 분모에서 SCH>=1
  Long64_t N_sch_none= 0; // 분모에서 SCH==0

  // 루프
  for(Long64_t ie=0; ie<N_total; ++ie){
    T->GetEntry(ie);
    if(!BH2||!HTOF||!SCH) continue;

    // BH2 유일 세그
    std::set<int> bh2_segs;
    unique_seg_from_branch(BH2, BH2_NSEG, ecutBH2MeV, bh2_segs);
    if(bh2_segs.empty()) continue;
    N_bh2_any++;

    // HTOF 전체 유일 세그
    std::set<int> htof_all;
    for(const auto& p : *HTOF){
      if(ecutHTOFMeV>0.0 && p.GetWeight()<ecutHTOFMeV) continue;
      int seg=-1;
      if(!mapHitToHTOFseg(p.Vx(),p.Vy(),p.Vz(), L_sign,plane_offset,flip_local, seg, R_tol_mm)) continue;
      if(seg<0 || seg>=kNSegHTOF) continue;
      htof_all.insert(seg);
    }
    if(htof_all.empty()) continue;

    // 이벤트-레벨 배제
    bool has_excl=false;
    for(int s : htof_all){ if(HTOF_EXCL.count(s)){ has_excl=true; break; } }
    if(has_excl) continue;

    N_pass1++;

    // 분모: HTOF mult>=2 (배제 없이)
    if((int)htof_all.size() >= 2){
      N_denom++;

      // SCH 유일 세그
      std::set<int> sch_unique;
      unique_seg_from_branch(SCH, SCH_NSEG, ecutSCHMeV, sch_unique);

      if(sch_unique.empty()){
        N_sch_none++;
      }else{
        N_sch_any++;
        for(int s : sch_unique){ hSCH->Fill(s); sch_counts[s]++; }
      }

      // BH2×SCH 2D (교차곱)
      for(int hb : bh2_segs){
        for(int sc : sch_unique){
          hBH2_vs_SCH->Fill(hb, sc);
          hBH2_vs_SCH_copy->Fill(hb, sc);
        }
      }
    }
  }

  // ===== (1) BH2×SCH에서 히트>=1인 좌표를 파일/터미널로 출력 =====
  std::ofstream tightOut("BH2&&SCH_tight_cut.txt");
  if(!tightOut.is_open()){
    std::cerr << "[ERR] cannot open output: BH2&&SCH_tight_cut.txt\n";
  }
  std::cout << "\n(1) Non-zero cells in BH2×SCH (tight):\n";
  int n_tight = 0;
  for(int bx=1; bx<=hBH2_vs_SCH->GetNbinsX(); ++bx){
    for(int by=1; by<=hBH2_vs_SCH->GetNbinsY(); ++by){
      if(hBH2_vs_SCH->GetBinContent(bx,by) > 0){
        int h = bx-1, s = by-1;
        ++n_tight;
        std::cout << "  (BH2,SCH)=(" << h << ", " << s << ")\n";
        if(tightOut.is_open()) tightOut << "(BH2,SCH)=(" << h << ", " << s << ")\n";
      }
    }
  }
  if(tightOut.is_open()) tightOut.close();
  std::cout << "  -> total " << n_tight << " pairs\n";
  std::cout << "  saved to: 'BH2&&SCH_tight_cut.txt'\n";

  // ===== (2) ±1 확장(dilation) 마스크 생성 & 빨간 박스 오버레이 & 파일 저장 =====
  // 마스크: false로 초기화
  std::vector<char> dilated(BH2_NSEG * SCH_NSEG, 0);
  auto idx = [SCH_NSEG](int h, int s){ return h*SCH_NSEG + s; };

  for(int bx=1; bx<=hBH2_vs_SCH->GetNbinsX(); ++bx){
    for(int by=1; by<=hBH2_vs_SCH->GetNbinsY(); ++by){
      if(hBH2_vs_SCH->GetBinContent(bx,by) > 0){
        int h = bx-1, s = by-1;
        // 양축 ±1 확장 (경계 체크)
        for(int dh=-1; dh<=1; ++dh){
          for(int ds=-1; ds<=1; ++ds){
            int hh = h + dh, ss = s + ds;
            if(0<=hh && hh<BH2_NSEG && 0<=ss && ss<SCH_NSEG){
              dilated[idx(hh,ss)] = 1;
            }
          }
        }
      }
    }
  }

  // 파일 저장 (확장된 좌표 전부)
  std::ofstream fitOut("BH2&&SCH_fit_cut.txt");
  if(!fitOut.is_open()){
    std::cerr << "[ERR] cannot open output: BH2&&SCH_fit_cut.txt\n";
  }
  int n_fit = 0;
  for(int h=0; h<BH2_NSEG; ++h){
    for(int s=0; s<SCH_NSEG; ++s){
      if(dilated[idx(h,s)]){
        ++n_fit;
        if(fitOut.is_open()) fitOut << "(BH2,SCH)=(" << h << ", " << s << ")\n";
      }
    }
  }
  if(fitOut.is_open()) fitOut.close();
  std::cout << "\n(2) Dilation(±1 in both axes) selected total " << n_fit << " pairs\n";
  std::cout << "  saved to: 'BH2&&SCH_fit_cut.txt'\n";

  // 캔버스: (왼) SCH 1D, (오) BH2×SCH 2D + 오버레이
  TCanvas* c = new TCanvas("cSCH_multi2","SCH & BH2xSCH (denominator events)",1200,500);
  c->Divide(2,1);

  // Left: SCH 1D
  c->cd(1);
  hSCH->SetFillColorAlpha(kOrange+7,0.35);
  hSCH->Draw("hist");

  // Right: 2D base + red overlay
  c->cd(2);
  gPad->SetRightMargin(0.14);
  hBH2_vs_SCH_copy->Draw("colz");
  gPad->Update();

  // 빨간 박스 오버레이
  for(int h=0; h<BH2_NSEG; ++h){
    for(int s=0; s<SCH_NSEG; ++s){
      if(!dilated[idx(h,s)]) continue;
      int bx = hBH2_vs_SCH_copy->GetXaxis()->FindBin(h);
      int by = hBH2_vs_SCH_copy->GetYaxis()->FindBin(s);
      double x1 = hBH2_vs_SCH_copy->GetXaxis()->GetBinLowEdge(bx);
      double x2 = hBH2_vs_SCH_copy->GetXaxis()->GetBinUpEdge  (bx);
      double y1 = hBH2_vs_SCH_copy->GetYaxis()->GetBinLowEdge(by);
      double y2 = hBH2_vs_SCH_copy->GetYaxis()->GetBinUpEdge  (by);
      TBox *box = new TBox(x1,y1,x2,y2);
      box->SetFillStyle(0);
      box->SetLineColor(kRed);
      box->SetLineWidth(2);
      box->Draw("same");
    }
  }
  gPad->RedrawAxis();

  // 요약 출력
  auto pct=[](Long64_t a,Long64_t b){ return (b>0)?100.0*double(a)/double(b):0.0; };
  std::cout<<std::fixed<<std::setprecision(2);
  std::cout<<"\n==================== SUMMARY (EVENT-LEVEL EXCLUDE) ====================\n";
  std::cout<<"Input file                         : "<<filename<<"\n";
  std::cout<<"Total events                       : "<<N_total<<"\n";
  std::cout<<"BH2 hit (any seg)                  : "<<N_bh2_any<<"\n";
  std::cout<<"BH2 hit & HTOF>=1 & !{excl}        : "<<N_pass1<<" ("<<pct(N_pass1,N_bh2_any)<<" %)\n";
  std::cout<<"Denominator (HTOF mult>=2, no excl): "<<N_denom<<" ("<<pct(N_denom,N_pass1)<<" %)\n";
  std::cout<<"  SCH>=1 / SCH==0 (denom)          : "<<N_sch_any<<" / "<<N_sch_none
           <<"  ("<<pct(N_sch_any,N_denom)<<" %  / "<<pct(N_sch_none,N_denom)<<" %)\n";
  std::cout<<"-----------------------------------------------------------------------\n";
  std::cout<<"(1) Non-zero pairs (tight)         : "<<n_tight<<"  -> BH2&&SCH_tight_cut.txt\n";
  std::cout<<"(2) Dilation(±1) pairs (fit)       : "<<n_fit  <<"  -> BH2&&SCH_fit_cut.txt\n";
  std::cout<<"=======================================================================\n";
}

// HTOF_trigger_study_2pi.C  (for forced-2pi sim)
// Selection pipeline:
//   0) load tree
//   1) require tgt_touch_flag == 1  (exclude all non-target events)
//   2) require BH2 in [bh2_lo, bh2_hi]  (denominator)
//   3) Trig1 = (HTOF mult >= 2)
//   4) BeamVeto = (20&21) OR (22&23)
//   5) Trig2 = Trig1 & (!BeamVeto)
//
// Policies:
//   - HTOF tile ID = TParticle::StatusCode()
//   - BH2 seg ID   = coordinate → segment (same mapper as before)
//   - edep source  = *_edep branch if present, else TParticle::Weight()
//   - "valid"      = Sum(edep per ID) >= threshold (no fallback to 'any hit')
// Denominator for % : #events that pass (tgt_touch_flag==1) AND (BH2 in range)
//
// Usage example:
//   root -l
//   .L HTOF_trigger_study_2pi.C+
//   HTOF_trigger_study_2pi("E45_with_SCH.root","g4hyptpc",4,10,0.1,2.0,5.0,10.0,true,"BH2_4_10");

#include "TFile.h"
#include "TTree.h"
#include "TParticle.h"
#include "TSystem.h"
#include "TInterpreter.h"
#include "TH1I.h"
#include "TCanvas.h"
#include "TString.h"
#include <vector>
#include <map>
#include <set>
#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>

namespace HTS2PI {

// geometry (E72-like)
static const int    kNBH2Seg   = 15;     // 0..14
static const double kBH2_x0    = -10.0;  // [mm]
static const double kBH2_pitch = 14.0;   // [mm]
static const double kBH2_total = kNBH2Seg * kBH2_pitch;
static const int    kNHTOF     = 34;     // 0..33

struct DictGuard {
  DictGuard(){ gSystem->Load("libPhysics");
               gInterpreter->GenerateDictionary("vector<TParticle>","TParticle.h;vector");}
};
static DictGuard _dg;

static inline double MIPThrMeV(double frac,double dEdx,double thk_mm){
  return frac * dEdx * (thk_mm*0.1); // mm→cm
}

// x(world) → BH2 seg id
static int MapBH2_WorldToSeg(double x, double, double){
  const double xloc = x - kBH2_x0;
  const double xmin = -0.5*kBH2_total, xmax = 0.5*kBH2_total;
  if (xloc < xmin || xloc >= xmax) return -1;
  int idx = int(std::floor((xloc - xmin)/kBH2_pitch));
  if(idx<0) idx=0; if(idx>=kNBH2Seg) idx=kNBH2Seg-1;
  return idx;
}

static inline std::string JoinTuple(const std::set<int>& S){
  std::string out="(";
  bool first=true; for(int v:S){ if(!first) out+=','; first=false; out+=std::to_string(v); }
  out+=')'; return out;
}

} // ns
void HTOF_trigger_study_2pi(const char* filename,
                            const char* treename="g4hyptpc",
                            int bh2_lo=4, int bh2_hi=10,
                            double mipFrac=0.1,
                            double mipMeVperCm=2.0,
                            double BH2_thickness_mm=5.0,
                            double HTOF_thickness_mm=10.0,
                            bool save=true,
                            const char* tag="BH2_4_10")
{
  using namespace HTS2PI;

  // (생략) --- 파일/브랜치 오픈 및 안내 출력, 스칼라 정의 동일 ---

  auto pct = [](long long a, long long b)->double{ return (b>0)? (100.0*(double)a/(double)b):0.0; };

  // counters
  Long64_t N_total=0;
  Long64_t N_tgt=0;        // tgt_touch_flag==1
  Long64_t N_bh2In=0;      // tgt_touch==1 AND BH2 in-range
  Long64_t N_trig1=0;      // above AND HTOF mult>=2
  Long64_t N_trig2=0;      // above AND !BeamVeto

  // NEW: BeamVeto counters
  Long64_t N_beamVeto_any   = 0; // (tgt&bh2In)에서 BeamVeto 패턴 만족한 모든 이벤트
  Long64_t N_beamVeto_inTrig1 = 0; // Trig1 내부에서 BeamVeto 패턴 만족

  std::map<std::string, Long64_t> comboCounts;

  const Long64_t N=T->GetEntries();
  for(Long64_t ie=0; ie<N; ++ie){
    T->GetEntry(ie); N_total++;

    // 1) require target touch
    if(tgt_touch_flag!=1) continue;
    N_tgt++;

    // --- BH2 유효 세그먼트 집계 (동일) ---
    std::map<int,double> bh2E;
    if(BH2){
      for(size_t i=0;i<BH2->size();++i){
        const TParticle& p=BH2->at(i);
        int sid = MapBH2_WorldToSeg(p.Vx(),p.Vy(),p.Vz());
        if(0<=sid && sid<kNBH2Seg){
          const double ed = (hasBH2edep && BH2_edep && i<BH2_edep->size()) ? BH2_edep->at(i) : p.GetWeight();
          bh2E[sid]+=ed;
        }
      }
    }
    std::set<int> bh2Valid;
    for(auto& kv:bh2E) if(kv.second>=thrBH2) bh2Valid.insert(kv.first);

    // 2) require BH2 in [lo,hi]
    bool bh2_in_range=false;
    for(int s:bh2Valid){ if(bh2_lo<=s && s<=bh2_hi){ bh2_in_range=true; break; } }
    if(!bh2_in_range) continue;
    N_bh2In++;

    // --- HTOF 유효 타일 집계 (동일) ---
    std::map<int,double> htofE;
    if(HTOF){
      for(size_t i=0;i<HTOF->size();++i){
        const TParticle& p=HTOF->at(i);
        int tid = p.GetStatusCode();
        if(0<=tid && tid<kNHTOF){
          const double ed = (hasHTOFedep && HTOF_edep && i<HTOF_edep->size()) ? HTOF_edep->at(i) : p.GetWeight();
          htofE[tid]+=ed;
        }
      }
    }
    std::set<int> htofValid;
    for(auto& kv:htofE) if(kv.second>=thrHTOF) htofValid.insert(kv.first);
    const int mult = (int)htofValid.size();

    // === NEW: BeamVeto 패턴을 Trig1 여부와 무관하게 먼저 판정 ===
    const bool has20 = htofValid.count(20);
    const bool has21 = htofValid.count(21);
    const bool has22 = htofValid.count(22);
    const bool has23 = htofValid.count(23);
    const bool beamVeto = (has20 && has21) || (has22 && has23);

    // (A) (tgt&bh2In) 모수에서의 BeamVeto 건수
    if(beamVeto) N_beamVeto_any++;

    // 3) Trig1: multiplicity >= 2
    if(mult>=2){
      N_trig1++;

      // (B) Trig1 내부에서의 BeamVeto 건수
      if(beamVeto) N_beamVeto_inTrig1++;

      // 5) Trig2: Veto 통과(=survive)
      if(!beamVeto){
        N_trig2++;
        comboCounts[ JoinTuple(htofValid) ]++;
      }
    }
  }

  // ---- summary ----
  std::cout<<"\n==== Trigger summary (2pi sim) | BH2 "<<bh2_lo<<"-"<<bh2_hi<<" ====\n";
  std::cout<<"Total events                         : "<<N_total<<"\n";
  std::cout<<"Target-touched (tgt_touch_flag==1)   : "<<N_tgt
           <<"  ("<<pct(N_tgt,N_total)<<" % of Total)\n";
  std::cout<<"BH2 in-range after tgt-touch (denom) : "<<N_bh2In
           <<"  ("<<pct(N_bh2In,N_tgt)<<" % of tgt-touched)\n";
  std::cout<<"Trig1 (m>=2)                         : "<<N_trig1
           <<"  ("<<pct(N_trig1,N_bh2In)<<" % of BH2-range)\n";
  std::cout<<"BeamVeto (any, after tgt&BH2)        : "<<N_beamVeto_any
           <<"  ("<<pct(N_beamVeto_any,N_bh2In)<<" % of BH2-range)\n";
  std::cout<<"BeamVeto (inside Trig1)              : "<<N_beamVeto_inTrig1
           <<"  ("<<pct(N_beamVeto_inTrig1,N_trig1)<<" % of Trig1)\n";
  std::cout<<"Trig2 (!BeamVeto)                    : "<<N_trig2
           <<"  ("<<pct(N_trig2,N_bh2In)<<" % of BH2-range)\n";

  // ---- key ratios ----
  const double overkill      = (N_trig1>0) ? (100.0* (double)N_beamVeto_inTrig1 / (double)N_trig1) : 0.0;  // %
  const double survivor     = (N_trig1>0) ? (100.0* (double)N_trig2 / (double)N_trig1) : 0.0;              // %
  const double beamlikeRate = (N_bh2In>0)? (100.0* (double)N_beamVeto_any   / (double)N_bh2In) : 0.0;      // %

  std::cout<<"\n---- Derived ratios ----\n";
  std::cout<<"Overkill (BeamVeto/Trig1)            : "<<overkill<<" %\n";
  std::cout<<"Survivor (Trig2/Trig1)               : "<<survivor<<" %  (check: Overkill ≈ 100 - Survivor)\n";
  std::cout<<"Beam-like pattern rate (BeamVeto/Beam): "<<beamlikeRate<<" %  (Beam≡BH2 in-range)\n";



  // ---- histogram of Trig2 HTOF combos ----
  const int nComb = (int)comboCounts.size();
  TH1I* hComb = new TH1I("hHTOFComb_trig2_2pi",
                         Form("HTOF combinations (Trig2 survivors) | BH2 %d-%d;Combo;Events",
                              bh2_lo,bh2_hi),
                         std::max(1,nComb), 0.5, std::max(1,nComb)+0.5);
  hComb->SetDirectory(nullptr);

  int bin=1;
  for(const auto& kv:comboCounts){
    hComb->GetXaxis()->SetBinLabel(bin, kv.first.c_str());
    hComb->SetBinContent(bin, (double)kv.second);
    ++bin;
  }

  TCanvas* c=new TCanvas("cTrig2Comb_2pi","Trig2 HTOF combinations (2pi sim)",
                         std::max(900, 18*std::max(1,nComb)), 600);
  hComb->LabelsOption("v","X");
  hComb->Draw("hist");

  if(save){
    TString t(tag);
    c->SaveAs("Trig2_HTOFcombo_2pi_"+t+".png");
  }
}

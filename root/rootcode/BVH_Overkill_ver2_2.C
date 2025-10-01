// BVH_Overkill_ver2_2.C  (with per-BH2 UxD maps, palette shown, 0~50k option, logZ option)
// Axis names & palette title are aligned to BVH_3D_ver2.C: X=#BVH_U Seg, Y=#BVH_D Seg, Z=Counts

#include "Rtypes.h"
#include "TFile.h"
#include "TTree.h"
#include "TParticle.h"
#include "TSystem.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TBox.h"
#include "TLine.h"
#include "TPaletteAxis.h"
#include "TGaxis.h"
#include "TString.h"
#include "TLatex.h"

#include <vector>
#include <tuple>
#include <unordered_set>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <iostream>
#include <iomanip>

// -------- Detector Segmentation (Ver2) --------
static const int N_BH2   = 15; // 0..14
static const int N_BVHU  = 22; // 0..21
static const int N_BVHD  = 32; // 0..31

// ===== Helper: branch auto-bind =====
static bool bind_branch(TTree* tr, const char* preferred,
                        const std::vector<const char*>& fallbacks,
                        std::vector<TParticle>*& ptr) {
  if (tr->GetBranch(preferred)) { tr->SetBranchAddress(preferred, &ptr); return true; }
  for (auto nm : fallbacks) {
    if (tr->GetBranch(nm)) { tr->SetBranchAddress(nm, &ptr);
      Warning("Overkill","using fallback branch '%s' for '%s'", nm, preferred);
      return true;
    }
  }
  Error("Overkill","branch '%s' not found (no fallback matched)", preferred);
  return false;
}

// ===== Helper: unique segment indices passing e-dep cut =====
static inline void get_unique_hits(const std::vector<TParticle>* v,
                                   int nmax, double cutMeV,
                                   std::vector<int>& out)
{
  out.clear();
  if(!v) return;
  std::unordered_set<int> s;
  s.reserve(8);
  for(const auto& p : *v){
    const double edepMeV = p.GetWeight();
    if(edepMeV < cutMeV) continue;
    int seg = p.GetMother(1);
    if(0<=seg && seg<nmax) s.insert(seg);
  }
  out.assign(s.begin(), s.end());
  std::sort(out.begin(), out.end());
}

static inline int key3(int h,int u,int d){ return (h<<16) | (u<<8) | d; }

// ===== Parser: overlay_triplets_ver2.txt → (h,u,d) set =====
static bool load_overlay_triplets(const char* path,
                                  std::unordered_set<int>& overlay_set,
                                  std::vector<std::tuple<int,int,int>>& overlay_list)
{
  overlay_set.clear();
  overlay_list.clear();

  std::ifstream fin(path);
  if(!fin.is_open()){
    Error("Overkill","Cannot open overlay list file: %s", path);
    return false;
  }

  auto trim = [](std::string& s){
    size_t a = s.find_first_not_of(" \t\r\n");
    size_t b = s.find_last_not_of(" \t\r\n");
    if(a==std::string::npos) { s.clear(); return; }
    s = s.substr(a, b-a+1);
  };

  std::string line;
  int ln = 0;
  while(std::getline(fin, line)){
    ln++;
    trim(line);
    if(line.empty()) continue;
    for(char& c : line){
      if(c=='{' || c=='}' || c=='[' || c==']' || c=='(' || c==')' || c==',') c=' ';
    }
    int h,u,d; std::istringstream iss(line);
    if(!(iss >> h >> u >> d)){
      Warning("Overkill","Cannot parse line %d: '%s' (skipped)", ln, line.c_str());
      continue;
    }
    if(!(0<=h && h<N_BH2 && 0<=u && u<N_BVHU && 0<=d && d<N_BVHD)){
      Warning("Overkill","Out-of-range triplet at line %d: {%d,%d,%d} (skipped)", ln, h,u,d);
      continue;
    }
    int k = key3(h,u,d);
    if(overlay_set.insert(k).second) overlay_list.emplace_back(h,u,d);
  }
  fin.close();

  std::sort(overlay_list.begin(), overlay_list.end(),
            [](const std::tuple<int,int,int>& a, const std::tuple<int,int,int>& b){
              if (std::get<0>(a) != std::get<0>(b)) return std::get<0>(a) < std::get<0>(b);
              if (std::get<1>(a) != std::get<1>(b)) return std::get<1>(a) < std::get<1>(b);
              return std::get<2>(a) < std::get<2>(b);
            });
  Info("Overkill","Loaded %zu overlay triplets from %s", overlay_list.size(), path);
  return true;
}

// ================== Main Overkill + Maps ==================
void ComputeOverkillFromOverlay(const char* data_root = "../E45_2pi_Ver2.root",
                                const char* overlay_txt = "overlay_triplets_ver2.txt",
                                double ecutBH2MeV   = 0.10,
                                double ecutUMeV     = 0.04,
                                double ecutDMeV     = 0.04,
                                bool   draw_maps    = true,
                                bool   save_png     = true,
                                bool   use_fixed_z  = true,   // 0~50k 고정 여부
                                bool   use_logz     = false)  // 로그 스케일 사용 여부
{
  gSystem->Load("libPhysics");
  gInterpreter->GenerateDictionary("vector<TParticle>", "TParticle.h;vector");

  // overlay triplets
  std::unordered_set<int> overlay_set;
  std::vector<std::tuple<int,int,int>> overlay_list;
  if(!load_overlay_triplets(overlay_txt, overlay_set, overlay_list)){
    Error("Overkill","Failed to load overlay triplets; abort.");
    return;
  }

  // open data
  TFile* f = TFile::Open(data_root, "READ");
  if(!f || f->IsZombie()){ Error("Overkill","Cannot open file: %s", data_root); return; }
  TTree* tr = (TTree*)f->Get("g4hyptpc");
  if(!tr){ Error("Overkill","Cannot find TTree 'g4hyptpc'"); f->Close(); return; }

  std::vector<TParticle> *BH2=nullptr, *BVHU=nullptr, *BVHDv=nullptr;
  bool ok = true;
  ok &= bind_branch(tr, "BH2",    {},              BH2);
  ok &= bind_branch(tr, "BVH_U",  {"BVH"},         BVHU);
  ok &= bind_branch(tr, "BVH_D",  {"BVH2","VD"},   BVHDv);
  if(!ok){ f->Close(); return; }

  // per-BH2 U×D histograms  (axis titles aligned to BVH_3D_ver2.C)
  std::vector<TH2F*> h_ud(N_BH2, nullptr);
  if (draw_maps) {
    for (int h = 0; h < N_BH2; ++h) {
      TString name  = TString::Format("h_ud_bh2_%d", h);
      TString title = TString::Format("BH2=%d;#BVH1;#BVH2", h);
      h_ud[h] = new TH2F(name, title, N_BVHU, -0.5, N_BVHU-0.5, N_BVHD, -0.5, N_BVHD-0.5);
    }
  }

  // counting + fill maps
  Long64_t denom_all = 0, num_all = 0;
  std::vector<Long64_t> denom_h(N_BH2, 0), num_h(N_BH2, 0);

  const Long64_t N = tr->GetEntries();
  std::vector<int> hitsH, hitsU, hitsD;
  hitsH.reserve(8); hitsU.reserve(8); hitsD.reserve(8);

  std::cout << "[Info] Computing overkill on " << N << " events..." << std::endl;
  for(Long64_t i=0;i<N;++i){
    if(i && (i%100000==0)) std::cout << "  " << i << " / " << N << "\r" << std::flush;
    tr->GetEntry(i);

    get_unique_hits(BH2,   N_BH2,   ecutBH2MeV, hitsH);
    get_unique_hits(BVHU,  N_BVHU,  ecutUMeV,   hitsU);
    get_unique_hits(BVHDv, N_BVHD,  ecutDMeV,   hitsD);

    // fill maps only when all three exist
    if (draw_maps && !hitsH.empty() && !hitsU.empty() && !hitsD.empty()){
      for (int h : hitsH){
        TH2F* H = h_ud[h]; if(!H) continue;
        for (int u : hitsU){
          for (int d : hitsD){
            H->Fill(u, d);
          }
        }
      }
    }

    if(hitsH.empty() || hitsU.empty()) continue; // denom = BH2 && BVH_U
    denom_all++;

    bool global_veto = false;
    for(int h : hitsH){
      bool veto_h = false;
      denom_h[h]++;
      if(!hitsD.empty()){
        for(int u : hitsU){
          for(int d : hitsD){
            if(overlay_set.count(key3(h,u,d))){ veto_h = true; global_veto = true; break; }
          }
          if(veto_h) break;
        }
      }
      if(veto_h) num_h[h]++;
    }
    if(global_veto) num_all++;
  }
  std::cout << "\n[Info] Done." << std::endl;

  // print results
  std::cout << "\n========== Overkill Summary (Ver2) ==========\n";
  std::cout << std::fixed << std::setprecision(3);
  std::cout << "Data file             : " << data_root << "\n";
  std::cout << "Overlay list          : " << overlay_txt << " (" << overlay_list.size() << " cells)\n";
  std::cout << "Energy cuts (MeV)     : BH2=" << ecutBH2MeV
            << ", BVH_U=" << ecutUMeV << ", BVH_D=" << ecutDMeV << "\n";
  std::cout << "Segments (U,D)        : " << N_BVHU << ", " << N_BVHD << "\n";
  std::cout << "Denominator (BH2&&U)  : " << denom_all << "\n";
  if(denom_all>0){
    double R = (double)num_all / (double)denom_all * 100.0;
    std::cout << "Overkill rate (global): " << num_all << " / " << denom_all
              << "  =  " << R << " %\n";
  }else{
    std::cout << "Overkill rate (global): n/a (denominator=0)\n";
  }
  std::cout << "\nPer-BH2 segment overkill:\n";
  std::cout << "  h :  numerator / denominator  =  rate(%)\n";
  for(int h=0; h<N_BH2; ++h){
    if(denom_h[h]>0){
      double r = (double)num_h[h] / (double)denom_h[h] * 100.0;
      std::cout << " " << std::setw(2) << h << " :  "
                << std::setw(9) << num_h[h] << " / " << std::setw(9) << denom_h[h]
                << "  =  " << std::setw(7) << r << "\n";
    }else{
      std::cout << " " << std::setw(2) << h << " :  "
                << std::setw(9) << 0 << " / " << std::setw(9) << 0
                << "  =      n/a\n";
    }
  }
  std::cout << "=============================================\n";

  // ===== Draw maps with palette & overlay boxes =====
  if (draw_maps) {
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kViridis);
    gStyle->SetNumberContours(50);
    TGaxis::SetMaxDigits(3);

    TLine grid; grid.SetLineStyle(kDotted); grid.SetLineColor(kGray+1);
    std::vector<TCanvas*> canvases;

    for (int h = 0; h < N_BH2; ++h) {
      if (h % 3 == 0) {
        TString c_name  = TString::Format("c_ud_bh2_%02d_to_%02d", h, std::min(h+2, N_BH2-1));
        TString c_title = TString::Format("BVH_U vs BVH_D | BH2 Segments %d to %d", h, std::min(h+2, N_BH2-1));
        TCanvas* c = new TCanvas(c_name, c_title, 1500, 500);
        c->Divide(3,1);
        canvases.push_back(c);
      }

      TCanvas* c = canvases.back();
      c->cd((h % 3) + 1);
      gPad->SetRightMargin(0.18);          // 팔레트 공간 확보
      gPad->SetTicks(1,1);

      if (use_logz) gPad->SetLogz(true);   // 1카운트도 잘 보이게

      TH2F* H = h_ud[h];
      if (!H) {
        TLatex t; t.SetTextAlign(22); t.DrawLatexNDC(0.5,0.5,"(no data)");
        continue;
      }

      // Z-range
      if (use_fixed_z) {
        H->SetMinimum(use_logz ? 0.5 : 0.0);  // 로그일 땐 0.5
        H->SetMaximum(5.0e4);
      } else {
        H->SetMinimum(use_logz ? 0.5 : 0.0);
      }

      // === Axis titles unified with BVH_3D_ver2.C ===
      H->GetXaxis()->SetTitle("#BVH1 Seg");
      H->GetYaxis()->SetTitle("#BVH2 Seg");
      H->GetZaxis()->SetTitle("Counts");    // ← palette title
      H->Draw("COLZ");
      gPad->Update();

      // palette cosmetics
      if (auto *pal = (TPaletteAxis*)H->GetListOfFunctions()->FindObject("palette")) {
        pal->SetX1NDC(pal->GetX1NDC() - 0.006);
        pal->SetX2NDC(pal->GetX2NDC() - 0.006);
        pal->SetLabelFont(42);
        pal->SetLabelSize(0.020);
        pal->SetTitleSize(0.020);
        pal->SetTitleOffset(0.9);
        H->GetZaxis()->SetLabelSize(0.020);
        H->GetZaxis()->SetTitleSize(0.020);
        H->GetZaxis()->SetTitleOffset(1.0);
      }

      // grid
      const double x_min = H->GetXaxis()->GetXmin(), x_max = H->GetXaxis()->GetXmax();
      const double y_min = H->GetYaxis()->GetXmin(), y_max = H->GetYaxis()->GetXmax();
      for (int i = 0; i <= N_BVHU; ++i) {
        const double x_pos = H->GetXaxis()->GetBinLowEdge(i+1);
        grid.DrawLine(x_pos, y_min, x_pos, y_max);
      }
      for (int j = 0; j <= N_BVHD; ++j) {
        const double y_pos = H->GetYaxis()->GetBinLowEdge(j+1);
        grid.DrawLine(x_min, y_pos, x_max, y_pos);
      }

      // overlay boxes
      for (const auto& tup : overlay_list) {
        int hh, uu, dd; std::tie(hh,uu,dd)=tup;
        if (hh != h) continue;
        const int bx = uu + 1, by = dd + 1;
        if (bx<1 || bx>H->GetNbinsX() || by<1 || by>H->GetNbinsY()) continue;
        const double x1 = H->GetXaxis()->GetBinLowEdge(bx);
        const double x2 = H->GetXaxis()->GetBinUpEdge (bx);
        const double y1 = H->GetYaxis()->GetBinLowEdge(by);
        const double y2 = H->GetYaxis()->GetBinUpEdge (by);
        auto *box = new TBox(x1, y1, x2, y2);
        box->SetFillStyle(0);
        box->SetLineColor(kRed);
        box->SetLineWidth(3);
        box->Draw("SAME");
      }

      gPad->RedrawAxis();

      if (save_png) {
        TString single = TString::Format("UDmap_bh2_%02d.png", h);
        gPad->Update();
        gPad->SaveAs(single);
      }
      if (h % 3 == 2 || h == N_BH2-1) {
        if (save_png) {
          TString out = TString::Format("UDmap_bh2_%02d_to_%02d.png",
                                        h - (h%3), std::min(h - (h%3) + 2,N_BH2-1));
          canvases.back()->SaveAs(out);
        }
      }
    }
  }

  tr->ResetBranchAddresses();
  f->Close();
}

// ================== wrapper ==================
void run_overkill() {
  const char* data_root   = "../E45_2pi_Ver2.root";
  const char* overlay_txt = "overlay_triplets_ver2.txt";
  double ecutBH2 = 0.10, ecutU = 0.04, ecutD = 0.04;

  // 마지막 두 옵션: use_fixed_z(0~50k), use_logz(로그)
  ComputeOverkillFromOverlay(data_root, overlay_txt, ecutBH2, ecutU, ecutD,
                             /*draw_maps=*/true, /*save_png=*/true,
                             /*use_fixed_z=*/true, /*use_logz=*/true);
}

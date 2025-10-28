// -*- C++ -*-
// StatError_10_28.C : Binomial errors & 95% CI (Clopper–Pearson) for Beam-through efficiency and Overkill
// requires ROOT (TEfficiency)

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <cstdio>
#include <TEfficiency.h>

struct DataPoint {         // <-- renamed (avoid clash with MacTypes.h::Point)
  double p_GeV;            // momentum
  long   N_beam;           // denominator: Beam (BH2 in-range)
  double val_pct;          // percent: efficiency(%) for beam-through, or overkill(%) for reactions
};

struct Series {
  std::string name;              // e.g., "tight", "fit", ...
  std::vector<DataPoint> pts;    // points over momenta
};

static void PrintBinomRow(double p_GeV, long N, double val_pct, bool is_efficiency)
{
  (void)is_efficiency; // same math either way; success=val_pct
  const double p = val_pct/100.0;
  const long   k = std::llround(p * N);

  const double se_pct = std::sqrt(p*(1.0-p)/N) * 100.0;

  const double alpha = 0.05; // 95% CL
  const double lo = TEfficiency::ClopperPearson(N, k, 1.0 - alpha, false);
  const double hi = TEfficiency::ClopperPearson(N, k, 1.0 - alpha, true );

  std::printf("  p=%.2f  N=%ld  value=%.3f %%  SE=%.3f %%  CI95=[%.3f, %.3f] %%\n",
              p_GeV, N, 100.0*p, se_pct, 100.0*lo, 100.0*hi);
}

void StatError_10_28()
{
  // (A) π− beam-through (efficiency %)
  std::vector<Series> minus_beam = {
    Series{"tight", {{0.98, 869012, 99.938}, {1.02, 864891, 99.959}, {1.05, 863179, 99.954}, {1.10, 866307, 99.945}, {1.15, 866191, 99.968}}},
    Series{"fit",   {{0.98, 869012, 99.945}, {1.02, 864891, 99.961}, {1.05, 863179, 99.957}, {1.10, 866307, 99.953}, {1.15, 866191, 99.968}}},
    Series{"wide",  {{0.98, 869012, 99.945}, {1.02, 864891, 99.961}, {1.05, 863179, 99.945}, {1.10, 866307, 99.959}, {1.15, 866191, 99.968}}},
    Series{"ultra", {{0.98, 869012, 99.945}, {1.02, 864891, 99.961}, {1.05, 863179, 99.957}, {1.10, 866307, 99.957}, {1.15, 866191, 99.968}}},
  };

  // (B) π+ beam-through (efficiency %)
  std::vector<Series> plus_beam = {
    Series{"tight", {{0.98, 868879, 99.958}, {1.02, 864368, 99.956}, {1.05, 866641, 99.949}, {1.10, 866086, 99.942}, {1.15, 864974, 99.951}}},
    Series{"fit",   {{0.98, 868879, 99.960}, {1.02, 864368, 99.961}, {1.05, 866641, 99.957}, {1.10, 866086, 99.953}, {1.15, 864974, 99.960}}},
    Series{"wide",  {{0.98, 868879, 99.960}, {1.02, 864368, 99.961}, {1.05, 866641, 99.957}, {1.10, 866086, 99.953}, {1.15, 864974, 99.960}}},
    Series{"ultra", {{0.98, 868879, 99.960}, {1.02, 864368, 99.961}, {1.05, 866641, 99.957}, {1.10, 866086, 99.953}, {1.15, 864974, 99.960}}},
  };

  // (C) π−p → π−π0p (Overkill %)
  std::vector<Series> minus_pi0p = {
    Series{"tight", {{0.98,865372,1.528},{1.02,865396,1.541},{1.05,866074,1.569},{1.10,865839,1.652},{1.15,864568,1.665}}},
    Series{"fit",   {{0.98,865372,2.370},{1.02,865396,2.337},{1.05,866074,2.364},{1.10,865839,2.466},{1.15,864568,2.458}}},
    Series{"wide",  {{0.98,865372,3.134},{1.02,865396,3.123},{1.05,866074,3.152},{1.10,865839,3.229},{1.15,864568,3.297}}},
    Series{"ultra", {{0.98,865372,3.602},{1.02,865396,3.620},{1.05,866074,3.630},{1.10,865839,3.719},{1.15,864568,3.777}}},
  };

  // (D) π−p → π−π+n (Overkill %)
  std::vector<Series> minus_piplusn = {
    Series{"tight", {{0.98,865535,0.826},{1.02,865843,0.828},{1.05,866032,0.839},{1.10,863873,0.880},{1.15,866719,0.925}}},
    Series{"fit",   {{0.98,865535,1.297},{1.02,865843,1.283},{1.05,866032,1.323},{1.10,863873,1.459},{1.15,866719,1.416}}},
    Series{"wide",  {{0.98,865535,1.719},{1.02,865843,1.736},{1.05,866032,1.769},{1.10,863873,1.937},{1.15,866719,1.977}}},
    Series{"ultra", {{0.98,865535,2.120},{1.02,865843,2.136},{1.05,866032,2.188},{1.10,863873,2.350},{1.15,866719,2.399}}},
  };

  // (E) π+p → π+π+n (Overkill %)
  std::vector<Series> plus_piplusn = {
    Series{"tight", {{0.98,840755,0.872},{1.02,842613,0.892},{1.05,844403,0.939},{1.10,847494,0.987},{1.15,848225,1.038}}},
    Series{"fit",   {{0.98,840755,1.251},{1.02,842613,1.286},{1.05,844403,1.341},{1.10,847494,1.389},{1.15,848225,1.419}}},
    Series{"wide",  {{0.98,840755,1.687},{1.02,842613,1.732},{1.05,844403,1.784},{1.10,847494,1.832},{1.15,848225,1.888}}},
    Series{"ultra", {{0.98,840755,1.938},{1.02,842613,1.990},{1.05,844403,2.056},{1.10,847494,2.109},{1.15,848225,2.195}}},
  };

  // (F) π+p → π+π0p (Overkill %)
  std::vector<Series> plus_pi0p = {
    Series{"tight", {{0.98,841829,1.492},{1.02,841769,1.423},{1.05,845908,1.470},{1.10,848408,1.516},{1.15,849493,1.559}}},
    Series{"fit",   {{0.98,841829,1.886},{1.02,841769,1.848},{1.05,845908,1.916},{1.10,848408,1.969},{1.15,849493,2.057}}},
    Series{"wide",  {{0.98,841829,2.721},{1.02,841769,2.737},{1.05,845908,2.839},{1.10,848408,2.877},{1.15,849493,2.946}}},
    Series{"ultra", {{0.98,841829,2.985},{1.02,841769,3.018},{1.05,845908,3.177},{1.10,848408,3.223},{1.15,849493,3.296}}},
  };

  auto run_block = [](const char* title, const std::vector<Series>& S, bool is_eff){
    std::cout << "\n==== " << title << " ====\n";
    for (const auto& s : S){
      std::cout << "[" << s.name << "]\n";
      for (const auto& pt : s.pts){
        PrintBinomRow(pt.p_GeV, pt.N_beam, pt.val_pct, is_eff);
      }
    }
  };

  run_block("π− beam-through (EFF %)", minus_beam, true);
  run_block("π+ beam-through (EFF %)", plus_beam,  true);
  run_block("π−p → π−π0p (Overkill %)",  minus_pi0p, false);
  run_block("π−p → π−π+n (Overkill %)",  minus_piplusn,false);
  run_block("π+p → π+π+n (Overkill %)",  plus_piplusn, false);
  run_block("π+p → π+π0p (Overkill %)",  plus_pi0p,    false);

  std::cout << "\n[Tip] SE는 대략적 시각화용, 보고/논문에는 Clopper–Pearson 95% CI 권장.\n";
}

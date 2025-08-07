// -*- C++ -*-

#include "DiffCrossSection.hh"

#include <Randomize.hh>

#include "FuncName.hh"
#include "Kinematics.hh"
#include "DiffCrossSectionMan.hh"

namespace
{
const auto& gDcsMan = DiffCrossSectionMan::GetInstance();
}

namespace DiffCrossSection
{

//______________________________________________________________________________
std::vector<G4double>
CoeffLinearInterpolation(const G4double mom_kaon, const std::vector<G4double>& mom_kaons, const std::unordered_map<G4double, std::vector<G4double>>& legendre_coeff)
{
  // -- check whick value mom_kaon is between -----
  auto it = std::lower_bound(mom_kaons.begin(), mom_kaons.end(), mom_kaon);
  G4int lower_index, upper_index;
  if (it == mom_kaons.begin()) {
    lower_index = 0;
    upper_index = 0;
  } else if (it == mom_kaons.end()) {
    lower_index = mom_kaons.size()-1;
    upper_index = mom_kaons.size()-1;
  } else {
    lower_index = std::distance(mom_kaons.begin(), it - 1);
    upper_index = std::distance(mom_kaons.begin(), it);
  }

  // -- interpolate and estimate coefficient -----
  std::vector<G4double> coeff = legendre_coeff.at(mom_kaons[lower_index]);
  if (lower_index != upper_index) {
    std::vector<G4double> coeff1 = legendre_coeff.at(mom_kaons[lower_index]);
    std::vector<G4double> coeff2 = legendre_coeff.at(mom_kaons[upper_index]);
    for (G4int order = 0, n_order = coeff1.size(); order < n_order; order++) {
      G4double a = (coeff2[order] - coeff1[order]) / (mom_kaons[upper_index] - mom_kaons[lower_index]);
      G4double b = coeff2[order] - a*mom_kaons[upper_index];
      coeff[order] = a*mom_kaon + b;
    }
  }

  return coeff;
}

//______________________________________________________________________________
std::vector<G4double>
CoeffSpline(const G4double mom_kaon)
{
  std::vector<G4double> coeff;
  G4int n_order = gDcsMan.GetNumSpline();
  for (G4int order = 0; order < n_order; order++) {
    coeff.push_back( gDcsMan.GetCoeff(order, mom_kaon) );
  }

  return coeff;
}

  
//______________________________________________________________________________
G4bool
RejectionSampling(const G4double cos_theta, const std::vector<G4double>& coeff)
{
  // -- search maximum value -----
  G4double maximum_value = 0.0;
  G4double range_min = -1.0;
  G4double range_max =  1.0;
  G4int n_iter   = 5;
  G4int n_points = 8;
  for (G4int _ = 0; _ < n_iter; _++) {
    G4double step = (range_max - range_min) / n_points;
    std::vector<G4double> diff_cs;
    for (G4int i = 0; i < n_points+1; i++) {
      G4double tmp_cos_theta = range_min + i*step;
      G4double tmp_diff_cs = 0.0;
      for (G4int order = 0, n_order = coeff.size(); order < n_order; order++) tmp_diff_cs += coeff[order]*Kinematics::Legendre(order, tmp_cos_theta);
      diff_cs.push_back( tmp_diff_cs );
    }
    auto max_it = std::max_element(diff_cs.begin(), diff_cs.end());
    G4int index = std::distance(diff_cs.begin(), max_it);
    maximum_value = *max_it;
    range_min = index==0 ? -1.0 : range_min + (index-1.0)*step;
    range_max = index==diff_cs.size()-1 ? 1.0 : range_min + (index+1.0)*step;
  }

  // -- calc. diff cross section and do rejection sampling -----
  G4double eval_diff_cs = 0.0;
  for (G4int order = 0, n_order = coeff.size(); order < n_order; order++) eval_diff_cs += coeff[order]*Kinematics::Legendre(order, cos_theta);
  return G4RandFlat::shoot(0.0, maximum_value) <= eval_diff_cs;
}


//______________________________________________________________________________
G4bool
EtaLambda(const G4double cos_theta, const G4double mom_kaon)
{
  std::vector<G4double> coeff;

  // -- Crystal Ball -----
  if (gDcsMan.IsReady()) {
    coeff = CoeffSpline(mom_kaon);
  } else {
    std::vector<G4double> mom_kaons = DiffCrossSection::mom_kaons_etaLambda_CB;
    std::unordered_map<G4double, std::vector<G4double>> legendre_coeff = DiffCrossSection::legendre_coeff_etaLambda_CB;
    coeff = CoeffLinearInterpolation(mom_kaon, mom_kaons, legendre_coeff);
  }
  
  return RejectionSampling(cos_theta, coeff);
}

//______________________________________________________________________________
G4bool
PiZeroLambda(const G4double cos_theta, const G4double mom_kaon)
{
  std::vector<G4double> coeff;

  // -- bubble chamber 1970 -----
  if (gDcsMan.IsReady()) {
    coeff = CoeffSpline(mom_kaon);
  } else {
    std::vector<G4double> mom_kaons = DiffCrossSection::mom_kaons_pi0Lambda_bubble1970;
    std::unordered_map<G4double, std::vector<G4double>> legendre_coeff = DiffCrossSection::legendre_coeff_pi0Lambda_bubble1970;
    coeff = CoeffLinearInterpolation(mom_kaon, mom_kaons, legendre_coeff);
  }
  
  return RejectionSampling(cos_theta, coeff);
}

//______________________________________________________________________________
G4bool
PiZeroSigmaZero(const G4double cos_theta, const G4double mom_kaon)
{
  std::vector<G4double> coeff;
    
  // -- bubble chamber 1970 -----
  if (gDcsMan.IsReady()) {
    coeff = CoeffSpline(mom_kaon);
  } else {
    std::vector<G4double> mom_kaons = DiffCrossSection::mom_kaons_pi0Sigma0_bubble1970;
    std::unordered_map<G4double, std::vector<G4double>> legendre_coeff = DiffCrossSection::legendre_coeff_pi0Sigma0_bubble1970;
    coeff = CoeffLinearInterpolation(mom_kaon, mom_kaons, legendre_coeff);
  }
  
  return RejectionSampling(cos_theta, coeff);
}

//______________________________________________________________________________
G4bool
PiPlusSigmaMinus(const G4double cos_theta, const G4double mom_kaon)
{
  std::vector<G4double> coeff;

  // -- bubble chamber 1970 -----
  if (gDcsMan.IsReady()) {
    coeff = CoeffSpline(mom_kaon);
  } else {
    std::vector<G4double> mom_kaons = DiffCrossSection::mom_kaons_piPlusSigmaMinus_bubble1970;
    std::unordered_map<G4double, std::vector<G4double>> legendre_coeff = DiffCrossSection::legendre_coeff_piPlusSigmaMinus_bubble1970;
    coeff = CoeffLinearInterpolation(mom_kaon, mom_kaons, legendre_coeff);
  }
  
  return RejectionSampling(cos_theta, coeff);
}

//______________________________________________________________________________
G4bool
PiMinusSigmaPlus(const G4double cos_theta, const G4double mom_kaon)
{
  std::vector<G4double> coeff;

  // -- bubble chamber 1970 -----
  if (gDcsMan.IsReady()) {
    coeff = CoeffSpline(mom_kaon);
  } else {
    std::vector<G4double> mom_kaons = DiffCrossSection::mom_kaons_piMinusSigmaPlus_bubble1970;
    std::unordered_map<G4double, std::vector<G4double>> legendre_coeff = DiffCrossSection::legendre_coeff_piMinusSigmaPlus_bubble1970;
    coeff = CoeffLinearInterpolation(mom_kaon, mom_kaons, legendre_coeff);
  }
  
  return RejectionSampling(cos_theta, coeff);
}

//______________________________________________________________________________
G4bool
KaonMinusProton(const G4double cos_theta, const G4double mom_kaon)
{
  std::vector<G4double> coeff;

  // -- bubble chamber 1970 -----
  if (gDcsMan.IsReady()) {
    coeff = CoeffSpline(mom_kaon);
  } else {
    std::vector<G4double> mom_kaons = DiffCrossSection::mom_kaons_Kp_bubble1970;
    std::unordered_map<G4double, std::vector<G4double>> legendre_coeff = DiffCrossSection::legendre_coeff_Kp_bubble1970;
    coeff = CoeffLinearInterpolation(mom_kaon, mom_kaons, legendre_coeff);
  }
  
  return RejectionSampling(cos_theta, coeff);
}

//______________________________________________________________________________
G4bool
KaonZeroNeutron(const G4double cos_theta, const G4double mom_kaon)
{
  std::vector<G4double> coeff;

  // -- bubble chamber 1970 -----
  if (gDcsMan.IsReady()) {
    coeff = CoeffSpline(mom_kaon);
  } else {
    std::vector<G4double> mom_kaons = DiffCrossSection::mom_kaons_K0n_bubble1970;
    std::unordered_map<G4double, std::vector<G4double>> legendre_coeff = DiffCrossSection::legendre_coeff_K0n_bubble1970;
    coeff = CoeffLinearInterpolation(mom_kaon, mom_kaons, legendre_coeff);
  }
  
  return RejectionSampling(cos_theta, coeff);
}

  
} //namespace DiffCrossSection

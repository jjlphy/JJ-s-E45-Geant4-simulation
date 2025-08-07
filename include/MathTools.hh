// -*- C++ -*-

#ifndef MATH_TOOLS_HH
#define MATH_TOOLS_HH

#include <limits>

#include <CLHEP/Units/SystemOfUnits.h>
#include <G4String.hh>

//_____________________________________________________________________________
namespace math
{
static const G4double Infinity = std::numeric_limits<G4double>::infinity();
static const G4double Epsilon = std::numeric_limits<G4double>::epsilon();

//___________________________________________________________________________
inline G4String
ClassName()
{
  static G4String s_name("MathTools");
  return s_name;
}

//___________________________________________________________________________
inline G4double
Deg2Rad()
{
  static const G4double s_v = CLHEP::pi/180.;
  return s_v;
}

//___________________________________________________________________________
inline G4double
Rad2Deg()
{
  static const G4double s_v = 180./CLHEP::pi;
  return s_v;
}

//___________________________________________________________________________
inline G4double
Sq(G4double x)
{
  return x*x;
}

//___________________________________________________________________________
inline G4double
pythag(G4double a, G4double b)
{
  G4double aa = std::abs(a);
  G4double ab = std::abs(b);
  if(aa>ab)
    return aa*std::sqrt(1. + (ab/aa)*(ab/aa));
  else if(ab!=0.)
    return ab*std::sqrt(1. + (aa/ab)*(aa/ab));
  else
    return 0.;
}

G4bool GaussElim(G4double **a, G4int n, G4double *b, G4int *indx, G4int *ipiv);

G4bool GaussJordan(G4double **a, G4int n, G4double *b,
                   G4int *indxc, G4int *indxd, G4int *ipiv);

G4bool InterpolateRatio(G4int n, const G4double *xa, const G4double *ya,
                        G4double *w1, G4double *w2,
                        G4double x, G4double &y, G4double &dy);

G4bool InterpolatePol(G4int n, const G4double *xa, const G4double *ya,
                      G4double *w1, G4double *w2,
                      G4double x, G4double &y, G4double &dy);

G4bool SVDksb(G4double **u, const G4double *w, G4double **v,
              G4int m, G4int n, const G4double *b, G4double *x, G4double *wv);

G4bool SVDcmp(G4double **a, G4int m, G4int n, G4double *w, G4double **v,
              G4double *wv);
}

#endif

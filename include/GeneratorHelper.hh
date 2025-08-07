// -*- C++ -*-

/**
 * for the EventGeneration for the E27 experiment
 */

#ifndef GENERATOR_HELPER_HH
#define GENERATOR_HELPER_HH

#include "G4ThreeVector.hh"

class AngDisGenerator;

G4ThreeVector UniformDirectionInUV(double u0, double v0,
                                   double hu=0, double hv=0);
G4ThreeVector GaussDirectionInUV(double u0, double v0,
                                 double su, double sv);
G4ThreeVector UniformPosition(double hx, double hy, double hz);
G4ThreeVector GaussPosition(double sx, double sy, double hz);
G4ThreeVector GaussPosition_LqTarg(double x0, double y0, double z0,
                                   double dx, double dy, double targ_r,
                                   double targ_height);
G4ThreeVector UniformDirectionInThetaPhi(double cost1, double cost2);
double BreitWigner(double mean, double gamma);
bool Decay2Body(double Mini, double Mf1, double Mf2,
                const G4ThreeVector & Pini,
                G4ThreeVector & Pf1,  G4ThreeVector & Pf2,
                const AngDisGenerator & generator);
bool Scattering2Body_theta(double Mi1, double Mi2, double Mf1, double Mf2,
                           const G4ThreeVector & Pini1,
                           const G4ThreeVector & Pini2,
                           G4ThreeVector & Pf1,  G4ThreeVector & Pf2,
                           double & theta_CM,
                           const AngDisGenerator & generator);
bool Scattering3Body_theta(double Mi1, double Mi2, double Mf1, double Mf2,
                           double Mf3,
                           const G4ThreeVector & Pini1,
                           const G4ThreeVector & Pini2,
                           G4ThreeVector & Pf1,  G4ThreeVector & Pf2,
                           G4ThreeVector & Pf3,
                           double & theta_CM,
                           const AngDisGenerator & generator);
bool Decay3BodyPhaseSpace(double Mini, double Mf1, double Mf2, double Mf3,
                          const G4ThreeVector & Pini,
                          G4ThreeVector & Pf1,  G4ThreeVector & Pf2,
                          G4ThreeVector & Pf3);
G4int Getnp_JAM(int i_num);
G4int GetPID_JAM(int i_num, int inp);
G4ThreeVector GetP_JAM(int i_num, int inp);

#endif

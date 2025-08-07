// -*- C++ -*-

#ifndef KINEMA_2BODY_HH
#define KINEMA_2BODY_HH

#include "KinemaCM2Lab.hh"

//_____________________________________________________________________________
class Kinema2Body
{
public:
  Kinema2Body();
  Kinema2Body(double m1, double m2, double m3, double m4);

private:
  KinemaCM2Lab kin;
  double       Phi;

public:
  double Beta2Gamma(double beta) const;
  int    CalcKinema();
  void   Dump() const;
  double GetEnergyCM(int i) const;
  double GetEnergyLab(int i) const;
  double GetMass(int i) const;
  double GetMomentumCM(int i) const;
  double GetMomentumLab(int i) const;
  double GetPhiLab() const;
  double GetTheta(int i) const;
  double GetThetaCM() const;
  double GetThetaLab() const;
  double E2p(double E, double m) const;
  double p2E(double p, double m) const;
  double pE2beta(double p,double E,double m_2) const;
  double pE2beta(double p1,double E1,double p2, double E2) const;
  double pE2beta(double p1,double E1,double m1, double p2,
                 double E2, double m2) const;
  void   SetMass(double *mass);
  void   SetMass(int i, double mass);
  void   SetMomentum(int i,double mom);
  void   SetThetaCM(double theta_cm);
  void   SetTheta(int i, double theta);
  double Theta2ThetaCM(double theta, double p, double E,
                       double kin_gamma, double beta) const;
  double ThetaCM2ThetaLab(double theta, double p,double E,
                          double gamma_cm,double beta_cm) const;
  double ThetaCM2PhiLab(double theta, double p, double E,
                        double gamma_cm, double beta_cm) const;
};

#endif

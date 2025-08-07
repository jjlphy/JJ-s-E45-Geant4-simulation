// -*- C++ -*-

#ifndef KINEMA_H_WEAK_NONRESO_HH
#define KINEMA_H_WEAK_NONRESO_HH

#include "Kinema2Body.hh"
#include "Kinema3Body.hh"

//_____________________________________________________________________________
struct KINEMA_HWEAKNONESO
{
  double E_1_lab;
  double p_1_lab;
  double M_1;

  double E_2_lab;
  double p_2_lab;
  double M_2;

  double E_3_lab;
  double p_3_lab;
  double M_3;
  double P_3_lab[3];
  double theta3, phi3;

  double E_4_lab;
  double p_4_lab;
  double M_4;
  double P_4_lab[3];
  double theta4, phi4;

  double E_5_lab;
  double p_5_lab;
  double M_5;
  double P_5_lab[3];
  double theta5, phi5;

  double E_6_lab;
  double p_6_lab;
  double M_6;
  double P_6_lab[3];
  double theta6, phi6;

  double E_res_lab;
  double p_res_lab;
  double M_res;
  double P_res_lab[3];
  double theta_res, phi_res;

  double Theta1CM, Theta2CM;
  double Phi1,Phi2;
};

//_____________________________________________________________________________
class KinemaHWeakNonReso
{
public:
  KinemaHWeakNonReso(double m1, double m2, double m3, double m4, double m5,
                     double m_res, double width, double p1, double p2);
  ~KinemaHWeakNonReso();

private:
  Kinema2Body kin1;
  // Kinema3Body kin2; // change.... shhwang
  KINEMA_HWEAKNONESO kin3;

public:
  double p2E(double p,double m);
  void   CalcDistoribution(double unitx, double unity, double unitz, double *theta, double *phi);
  double deg2rad(double theta);
  double rag2deg(double rag);
  double RandSin();
  void   Dump();
  double GetEnergy(int i);
  double GetMomentum(int i);
  void   GetMomentum(int i, double *mom);
  double GetTheta(int i);
  double GetPhi(int i);
  double GetThetaCM(int i);
  double GetPhiCM(int i);
  void   RotateMom(int i, double deg, double *mom);
  double GetResMass();
};

#endif

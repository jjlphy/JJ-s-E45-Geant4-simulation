// -*- C++ -*-

#ifndef KINEMA_3BODY_HH
#define KINEMA_3BODY_HH

#include "Kinema2Body.hh"

//_____________________________________________________________________________
struct KINEMA_3BODY
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

  double Theta1CM, Theta2CM;
  double Phi1,Phi2;

  double m34;
};

//_____________________________________________________________________________
class Kinema3Body
{
public:
  Kinema3Body();
  Kinema3Body(double m1, double m2, double m3, double m4, double m5,
	       double p1, double p2);
  ~Kinema3Body();

private:
  Kinema2Body kin1;
  Kinema2Body kin2;
  KINEMA_3BODY kin3;

public:
  double p2E(double p,double m);
  void CalcDistoribution(double unitx, double unity, double unitz, double *theta, double *phi);
  double deg2rad(double theta);
  double rag2deg(double rag);
  double CalcM12(double m1, double m2, double m3, double M);
  double GetMaxFunc(double m1, double m2, double m3, double M);
  double ProFunction(double m12, double m1, double m2, double m3, double M);
  double RandSin();
  void Dump();
  double GetEnergy(int i);
  double GetMomentum(int i);
  double GetMomentum(int i, double *mom);
  double GetTheta(int i);
  double GetPhi(int i);
  double GetThetaCM(int i);
  double GetPhiCM(int i);
  double GetM34();
  void  RotateMom(int i, double deg, double *mom);
};

#endif

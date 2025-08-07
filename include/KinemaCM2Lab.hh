// -*- C++ -*-

#ifndef KINEMA_CM2LAB_HH
#define KINEMA_CM2LAB_HH

//_____________________________________________________________________________
struct KinemaCM2Lab
{
  double E_1_lab;
  double p_1_lab;
  double p_1_theta;
  double M_1;
  double E_1_cm;

  double p_12_cm;

  double E_2_lab;
  double p_2_lab;
  double p_2_theta;
  double M_2;
  double E_2_cm;

  double beta_cm;
  double gamma_cm;

  double E_3_cm;
  double E_4_cm;
  double p_34_cm;
  double M_3;
  double M_4;

  double theta_cm;

  double E_3_lab;
  double p_3_lab;

  double E_4_lab;
  double p_4_lab;

  double theta_lab;
  double phi_lab;
};


#endif

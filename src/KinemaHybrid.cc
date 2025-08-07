// -*- C++ -*-

#include "KinemaHybrid.hh"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <CLHEP/Units/SystemOfUnits.h>

#include <G4ios.hh>
#include <Randomize.hh>

//_____________________________________________________________________________
KinemaHybrid::KinemaHybrid(double m1, double m2, double m3, double m4,
                           double m5, double m_res, double width,
                           double p1, double p2)
{
  double ECM;
  double vx_res, vy_res, vz_res;   /* unit vector */
  // double vx3, vy3, vz3;            /* unit vector */
  // double vx4, vy4, vz4;            /* unit vector */
  // double vx5, vy5, vz5;            /* unit vector */
  // double theta1, theta2;           /* tmpolary; theta1 represents
  // 				      1st kinematics of m1, m2, m_res, and m5.
  // 				      theta2 represents 2nd kinematics of
  // 				      m_res, m3, and m4 */
  // double phi3;                     /* 2phi(CM system)*/
  double phi6;                     /* log note p.43 */
  double theta_res;                /* log note p.43 */
  // double Theta3,Phi3;              /*heta,Ph*/
  // double Theta4,Phi4;              /* eta,Ph*/
  // double Theta5,Phi5;              /* Ph*/
  double Theta_res,Phi_res;        /* Theta,Ph*/

  // double m12;                      /* */
  double mom[3];                      /* */

  //  G4cout<<"m1"<<m1<<G4endl;
  //  G4cout<<"m2"<<m2<<G4endl;
  //  G4cout<<"m3"<<m3<<G4endl;
  //  G4cout<<"m4"<<m4<<G4endl;
  //  G4cout<<"m5"<<m5<<G4endl;
  //  G4cout<<"m_res"<<m_res<<G4endl;
  //  G4cout<<"shhwang test test 1"<<G4endl;

  kin3.M_1 = m1; //--> beam
  kin3.M_2 = m2; //--> target
  kin3.M_3 = m3; //--> pro decay 1
  kin3.M_4 = m4; //--> pi1 decay 2
  kin3.M_5 = m5; //--> pi2 decay 3

  kin3.p_1_lab = p1;
  kin3.p_2_lab = p2;

  kin3.E_1_lab = p2E(kin3.p_1_lab, kin3.M_1);
  kin3.E_2_lab = p2E(kin3.p_2_lab, kin3.M_2);

  ECM = sqrt((kin3.E_1_lab+kin3.E_2_lab)*(kin3.E_1_lab+kin3.E_2_lab)
	     -(kin3.p_1_lab+kin3.p_2_lab)*(kin3.p_1_lab+kin3.p_2_lab));
  do {
    kin3.M_res = CLHEP::RandBreitWigner::shoot(m_res, width);
  } while (kin3.M_3+kin3.M_4 > kin3.M_res || kin3.M_res > ECM-kin3.M_5);
  //  kin3.M_res = m_res;
  kin1 = Kinema2Body(m1, m2, kin3.M_res, 0);
  kin1.SetMomentum(1, p1);
  kin1.SetMomentum(2, p2);
  //  kin1.SetThetaCM((double)RandSin()); //set to zero
  kin1.SetThetaCM(0.); //set to zero
  kin1.CalcKinema();
  phi6 = (360.0*(double)CLHEP::RandFlat::shoot());

  kin3.Theta1CM = kin1.GetThetaCM();
  kin3.Phi1     = phi6;
  /* calculate m_res */
  theta_res = kin1.GetThetaLab();

  vx_res = sin(CLHEP::pi*theta_res/180.0)*cos(CLHEP::pi*phi6/180.0);
  vy_res = -sin(deg2rad(theta_res))*sin(deg2rad(phi6));
  vz_res = cos(CLHEP::pi*theta_res/180.0);

  CalcDistoribution(vx_res, vy_res, vz_res, &Theta_res, &Phi_res);

  kin3.E_res_lab = kin1.GetEnergyLab(3);
  kin3.p_res_lab = kin1.GetMomentumLab(3);

  kin3.P_res_lab[0] = kin3.p_res_lab*vx_res;
  kin3.P_res_lab[1] = kin3.p_res_lab*vy_res;
  kin3.P_res_lab[2] = kin3.p_res_lab*vz_res;

  kin3.theta_res = Theta_res;
  kin3.phi_res = Phi_res;

  /* calculate m3,m4, m5 */
  if (kin3.M_res == m3 && m4 == 0.0 && m5 == 0.0) {
    kin3.E_4_lab = 0.0;
    kin3.p_4_lab = 0.0;
    kin3.P_4_lab[0] = 0.0;
    kin3.P_4_lab[1] = 0.0;
    kin3.P_4_lab[2] = 0.0;
    kin3.theta4 = 0.0;
    kin3.phi4 = 0.0;

    kin3.E_5_lab = 0.0;
    kin3.p_5_lab = 0.0;
    kin3.P_5_lab[0] = 0.0;
    kin3.P_5_lab[1] = 0.0;
    kin3.P_5_lab[2] = 0.0;
    kin3.theta5 = 0.0;
    kin3.phi5 = 0.0;

    kin3.E_3_lab = kin3.E_res_lab;
    kin3.p_3_lab =  kin3.p_res_lab;
    kin3.P_3_lab[0] = kin3.P_res_lab[0];
    kin3.P_3_lab[1] = kin3.P_res_lab[1];
    kin3.P_3_lab[2] = kin3.P_res_lab[2];
    kin3.theta3 = Theta_res;
    kin3.phi3 = Phi_res;
  } else {
    //// I need to change the 3body decay mode;;
    //kin3.p_res_lab = kin1.GetMomentumLab(3);
    double temp= kin1.GetMomentumLab(3);
    //    kin2 = Kinema2Body(kin3.M_res, 0.0, m3, m4);
    Kinema3Body kin2(kin3.M_res, 0.0, m3, m4, m5, temp, 0.0);
    //    Kinema3Body kin2(kin3.M_res, 0.0, m3, m4, m5, temp, 0.0);

    /* m3 */
    kin3.E_3_lab = kin2.GetEnergy(3);
    kin3.p_3_lab = kin2.GetMomentum(3);

    kin3.theta3 = kin2.GetTheta(3);
    kin3.phi3 = kin2.GetPhi(3);

    kin2.GetMomentum(3,mom);
    kin3.P_3_lab[0] = mom[0];
    kin3.P_3_lab[1] = mom[1];
    kin3.P_3_lab[2] = mom[2];

    /* m4 */
    kin3.E_4_lab = kin2.GetEnergy(4);
    kin3.p_4_lab = kin2.GetMomentum(4);
    kin2.GetMomentum(4,mom);
    kin3.P_4_lab[0] = mom[0];
    kin3.P_4_lab[1] = mom[1];
    kin3.P_4_lab[2] = mom[2];

    kin3.theta4 = kin2.GetTheta(4);
    kin3.phi4 = kin2.GetPhi(4);


    /* m5 */
    kin3.E_5_lab = kin2.GetEnergy(5);
    kin3.p_5_lab = kin2.GetMomentum(5);
    kin2.GetMomentum(5,mom);
    kin3.P_5_lab[0] = mom[0];
    kin3.P_5_lab[1] = mom[1];
    kin3.P_5_lab[2] = mom[2];

    kin3.theta5 = kin2.GetTheta(5);
    kin3.phi5 = kin2.GetPhi(5);

    /*
    // m3 //
    theta1 = kin1.GetThetaLab();
    theta2 = kin2.GetThetaLab();

    vx3 = cos(deg2rad(theta2))*cos(deg2rad(theta1)) -
    sin(deg2rad(theta1))*cos(deg2rad(phi3))*sin(deg2rad(theta2));

    vy3 = cos(deg2rad(phi6))*cos(deg2rad(theta2))*sin(deg2rad(theta1)) +
    cos(deg2rad(theta1))*cos(deg2rad(phi6))*cos(deg2rad(phi3))*sin(deg2rad(theta2)) -
    sin(deg2rad(phi6))*sin(deg2rad(phi3))*sin(deg2rad(theta2));

    vz3 = -sin(deg2rad(phi6))*cos(deg2rad(theta2))*sin(deg2rad(theta1)) -
    cos(deg2rad(theta1))*sin(deg2rad(phi6))*cos(deg2rad(phi3))*sin(deg2rad(theta2)) -
    cos(deg2rad(phi6))*sin(deg2rad(phi3))*sin(deg2rad(theta2));

    CalcDistoribution(vx3, vy3, vz3, &Theta3, &Phi3);

    kin3.E_3_lab = kin2.GetEnergyLab(3);
    kin3.p_3_lab = kin2.GetMomentumLab(3);
    kin3.P_3_lab[0] = kin3.p_3_lab*vx3;
    kin3.P_3_lab[1] = kin3.p_3_lab*vy3;
    kin3.P_3_lab[2] = kin3.p_3_lab*vz3;
    kin3.theta3 = Theta3;
    kin3.phi3 = Phi3;

    // m4 //
    theta1 = kin1.GetThetaLab();
    theta2 = -kin2.GetPhiLab();

    vx4 = cos(deg2rad(theta2))*cos(deg2rad(theta1)) -
    sin(deg2rad(theta1))*cos(deg2rad(phi3))*sin(deg2rad(theta2));

    vy4 = cos(deg2rad(phi6))*cos(deg2rad(theta2))*sin(deg2rad(theta1)) +
    cos(deg2rad(theta1))*cos(deg2rad(phi6))*cos(deg2rad(phi3))*sin(deg2rad(theta2)) -
    sin(deg2rad(phi6))*sin(deg2rad(phi3))*sin(deg2rad(theta2));

    vz4 = -sin(deg2rad(phi6))*cos(deg2rad(theta2))*sin(deg2rad(theta1)) -
    cos(deg2rad(theta1))*sin(deg2rad(phi6))*cos(deg2rad(phi3))*sin(deg2rad(theta2)) -
    cos(deg2rad(phi6))*sin(deg2rad(phi3))*sin(deg2rad(theta2));

    CalcDistoribution(vx4, vy4, vz4, &Theta4, &Phi4);

    kin3.E_4_lab = kin2.GetEnergyLab(4);
    kin3.p_4_lab = kin2.GetMomentumLab(4);
    kin3.P_4_lab[0] = kin3.p_4_lab*vx4;
    kin3.P_4_lab[1] = kin3.p_4_lab*vy4;
    kin3.P_4_lab[2] = kin3.p_4_lab*vz4;
    kin3.theta4 = Theta4;
    kin3.phi4 = Phi4;
    */
  }
  //Dump();
}

//_____________________________________________________________________________
KinemaHybrid::~KinemaHybrid()
{
}

//_____________________________________________________________________________
double
KinemaHybrid::p2E(double p,double m)
{
  return sqrt(p*p + m*m);
}

//_____________________________________________________________________________
void
KinemaHybrid::CalcDistoribution(double unitx, double unity, double unitz,
                                double *theta, double *phi)
{
  *theta = rag2deg(acos(unitx));

  if (unity>=0.0 && unitz>0.0)
    *phi = rag2deg(acos(unity/sin(deg2rad(*theta))));
  else if (unity<0.0 && unitz>=0.0)
    *phi = rag2deg(acos(unity/sin(deg2rad(*theta))));
  else if (unity<=0.0 && unitz<0.0)
    *phi = 360.0-rag2deg(acos(unity/sin(deg2rad(*theta))));
  else if (unity>0.0 && unitz<=0.0)
    *phi = 360.0-rag2deg(acos(unity/sin(deg2rad(*theta))));
  else {
    fprintf(stderr,
            "KinemaHybrid::CalcDistribution No such reagion unity=%f, unitz=%f\n",
	    unity, unitz);
    Dump();
    exit(1);
  }
  return;
}

//_____________________________________________________________________________
double KinemaHybrid::deg2rad(double theta)
{
  return CLHEP::pi*theta/180.0;
}

//_____________________________________________________________________________
double
KinemaHybrid::rag2deg(double rag)
{
  return 360.0 * rag/ (2.0 * CLHEP::pi);
}

//_____________________________________________________________________________
double
KinemaHybrid::RandSin()
{
  int success=0;
  double x,fx;

  do {
    x = 180.0 * (double)CLHEP::RandFlat::shoot();
    //x = 180.0 * (double)rand()/(RAND_MAX+1.0);
    fx = sin(deg2rad(x));
    if (fx >= (double)CLHEP::RandFlat::shoot())
      success = 1;
  } while (success==0);

  return x;
}

//_____________________________________________________________________________
void
KinemaHybrid::Dump()
{
  printf("======KinemaHybrid Dump======\n");
  printf("--Particle1--\n");
  printf("mass=%f, p_lab=%f, E_lab=%f\n",kin3.M_1, kin3.p_1_lab, kin3.E_1_lab);
  printf("--Particle2--\n");
  printf("mass=%f, p_lab=%f, E_lab=%f\n",kin3.M_2, kin3.p_2_lab, kin3.E_2_lab);
  printf("--Resonance--\n");
  printf("mass=%f\n",kin3.M_res);
  printf("--Particle3--\n");
  printf("mass=%f, p_lab=%f, E_lab=%f\n",kin3.M_3, kin3.p_3_lab, kin3.E_3_lab);
  printf("momentum=(%f, %f, %f), (theta, phi)=(%f, %f)\n",
	 kin3.P_3_lab[0], kin3.P_3_lab[1], kin3.P_3_lab[2],
	 kin3.theta3, kin3.phi3);
  printf("--Particle4--\n");
  printf("mass=%f, p_lab=%f, E_lab=%f\n",kin3.M_4, kin3.p_4_lab, kin3.E_4_lab);
  printf("momentum=(%f, %f, %f), (theta, phi)=(%f, %f)\n",
	 kin3.P_4_lab[0], kin3.P_4_lab[1], kin3.P_4_lab[2],
	 kin3.theta4, kin3.phi4);
  printf("--Particle5--\n");
  printf("mass=%f, p_lab=%f, E_lab=%f\n",kin3.M_5, kin3.p_5_lab, kin3.E_5_lab);
  printf("momentum=(%f, %f, %f), (theta, phi)=(%f, %f)\n",
	 kin3.P_5_lab[0], kin3.P_5_lab[1], kin3.P_5_lab[2],
	 kin3.theta5, kin3.phi6);

  printf("Energy:E1+E2=%f, E1+E2+E3=%f\n", kin3.E_1_lab+kin3.E_2_lab,
	 kin3.E_3_lab+kin3.E_4_lab+kin3.E_5_lab);
  printf("Momentum: x-> p1+p2=%f, p3+p4+p5=%f\n", kin3.p_1_lab+kin3.p_2_lab,
	 kin3.P_3_lab[0]+kin3.P_4_lab[0]+kin3.P_5_lab[0]);
  printf("          y-> p3+p4+p5=%f\n",
	 kin3.P_3_lab[1]+kin3.P_4_lab[1]+kin3.P_5_lab[1]);
  printf("          z-> p3+p4+p5=%f\n",
	 kin3.P_3_lab[2]+kin3.P_4_lab[2]+kin3.P_5_lab[2]);

  return;
}

//_____________________________________________________________________________
double
KinemaHybrid::GetEnergy(int i)
{
  switch (i) {
  case 1:
    return kin3.E_1_lab;
    break;
  case 2:
    return kin3.E_2_lab;
    break;
  case 3:
    return kin3.E_3_lab;
    break;
  case 4:
    return kin3.E_4_lab;
    break;
  case 5:
    return kin3.E_5_lab;
    break;
  default:
    fprintf(stderr, "KinemaHybrid::GetEnergy No such particle %d\n", i);
    exit(1);
  }
}

//_____________________________________________________________________________
double
KinemaHybrid::GetMomentum(int i)
{
  switch (i) {
  case 1:
    return kin3.p_1_lab;
    break;
  case 2:
    return kin3.p_2_lab;
    break;
  case 3:
    return kin3.p_3_lab;
    break;
  case 4:
    return kin3.p_4_lab;
    break;
  case 5:
    return kin3.p_5_lab;
    break;
  default:
    fprintf(stderr, "KinemaHybrid::GetMomentum No such particle %d\n", i);
    exit(1);
  }
}

//_____________________________________________________________________________
void
KinemaHybrid::GetMomentum(int i, double *mom)
{
  switch (i) {
  case 1:
    mom[0] = kin3.p_1_lab;
    mom[1] = 0.0;
    mom[2] = 0.0;
    break;
  case 2:
    mom[0] = kin3.p_2_lab;
    mom[1] = 0.0;
    mom[2] = 0.0;
    break;
  case 3:
    mom[0] = kin3.P_3_lab[0];
    mom[1] = kin3.P_3_lab[1];
    mom[2] = kin3.P_3_lab[2];
    break;
  case 4:
    mom[0] = kin3.P_4_lab[0];
    mom[1] = kin3.P_4_lab[1];
    mom[2] = kin3.P_4_lab[2];
    break;
  case 5:
    mom[0] = kin3.P_5_lab[0];
    mom[1] = kin3.P_5_lab[1];
    mom[2] = kin3.P_5_lab[2];
    break;
  default:
    fprintf(stderr, "KinemaHybrid::GetMomentum No such particle %d\n", i);
    exit(1);
  }
}

//_____________________________________________________________________________
double
KinemaHybrid::GetTheta(int i)
{
  switch (i) {
  case 1:
    return 0.0;
    break;
  case 2:
    return 0.0;
    break;
  case 3:
    return kin3.theta3;
    break;
  case 4:
    return kin3.theta4;
    break;
  case 5:
    return kin3.theta5;
    break;
  default:
    fprintf(stderr, "KinemaHybrid::GetTheta No such particle %d\n", i);
    exit(1);
  }
}

//_____________________________________________________________________________
double
KinemaHybrid::GetPhi(int i)
{
  switch (i) {
  case 1:
    return 0.0;
    break;
  case 2:
    return 0.0;
    break;
  case 3:
    return kin3.phi3;
    break;
  case 4:
    return kin3.phi4;
    break;
  case 5:
    return kin3.phi5;
    break;
  default:
    fprintf(stderr, "KinemaHybrid::GetPhi No such particle %d\n", i);
    exit(1);
  }
}

//_____________________________________________________________________________
double
KinemaHybrid::GetThetaCM(int i)
{
  switch (i) {
  case 1:
    return kin3.Theta1CM;
    break;
  case 2:
    return kin3.Theta2CM;
    break;
  default:
    fprintf(stderr, "KinemaHybrid::GetThetaCM index should be 1 or 2 ->%d\n", i);
    exit(1);
  }
}

//_____________________________________________________________________________
double
KinemaHybrid::GetPhiCM(int i)
{
  switch (i) {
  case 1:
    return kin3.Phi1;
    break;
  case 2:
    return kin3.Phi2;
    break;
  default:
    fprintf(stderr, "KinemaHybrid::GetPhiCM index should be 1 or 2 ->%d\n", i);
    exit(1);
  }
}

//_____________________________________________________________________________
void
KinemaHybrid::RotateMom(int i, double deg, double *mom)
{
  double Sin,Cos;

  Sin=sin(deg2rad(deg));
  Cos=cos(deg2rad(deg));
  switch (i) {
  case 3:
    mom[0] = Cos*kin3.P_3_lab[0] - Sin*kin3.P_3_lab[1];
    mom[1] = Sin*kin3.P_3_lab[0] + Cos*kin3.P_3_lab[1];
    mom[2] = kin3.P_3_lab[2];
    break;
  case 4:
    mom[0] = Cos*kin3.P_4_lab[0] - Sin*kin3.P_4_lab[1];
    mom[1] = Sin*kin3.P_4_lab[0] + Cos*kin3.P_4_lab[1];
    mom[2] = kin3.P_4_lab[2];
    break;
  case 5:
    mom[0] = Cos*kin3.P_5_lab[0] - Sin*kin3.P_5_lab[1];
    mom[1] = Sin*kin3.P_5_lab[0] + Cos*kin3.P_5_lab[1];
    mom[2] = kin3.P_5_lab[2];
    break;
  default:
    fprintf(stderr, "KinemaHybrid::RotateMom should be 3,4,5 ->%d\n",i);
    exit(1);
  }
}

//_____________________________________________________________________________
double
KinemaHybrid::GetResMass()
{
  return kin3.M_res;
}

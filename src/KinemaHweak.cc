// -*- C++ -*-

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "KinemaHweak.hh"
#include "Randomize.hh"
#include "G4ios.hh"

//_____________________________________________________________________________
KinemaHweak::KinemaHweak(double m1, double m2, double m3,
                         double m4, double p1, double p2)
{
  double ECM;
  double vx3, vy3, vz3;            /* unit vector */
  double vx4, vy4, vz4;            /* unit vector */
  double theta3, theta4;
  double phi3;                     /* 2phi(CM system)*/
  double Theta3,Phi3;              /*heta,Ph*/
  double Theta4,Phi4;              /* eta,Ph*/

  /*
    G4cout<<"m1 :"<<m1<<G4endl;
    G4cout<<"m2 :"<<m2<<G4endl;
    G4cout<<"m3 :"<<m3<<G4endl;
    G4cout<<"m4 :"<<m4<<G4endl;
    G4cout<<"shhwang test test 1"<<G4endl;
  */

  kin3.M_1 = m1; //--> beam
  kin3.M_2 = m2; //--> target
  kin3.M_3 = m3; //--> K+
  kin3.M_4 = m4; //--> H-dibaryon

  kin3.p_1_lab = p1;
  kin3.p_2_lab = p2;

  kin3.E_1_lab = p2E(kin3.p_1_lab, kin3.M_1);
  kin3.E_2_lab = p2E(kin3.p_2_lab, kin3.M_2);

  ECM = sqrt((kin3.E_1_lab+kin3.E_2_lab)*(kin3.E_1_lab+kin3.E_2_lab)
	     -(kin3.p_1_lab+kin3.p_2_lab)*(kin3.p_1_lab+kin3.p_2_lab));
  //  do {
  //    kin3.M_res = CLHEP::RandBreitWigner::shoot(m_res, width);
  //  } while (kin3.M_3+kin3.M_4 > kin3.M_res || kin3.M_res > ECM-kin3.M_5);
  //  kin3.M_res = m_res;

  if(ECM<m3+m4){
    G4cout<<"Center of energy less than the mass sum for the produced particles"<<G4endl;
    return;
  }

  kin1 = Kinema2Body(m1, m2, m3, m4);
  kin1.SetMomentum(1, p1);
  kin1.SetMomentum(2, p2);
  kin1.SetThetaCM((double)RandSin()); //set to zero
  kin1.CalcKinema();
  phi3 = (-180+360.0*(double)CLHEP::RandFlat::shoot());

  kin3.Theta1CM = kin1.GetThetaCM();
  kin3.Phi1     = phi3;
  /* calculate m4 */
  //  theta4 = kin1.GetThetaLab();
  theta3 = kin1.GetThetaLab();

  vx3 = sin(deg2rad(theta3))*cos(deg2rad(phi3));
  vy3 = sin(deg2rad(theta3))*sin(deg2rad(phi3));
  vz3 = cos(deg2rad(theta3));

  CalcDistoribution(vx3, vy3, vz3, &Theta3, &Phi3);

  kin3.E_3_lab = kin1.GetEnergyLab(3);
  kin3.p_3_lab = kin1.GetMomentumLab(3);

  kin3.P_3_lab[0] = kin3.p_3_lab*vx3;
  kin3.P_3_lab[1] = kin3.p_3_lab*vy3;
  kin3.P_3_lab[2] = kin3.p_3_lab*vz3;

  kin3.theta3 = Theta3;
  kin3.phi3 = Phi3;

  theta4=-kin1.GetPhiLab();

  vx4 = sin(deg2rad(theta4))*cos(deg2rad(phi3));
  vy4 = sin(deg2rad(theta4))*sin(deg2rad(phi3));
  vz4 = cos(deg2rad(theta4));
  CalcDistoribution(vx4, vy4, vz4, &Theta4, &Phi4);
  //  G4cout<<"atan2"<<atan2(vy5,vx5)/PI*180<<G4endl;
  kin3.E_4_lab = kin1.GetEnergyLab(4);
  kin3.p_4_lab = kin1.GetMomentumLab(4);
  kin3.P_4_lab[0] = kin3.p_4_lab*vx4;
  kin3.P_4_lab[1] = kin3.p_4_lab*vy4;
  kin3.P_4_lab[2] = kin3.p_4_lab*vz4;
  kin3.theta4 = Theta4;
  kin3.phi4 = Phi4;
  //Dump();
}

//_____________________________________________________________________________
double
KinemaHweak::p2E(double p,double m)
{
  return sqrt(p*p + m*m);
}

//_____________________________________________________________________________
void
KinemaHweak::CalcDistoribution(double unitx, double unity, double /* unitz */,
                               double *theta, double *phi)
{
  *theta = rag2deg(acos(unitx));
  *phi=rag2deg(atan2(unity,unitx));
  /*  if (unity>=0.0 && unitz>0.0)
   *phi = rag2deg(acos(unity/sin(deg2rad(*theta))));
   else if (unity<0.0 && unitz>=0.0)
   *phi = rag2deg(acos(unity/sin(deg2rad(*theta))));
   else if (unity<=0.0 && unitz<0.0)
   *phi = 360.0-rag2deg(acos(unity/sin(deg2rad(*theta))));
   else if (unity>0.0 && unitz<=0.0)
   *phi = 360.0-rag2deg(acos(unity/sin(deg2rad(*theta))));
   else {
   fprintf(stderr,
   "KinemaHweak::CalcDistribution No such reagion unity=%f, unitz=%f\n",
   unity, unitz);
   Dump();
   exit(1);
   }
  */
  return;
}

//_____________________________________________________________________________
double
KinemaHweak::deg2rad(double theta)
{
  return 3.141592654*theta/180.0;
}

//_____________________________________________________________________________
double
KinemaHweak::rag2deg(double rag)
{
  return 360.0 * rag/ (2.0 * 3.141592654);
}

//_____________________________________________________________________________
double
KinemaHweak::RandSin()
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
KinemaHweak::Dump()
{
  printf("======KinemaHweak Dump======\n");
  printf("--Particle1--\n");
  printf("mass=%f, p_lab=%f, E_lab=%f\n",kin3.M_1, kin3.p_1_lab, kin3.E_1_lab);
  printf("--Particle2--\n");
  printf("mass=%f, p_lab=%f, E_lab=%f\n",kin3.M_2, kin3.p_2_lab, kin3.E_2_lab);
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

  printf("Energy:E1+E2=%f, E3+E4=%f\n", kin3.E_1_lab+kin3.E_2_lab,
	 kin3.E_3_lab+kin3.E_4_lab);
  printf("Momentum: x-> p1+p2=%f, p3+p4=%f\n", kin3.p_1_lab+kin3.p_2_lab,
	 kin3.P_3_lab[0]+kin3.P_4_lab[0]);
  printf("          y-> p3+p4=%f\n",
	 kin3.P_3_lab[1]+kin3.P_4_lab[1]);
  printf("          z-> p3+p4=%f\n",
	 kin3.P_3_lab[2]+kin3.P_4_lab[2]);

  return;
}

//_____________________________________________________________________________
double
KinemaHweak::GetEnergy(int i)
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
  default:
    fprintf(stderr, "KinemaHweak::GetEnergy No such particle %d\n", i);
    exit(1);
  }
}

//_____________________________________________________________________________
double
KinemaHweak::GetMomentum(int i)
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
  default:
    fprintf(stderr, "KinemaHweak::GetMomentum No such particle %d\n", i);
    exit(1);
  }
}

//_____________________________________________________________________________
void
KinemaHweak::GetMomentum(int i, double *mom)
{
  switch (i) {
  case 1:
    mom[0] = 0.0;
    mom[1] = 0.0;
    mom[2] = kin3.p_1_lab;
    break;
  case 2:
    mom[0] = 0.0;
    mom[1] = 0.0;
    mom[2] = kin3.p_2_lab;
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
  default:
    fprintf(stderr, "KinemaHweak::GetMomentum No such particle %d\n", i);
    exit(1);
  }
}

//_____________________________________________________________________________
double
KinemaHweak::GetTheta(int i)
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
  default:
    fprintf(stderr, "KinemaHweak::GetTheta No such particle %d\n", i);
    exit(1);
  }
}

//_____________________________________________________________________________
double
KinemaHweak::GetPhi(int i)
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
  default:
    fprintf(stderr, "KinemaHweak::GetPhi No such particle %d\n", i);
    exit(1);
  }
}

//_____________________________________________________________________________
double
KinemaHweak::GetThetaCM(int i)
{
  switch (i) {
  case 1:
    return kin3.Theta1CM;
    break;
  default:
    fprintf(stderr, "KinemaHweak::GetThetaCM index should be 1 or 2 ->%d\n", i);
    exit(1);
  }
}

//_____________________________________________________________________________
double
KinemaHweak::GetPhiCM(int i)
{
  switch (i) {
  case 1:
    return kin3.Phi1;
    break;
  default:
    fprintf(stderr, "KinemaHweak::GetPhiCM index should be 1 or 2 ->%d\n", i);
    exit(1);
  }
}

//_____________________________________________________________________________
void
KinemaHweak::RotateMom(int i, double deg, double *mom)
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
  default:
    fprintf(stderr, "KinemaHweak::RotateMom should be 3,4->%d\n",i);
    exit(1);
  }
}

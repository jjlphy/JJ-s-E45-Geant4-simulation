// -*- C++ -*-

#include "KinemaFermi.hh"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <G4ios.hh>
#include <Randomize.hh>

#include <CLHEP/Units/SystemOfUnits.h>

#include "FuncName.hh"
#include "PrintHelper.hh"

//_____________________________________________________________________________
KinemaFermi::KinemaFermi(G4double m1, G4double m2, G4double m3, G4double m4,
                         const G4ThreeVector& p1,
                         const G4ThreeVector& p2, G4double cos_theta)
  : m_lv(NumOfParticles)
{
  Calculate(m1, m2, m3, m4, p1, p2, cos_theta);
}

//_____________________________________________________________________________
KinemaFermi::KinemaFermi(G4double m1, G4double m2, G4double m3, G4double m4,
                         G4double *p1, G4double *p2,G4double cos_theta)
  : m_lv(NumOfParticles)
{
  Calculate(m1, m2, m3, m4,
            G4ThreeVector(p1[0], p1[1], p1[2]),
            G4ThreeVector(p2[0], p2[1], p2[2]), cos_theta);
}

//_____________________________________________________________________________
// p1 + p2 -> p3 + p4
void
KinemaFermi::Calculate(G4double m1, G4double m2, G4double m3, G4double m4,
                       const G4ThreeVector& p1,
                       const G4ThreeVector& p2, G4double cos_theta)
{
  m_lv[0].setVectM(p1, m1);
  m_lv[1].setVectM(p2, m2);
  m_lv[2].setVectM(G4ThreeVector(), m3);
  m_lv[3].setVectM(G4ThreeVector(), m4);
  m_lv[4] = m_lv[0] + m_lv[1];

  G4double theta3, theta4;
  G4double phi4;
  G4double phi3;                     /* 2phi(CM system)*/

  ////////////shhwang ///////////////////////////////////
  /////production angle, total angle will be rotate
  G4double theta1 = m_lv[4].theta();
  G4double phi1 = m_lv[4].phi();

  if(m_lv[4].mag() < m3 + m4){
    G4Exception(FUNC_NAME,
                "CM energy less than the total mass.",
                RunMustBeAborted,
                "");
    return;
  }

  m_kinema2body = Kinema2Body(m1, m2, m3, m4);

  G4double theta_p1_psum;
  G4double theta_p2_psum;
  // G4double theta_p1_p2;
  if(m_lv[1].v().mag() == 0){
    theta_p1_psum=0.;
    theta_p2_psum=0.;
    // theta_p1_p2=0.;
  }else {
    theta_p1_psum = m_lv[0].v().theta(m_lv[4]); // beam
    theta_p2_psum = m_lv[1].v().theta(m_lv[4]); // proton
    // theta_p1_p2 = m_lv[0].v().theta(m_lv[1]);
  }

  m_kinema2body.SetMomentum(1, cos(theta_p1_psum)*m_lv[0].v().mag());
  m_kinema2body.SetMomentum(2, cos(theta_p2_psum)*m_lv[1].v().mag());
  m_kinema2body.SetTheta(1, theta_p1_psum);
  m_kinema2body.SetTheta(2, theta_p2_psum);
  m_kinema2body.SetThetaCM(acos(cos_theta)*180./CLHEP::pi);
  m_kinema2body.CalcKinema();

  // v3
  theta3 = m_kinema2body.GetThetaLab();
  phi3 = (360.0*CLHEP::RandFlat::shoot());
  G4ThreeVector v3(
		   (cos(phi1)*cos(deg2rad(theta3))*sin(theta1) +
		    cos(theta1)*cos(phi1)*cos(deg2rad(phi3))*sin(deg2rad(theta3)) -
		    sin(phi1)*sin(deg2rad(phi3))*sin(deg2rad(theta3))),
		   (sin(phi1)*cos(deg2rad(theta3))*sin(theta1) +
		    cos(theta1)*sin(phi1)*cos(deg2rad(phi3))*sin(deg2rad(theta3))+
		    cos(phi1)*sin(deg2rad(phi3))*sin(deg2rad(theta3))),
		   cos(deg2rad(theta3))*cos(theta1) -
		   sin(theta1)*cos(deg2rad(phi3))*sin(deg2rad(theta3)));
  v3.setMag(m_kinema2body.GetMomentumLab(3));
  m_lv[2].setVectM(v3, m_lv[2].m());
  // v4
  theta4 = -m_kinema2body.GetPhiLab();
  phi4 = phi3;
  G4ThreeVector v4(
		   (cos(phi1)*cos(deg2rad(theta4))*sin(theta1) +
		    cos(theta1)*cos(phi1)*cos(deg2rad(phi4))*sin(deg2rad(theta4)) -
		    sin(phi1)*sin(deg2rad(phi4))*sin(deg2rad(theta4))),
		   (sin(phi1)*cos(deg2rad(theta4))*sin(theta1) +
		    cos(theta1)*sin(phi1)*cos(deg2rad(phi4))*sin(deg2rad(theta4)) +
		    cos(phi1)*sin(deg2rad(phi4))*sin(deg2rad(theta4))),
		   cos(deg2rad(theta4))*cos(theta1) -
		   sin(theta1)*cos(deg2rad(phi4))*sin(deg2rad(theta4)));
  v4.setMag(m_kinema2body.GetMomentumLab(4));
  m_lv[3].setVectM(v4, m_lv[3].m());
  m_theta_cm = m_kinema2body.GetThetaCM();
  m_phi_cm = phi3;
}

//_____________________________________________________________________________
KinemaFermi::~KinemaFermi()
{
}

//_____________________________________________________________________________
G4double
KinemaFermi::deg2rad(G4double theta)
{
  return CLHEP::pi*theta/180.0;
}

//_____________________________________________________________________________
G4double
KinemaFermi::rag2deg(G4double rag)
{
  return 360.0 * rag/ (2.0 * CLHEP::pi);
}

//_____________________________________________________________________________
G4double
KinemaFermi::RandSin()
{
  G4int success=0;
  G4double x,fx;

  do {
    x = 180.0 * CLHEP::RandFlat::shoot();
    //x = 180.0 * rand()/(RAND_MAX+1.0);
    fx = sin(deg2rad(x));
    if (fx >= CLHEP::RandFlat::shoot())
      success = 1;
  } while (success==0);

  return x;
}

//_____________________________________________________________________________
G4double
KinemaFermi::GetEnergy(G4int i)
{
  switch (i) {
  case 1:
    return m_lv[0].e();
    break;
  case 2:
    return m_lv[1].e();
    break;
  case 3:
    return m_lv[2].e();
    break;
  case 4:
    return m_lv[3].e();
    break;
  default:
    fprintf(stderr, "KinemaFermi::GetEnergy No such particle %d\n", i);
    exit(1);
  }
}

//_____________________________________________________________________________
const G4LorentzVector&
KinemaFermi::GetLorentzVector(G4int i) const
{
  return m_lv[i];
}

//_____________________________________________________________________________
G4double
KinemaFermi::GetMomentum(G4int i)
{
  switch (i) {
  case 1:
    return m_lv[0].v().mag();
    break;
  case 2:
    return m_lv[1].v().mag();
    break;
  case 3:
    return m_lv[2].v().mag();
    break;
  case 4:
    return m_lv[3].v().mag();
    break;
  default:
    fprintf(stderr, "KinemaFermi::GetMomentum No such particle %d\n", i);
    exit(1);
  }
}

//_____________________________________________________________________________
void
KinemaFermi::GetMomentum(G4int i, G4double *mom)
{
  switch (i) {
  case 1:
    mom[0] = m_lv[0].px();
    mom[1] = m_lv[0].py();
    mom[2] = m_lv[0].pz();
    break;
  case 2:
    mom[0] = m_lv[1].px();
    mom[1] = m_lv[1].py();
    mom[2] = m_lv[1].pz();
    break;
  case 3:
    mom[0] = m_lv[2].px();
    mom[1] = m_lv[2].py();
    mom[2] = m_lv[2].pz();
    break;
  case 4:
    mom[0] = m_lv[3].px();
    mom[1] = m_lv[3].py();
    mom[2] = m_lv[3].pz();
    break;
  default:
    fprintf(stderr, "KinemaFermi::GetMomentum No such particle %d\n", i);
    exit(1);
  }
}

//_____________________________________________________________________________
G4double
KinemaFermi::GetTheta(G4int i)
{
  switch (i) {
  case 1:
    return m_lv[0].theta();
    break;
  case 2:
    return m_lv[1].theta();
    break;
  case 3:
    return m_lv[2].theta();
    break;
  case 4:
    return m_lv[3].theta();
    break;
  default:
    fprintf(stderr, "KinemaFermi::GetTheta No such particle %d\n", i);
    exit(1);
  }
}

//_____________________________________________________________________________
G4double
KinemaFermi::GetPhi(G4int i)
{
  switch (i) {
  case 1:
    return m_lv[0].phi();
    break;
  case 2:
    return m_lv[1].phi();
    break;
  case 3:
    return m_lv[2].phi();
    break;
  case 4:
    return m_lv[3].phi();
    break;
  default:
    fprintf(stderr, "KinemaFermi::GetPhi No such particle %d\n", i);
    exit(1);
  }
}

//_____________________________________________________________________________
G4double
KinemaFermi::GetThetaCM(G4int i)
{
  switch (i) {
  case 1:
    return m_theta_cm;
    break;
  default:
    fprintf(stderr, "KinemaFermi::GetThetaCM index should be 1 or 2 ->%d\n", i);
    exit(1);
  }
}

//_____________________________________________________________________________
G4double
KinemaFermi::GetPhiCM(G4int i)
{
  switch (i) {
  case 1:
    return m_phi_cm;
    break;
  default:
    fprintf(stderr, "KinemaFermi::GetPhiCM index should be 1 or 2 ->%d\n", i);
    exit(1);
  }
}

//_____________________________________________________________________________
void
KinemaFermi::Print() const
{
  PrintHelper helper(7, std::ios::scientific, G4cout);
  static const G4int w = 20;
  G4cout << FUNC_NAME << G4endl
	 << std::setw(w) << "Px"
	 << " " << std::setw(w) << "Py"
	 << " " << std::setw(w) << "Pz"
	 << " " << std::setw(w) << "P"
	 << " " << std::setw(w) << "E"
	 << " " << std::setw(w) << "M"
	 << " " << std::setw(w) << "Theta"
	 << " " << std::setw(w) << "Phi"
	 << G4endl;
  for(const auto& lv : m_lv){
    G4cout << std::setw(w) << lv.px()
	   << " " << std::setw(w) << lv.py()
	   << " " << std::setw(w) << lv.pz()
	   << " " << std::setw(w) << lv.v().mag()
	   << " " << std::setw(w) << lv.e()
	   << " " << std::setw(w) << lv.mag()
	   << " " << std::setw(w) << lv.theta()/CLHEP::deg
	   << " " << std::setw(w) << lv.phi()/CLHEP::deg
	   << G4endl;
  }
}

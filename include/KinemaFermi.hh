// -*- C++ -*-

#ifndef KINEMA_FERMI_HH
#define KINEMA_FERMI_HH

#include <vector>
#include <G4LorentzVector.hh>
#include <G4ThreeVector.hh>

#include "Kinema2Body.hh"

//_____________________________________________________________________________
class KinemaFermi
{
public:
  static G4String ClassName();
  KinemaFermi(G4double m1, G4double m2, G4double m3, G4double m4,
              const G4ThreeVector& p1,
              const G4ThreeVector& p2, G4double cos_theta);
  KinemaFermi(G4double m1, G4double m2, G4double m3, G4double m4,
              G4double *p1, G4double *p2, G4double cos_theta);
  ~KinemaFermi();

private:
  static const G4int NumOfParticles = 2 + 2 + 1; // 2 in, 2 out, and 1 total
  Kinema2Body                  m_kinema2body;
  std::vector<G4LorentzVector> m_lv;
  G4double                     m_theta_cm;
  G4double                     m_phi_cm;

public:
  void                   Calculate(G4double m1, G4double m2,
                                   G4double m3, G4double m4,
                                   const G4ThreeVector& p1,
                                   const G4ThreeVector& p2,
                                   G4double cos_theta);
  G4double               deg2rad(G4double theta);
  G4double               rag2deg(G4double rag);
  G4double               RandSin();
  void                   Dump();
  G4double               GetEnergy(G4int i);
  const G4LorentzVector& GetLorentzVector(G4int i) const;
  G4double               GetMomentum(G4int i);
  void                   GetMomentum(G4int i, G4double *mom);
  G4double               GetTheta(G4int i);
  G4double               GetPhi(G4int i);
  G4double               GetThetaCM(G4int i);
  G4double               GetPhiCM(G4int i);
  void                   Print() const;
};

//_____________________________________________________________________________
inline G4String
KinemaFermi::ClassName()
{
  static G4String s_name("KinemaFermi");
  return s_name;
}

#endif

// -*- C++ -*-

/**
 *  for the angular distribution for the E27 experiment
 */

#include "AngDisGenerator.hh"

#include <cmath>

#include <CLHEP/Units/SystemOfUnits.h>
#include <Randomize.hh>

//_____________________________________________________________________________
AngDisGenerator::AngDisGenerator(G4double cost1, G4double cost2)
  : m_cost1(cost1),
    m_cost2(cost2)
{
}

//_____________________________________________________________________________
AGSWave::AGSWave(G4double cost1, G4double cost2)
  : AngDisGenerator(cost1, cost2)
{
}

//_____________________________________________________________________________
G4ThreeVector
AGSWave::GenerateDirection() const
{
  G4double cost = m_cost1 + G4UniformRand()*(m_cost2 - m_cost1);
  G4double sint = std::sqrt(1. - cost*cost);
  G4double phi  = G4UniformRand()*CLHEP::pi*2.;
  G4double cosp = std::cos(phi);
  G4double sinp = std::sin(phi);
  return G4ThreeVector(sint*cosp, sint*sinp, cost);
}

//_____________________________________________________________________________
AGPWaveFP::AGPWaveFP(G4double cost1, G4double cost2)
  : AngDisGenerator(cost1, cost2)
{
}

//_____________________________________________________________________________
G4ThreeVector
AGPWaveFP::GenerateDirection() const
{
  G4double a    = 2./(m_cost2 - m_cost1)/(2. + m_cost2 + m_cost1);
  G4double cost = -1. + std::sqrt((m_cost1 + 1.)*(m_cost1 + 1.) +
				   2.*G4UniformRand()/a);
  G4double sint = std::sqrt(1. - cost*cost);
  G4double phi  = G4UniformRand()*CLHEP::pi*2.;
  G4double cosp = std::cos(phi);
  G4double sinp = std::sin(phi);
  return G4ThreeVector(sint*cosp, sint*sinp, cost);
}

//_____________________________________________________________________________
AGPWaveBP::AGPWaveBP(G4double cost1, G4double cost2)
  : AngDisGenerator(cost1, cost2)
{
}

//_____________________________________________________________________________
G4ThreeVector
AGPWaveBP::GenerateDirection() const
{
  G4double a    = 2./(m_cost2 - m_cost1)/(2. - m_cost2 - m_cost1);
  G4double cost = 1. - std::sqrt((m_cost1 - 1.)*(m_cost1 - 1.) -
				  2.*G4UniformRand()/a);
  G4double sint = std::sqrt(1. - cost*cost);
  G4double phi  = G4UniformRand()*CLHEP::pi*2.;
  G4double cosp = std::cos(phi);
  G4double sinp = std::sin(phi);
  return G4ThreeVector(sint*cosp, sint*sinp, cost);
}

//_____________________________________________________________________________
AGDWave1::AGDWave1(G4double cost1, G4double cost2)
  : AngDisGenerator(cost1, cost2)
{
}

//_____________________________________________________________________________
G4ThreeVector
AGDWave1::GenerateDirection() const
{
  G4double cost = (m_cost2 - m_cost1)*G4UniformRand() + m_cost1;
  G4double p    = G4UniformRand();
  while (p > Dfunc(cost)){
    cost = (m_cost2 - m_cost1)*G4UniformRand() + m_cost1;
    p    = G4UniformRand();
  }
  G4double sint = std::sqrt(1. - cost*cost);
  G4double phi  = G4UniformRand()*CLHEP::pi*2.;
  G4double cosp = std::cos(phi);
  G4double sinp = std::sin(phi);
  return G4ThreeVector(sint*cosp, sint*sinp, cost);
}

//_____________________________________________________________________________
AGSigma1385Zero::AGSigma1385Zero(G4double cost1, G4double cost2)
  : AngDisGenerator(cost1, cost2)
{
}

//_____________________________________________________________________________
G4ThreeVector
AGSigma1385Zero::GenerateDirection() const
{
  G4double cost = (m_cost2 - m_cost1)*G4UniformRand() + m_cost1;
  G4double p    = G4UniformRand();
  while (p > Dfunc(cost)){
    cost = (m_cost2 - m_cost1)*G4UniformRand() + m_cost1;
    p    = G4UniformRand();
  }
  G4double sint = std::sqrt(1. - cost*cost);
  G4double phi  = G4UniformRand()*CLHEP::pi*2.;
  G4double cosp = std::cos(phi);
  G4double sinp = std::sin(phi);
  return G4ThreeVector(sint*cosp, sint*sinp, cost);
}

//_____________________________________________________________________________
AGSigma1385Plus::AGSigma1385Plus(G4double cost1, G4double cost2)
  : AngDisGenerator(cost1, cost2)
{
}

//_____________________________________________________________________________
G4ThreeVector
AGSigma1385Plus::GenerateDirection() const
{
  G4double cost = (m_cost2 - m_cost1)*G4UniformRand() + m_cost1;
  G4double p    = G4UniformRand();
  while(p > Dfunc(cost)){
    cost = (m_cost2 - m_cost1)*G4UniformRand() + m_cost1;
    p    = G4UniformRand();
  }
  G4double sint = std::sqrt(1. - cost*cost);
  G4double phi  = G4UniformRand()*CLHEP::pi*2.;
  G4double cosp = std::cos(phi);
  G4double sinp = std::sin(phi);
  return G4ThreeVector(sint*cosp, sint*sinp, cost);
}

//_____________________________________________________________________________
AGLambda1405::AGLambda1405(G4double cost1, G4double cost2)
  : AngDisGenerator(cost1, cost2)
{
}

//_____________________________________________________________________________
G4ThreeVector
AGLambda1405::GenerateDirection() const
{
  G4double cost = (m_cost2 - m_cost1)*G4UniformRand() + m_cost1;
  G4double p    = G4UniformRand();
  while (p > Dfunc(cost)){
    cost = (m_cost2 - m_cost1)*G4UniformRand() + m_cost1;
    p    = G4UniformRand();
  }
  G4double sint = std::sqrt(1. - cost*cost);
  G4double phi  = G4UniformRand()*CLHEP::pi*2.;
  G4double cosp = std::cos(phi);
  G4double sinp = std::sin(phi);
  return G4ThreeVector(sint*cosp, sint*sinp, cost);
}

//_____________________________________________________________________________
AGLambda::AGLambda(G4double cost1, G4double cost2)
  : AngDisGenerator(cost1, cost2)
{
}

//_____________________________________________________________________________
G4ThreeVector
AGLambda::GenerateDirection() const
{
  G4double cost = (m_cost2 - m_cost1)*G4UniformRand() + m_cost1;
  G4double p    = G4UniformRand();
  while(p > Dfunc(cost)){
    cost = (m_cost2 - m_cost1)*G4UniformRand() + m_cost1;
    p    = G4UniformRand();
  }
  G4double sint = std::sqrt(1. - cost*cost);
  G4double phi  = G4UniformRand()*CLHEP::pi*2.;
  G4double cosp = std::cos(phi);
  G4double sinp = std::sin(phi);
  return G4ThreeVector(sint*cosp, sint*sinp, cost);
}

//_____________________________________________________________________________
AGSigmaZ::AGSigmaZ(G4double cost1, G4double cost2)
  : AngDisGenerator(cost1, cost2)
{
}

//_____________________________________________________________________________
G4ThreeVector
AGSigmaZ::GenerateDirection() const
{
  G4double cost = (m_cost2 - m_cost1)*G4UniformRand() + m_cost1;
  G4double p    = G4UniformRand();
  while (p > Dfunc(cost)){
    cost = (m_cost2 - m_cost1)*G4UniformRand() + m_cost1;
    p    = G4UniformRand();
  }
  G4double sint = std::sqrt(1. - cost*cost);
  G4double phi  = G4UniformRand()*CLHEP::pi*2.;
  G4double cosp = std::cos(phi);
  G4double sinp = std::sin(phi);
  return G4ThreeVector(sint*cosp, sint*sinp, cost);
}

//_____________________________________________________________________________
AGSigmaP::AGSigmaP(G4double cost1, G4double cost2)
  : AngDisGenerator(cost1, cost2)
{
}

//_____________________________________________________________________________
G4ThreeVector
AGSigmaP::GenerateDirection() const
{
  G4double cost = (m_cost2 - m_cost1)*G4UniformRand() + m_cost1;
  G4double p    = G4UniformRand();
  while(p > Dfunc(cost)){
    cost = (m_cost2 - m_cost1)*G4UniformRand() + m_cost1;
    p    = G4UniformRand();
  }
  G4double sint = std::sqrt(1. - cost*cost);
  G4double phi  = G4UniformRand()*CLHEP::pi*2.;
  G4double cosp = std::cos(phi);
  G4double sinp = std::sin(phi);
  return G4ThreeVector(sint*cosp, sint*sinp, cost);
}

//_____________________________________________________________________________
AGPol::AGPol(G4double cost1, G4double cost2)
  : AngDisGenerator(cost1, cost2)
{
}

//_____________________________________________________________________________
G4ThreeVector
AGPol::GenerateDirection() const
{
  G4double a    = 2./(m_cost2 - m_cost1)/(2. + m_cost2 + m_cost1);
  G4double cost = -1. + std::sqrt((m_cost1 + 1.)*(m_cost1 + 1.) +
				   2.*G4UniformRand()/a);
  G4double sint = std::sqrt(1. - cost*cost);
  G4double phi  = G4UniformRand()*CLHEP::pi*2.;
  G4double cosp = std::cos(phi);
  G4double sinp = std::sin(phi);
  return G4ThreeVector(sint*cosp, sint*sinp, cost);
}

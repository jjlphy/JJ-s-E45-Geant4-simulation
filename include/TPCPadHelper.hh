// -*- C++ -*-

#ifndef TPCPADHELPER_HH
#define TPCPADHELPER_HH

#include <G4String.hh>
#include <G4ThreeVector.hh>
#include <CLHEP/Units/PhysicalConstants.h>
#include <CLHEP/Units/SystemOfUnits.h>




namespace TPCPadHelper
{
G4int GetPadID(const G4int layerID, const G4int rowID);
G4int GetLayerID(G4int padID);
G4int GetRowID(G4int padID);
G4double GetTheta(G4int padID);
G4double GetTheta(const G4int layerID, const G4double m_row);
G4double GetMrow(const G4int layerID, const G4double m_phi);
G4double GetRadius(const G4int layerID);
G4double GetR(G4int padID);
G4ThreeVector GetPosition(G4int padID);
G4ThreeVector GetPosition(const G4int layerID, const G4double m_row);
G4int FindPadID(G4double z, G4double x);
G4double ArcLength(const G4int layerID, const G4double row1, const G4double row2);
G4bool GetDeadCon(const G4int padID);
G4bool GetDeadCon(const G4int layerID, const G4int rowID);

inline G4String ClassName()
{
  static G4String s_name("TPCPadHelper");
  return s_name;
}
}

#endif

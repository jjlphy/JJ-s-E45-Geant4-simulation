// -*- C++ -*-

#ifndef BH2_SD_HH
#define BH2_SD_HH

#include <G4VSensitiveDetector.hh>

#include "BH2Hit.hh"

class G4HCofThisEvent;
class G4Step;
class G4TouchableHistory;

//_____________________________________________________________________________
class BH2SD : public G4VSensitiveDetector
{
public:
  static G4String ClassName();
  BH2SD(const G4String& name);
  virtual ~BH2SD();

private:
  G4THitsCollection<BH2Hit>* m_hits_collection;

public:
  virtual G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist);
  virtual void   Initialize(G4HCofThisEvent* HCTE);
  virtual void   EndOfEvent(G4HCofThisEvent* HCTE);
  virtual void   DrawAll();
  virtual void   PrintAll();
};

//_____________________________________________________________________________
inline G4String
BH2SD::ClassName()
{
  static G4String s_name("BH2SD");
  return s_name;
}

#endif

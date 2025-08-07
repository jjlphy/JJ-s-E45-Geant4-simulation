// -*- C++ -*-

#ifndef VP_SD_HH
#define VP_SD_HH

#include <G4VSensitiveDetector.hh>

#include "VPHit.hh"

class G4HCofThisEvent;
class G4Step;
class G4TouchableHistory;

//_____________________________________________________________________________
class VPSD : public G4VSensitiveDetector
{
public:
  static G4String ClassName();
  VPSD(const G4String& name);
  virtual ~VPSD();

private:
  G4THitsCollection<VPHit>* m_hits_collection;

public:
  virtual G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist);
  virtual void   Initialize(G4HCofThisEvent* HCTE);
  virtual void   EndOfEvent(G4HCofThisEvent* HCTE);
  virtual void   DrawAll();
  virtual void   PrintAll();
};

//_____________________________________________________________________________
inline G4String
VPSD::ClassName()
{
  static G4String s_name("VPSD");
  return s_name;
}

#endif

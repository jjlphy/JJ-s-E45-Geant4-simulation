// -*- C++ -*-

#ifndef SDC_SD_HH
#define SDC_SD_HH

#include <G4VSensitiveDetector.hh>

#include "SDCHit.hh"

class G4HCofThisEvent;
class G4Step;
class G4TouchableHistory;

//_____________________________________________________________________________
class SDCSD : public G4VSensitiveDetector
{
public:
  static G4String ClassName();
  SDCSD(const G4String& name);
  virtual ~SDCSD();

private:
  G4THitsCollection<SDCHit>* m_hits_collection;

public:
  virtual G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist);
  virtual void   Initialize(G4HCofThisEvent* HCTE);
  virtual void   EndOfEvent(G4HCofThisEvent* HCTE);
  virtual void   DrawAll();
  virtual void   PrintAll();
};

//_____________________________________________________________________________
inline G4String
SDCSD::ClassName()
{
  static G4String s_name("SDCSD");
  return s_name;
}

#endif

// -*- C++ -*-

#ifndef SCH_SD_HH
#define SCH_SD_HH

#include <G4VSensitiveDetector.hh>

#include "SCHHit.hh"

class G4HCofThisEvent;
class G4Step;
class G4TouchableHistory;

//_____________________________________________________________________________
class SCHSD : public G4VSensitiveDetector
{
public:
  static G4String ClassName();
  SCHSD(const G4String& name);
  virtual ~SCHSD();

private:
  G4THitsCollection<SCHHit>* m_hits_collection;

public:
  virtual G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist);
  virtual void   Initialize(G4HCofThisEvent* HCTE);
  virtual void   EndOfEvent(G4HCofThisEvent* HCTE);
  virtual void   DrawAll();
  virtual void   PrintAll();
};

//_____________________________________________________________________________
inline G4String
SCHSD::ClassName()
{
  static G4String s_name("SCHSD");
  return s_name;
}

#endif

// -*- C++ -*-

#ifndef TARGET_SD_HH
#define TARGET_SD_HH

#include <G4VSensitiveDetector.hh>

#include "TargetHit.hh"

class G4HCofThisEvent;
class G4Step;
class G4TouchableHistory;

//_____________________________________________________________________________
class TargetSD : public G4VSensitiveDetector
{
public:
  static G4String ClassName();
  TargetSD(const G4String& name);
  virtual ~TargetSD();

private:
  G4THitsCollection<TargetHit>* m_hits_collection;

public:
  G4int ntrk;
  virtual G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist);
  virtual void   Initialize(G4HCofThisEvent* HCTE);
  virtual void   EndOfEvent(G4HCofThisEvent* HCTE);
  virtual void   DrawAll();
  virtual void   PrintAll();

};

//_____________________________________________________________________________
inline G4String
TargetSD::ClassName()
{
  static G4String s_name("TargetSD");
  return s_name;
}

#endif

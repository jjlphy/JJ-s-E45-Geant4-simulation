// -*- C++ -*-

#ifndef HTOF_SD_HH
#define HTOF_SD_HH

#include <G4VSensitiveDetector.hh>

#include "HTOFHit.hh"

class G4HCofThisEvent;
class G4Step;
class G4TouchableHistory;

//_____________________________________________________________________________
class HTOFSD : public G4VSensitiveDetector
{
public:
  static G4String ClassName();
  HTOFSD(const G4String& name);
  virtual ~HTOFSD();

private:
  G4THitsCollection<HTOFHit>* m_hits_collection;

public:
  virtual G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist);
  virtual void   Initialize(G4HCofThisEvent* HCTE);
  virtual void   EndOfEvent(G4HCofThisEvent* HCTE);
  virtual void   DrawAll();
  virtual void   PrintAll();
};

//_____________________________________________________________________________
inline G4String
HTOFSD::ClassName()
{
  static G4String s_name("HTOFSD");
  return s_name;
}

#endif

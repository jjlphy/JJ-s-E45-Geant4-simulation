// -*- C++ -*-

#ifndef FTOF_SD_HH
#define FTOF_SD_HH

#include <G4VSensitiveDetector.hh>

#include "FTOFHit.hh"

class G4HCofThisEvent;
class G4Step;
class G4TouchableHistory;

//_____________________________________________________________________________
class FTOFSD : public G4VSensitiveDetector
{
public:
  static G4String ClassName();
  FTOFSD(const G4String& name);
  virtual ~FTOFSD();

private:
  G4THitsCollection<FTOFHit>* m_hits_collection;

public:
  virtual G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist);
  virtual void   Initialize(G4HCofThisEvent* HCTE);
  virtual void   EndOfEvent(G4HCofThisEvent* HCTE);
  virtual void   DrawAll();
  virtual void   PrintAll();
};

//_____________________________________________________________________________
inline G4String
FTOFSD::ClassName()
{
  static G4String s_name("FTOFSD");
  return s_name;
}

#endif

// -*- C++ -*-

#ifndef TPC_SD_HH
#define TPC_SD_HH

#include <G4VSensitiveDetector.hh>

#include "TPCHit.hh"

class G4HCofThisEvent;
class G4Step;
class G4TouchableHistory;

//_____________________________________________________________________________
class TPCSD : public G4VSensitiveDetector
{
public:
  static G4String ClassName();
  TPCSD(const G4String& name);
  virtual ~TPCSD();

private:
  G4THitsCollection<TPCHit>* m_hits_collection;

public:
  virtual G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist);
  virtual void   Initialize(G4HCofThisEvent* HCTE);
  virtual void   EndOfEvent(G4HCofThisEvent* HCTE);
  virtual void   DrawAll();
  virtual void   PrintAll();
};

//_____________________________________________________________________________
inline G4String
TPCSD::ClassName()
{
  static G4String s_name("TPCSD");
  return s_name;
}
#endif

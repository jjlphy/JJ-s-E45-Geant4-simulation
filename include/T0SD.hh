// -*- C++ -*-
#ifndef T0_SD_HH
#define T0_SD_HH

#include <G4VSensitiveDetector.hh>
#include <G4THitsCollection.hh>

#include "T0Hit.hh"

class G4Step;
class G4TouchableHistory;

class T0SD : public G4VSensitiveDetector
{
public:
  T0SD(const G4String& name);
  virtual ~T0SD();

  virtual void Initialize(G4HCofThisEvent* HCTE);
  virtual G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist);
  virtual void EndOfEvent(G4HCofThisEvent* HCTE);
  virtual void DrawAll();
  virtual void PrintAll();

private:
  G4THitsCollection<T0Hit>* m_hits_collection;
};

#endif

// -*- C++ -*-

#ifndef LAC_SD_HH
#define LAC_SD_HH

#include <G4VSensitiveDetector.hh>

#include "LACHit.hh"

class G4HCofThisEvent;
class G4Step;
class G4TouchableHistory;

//_____________________________________________________________________________
class LACSD : public G4VSensitiveDetector
{
public:
  static G4String ClassName();
  LACSD(const G4String& name);
  virtual ~LACSD();

private:
  G4THitsCollection<LACHit>* m_hits_collection;
  G4double                   m_refractive_index;

public:
  G4double GetRefractiveIndex() const { return m_refractive_index; }
  void     SetRefractiveIndex(G4double index){ m_refractive_index = index; }

public:
  virtual G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist);
  virtual void   Initialize(G4HCofThisEvent* HCTE);
  virtual void   EndOfEvent(G4HCofThisEvent* HCTE);
  virtual void   DrawAll();
  virtual void   PrintAll();
};

//_____________________________________________________________________________
inline G4String
LACSD::ClassName()
{
  static G4String s_name("LACSD");
  return s_name;
}

#endif

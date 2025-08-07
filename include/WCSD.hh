// -*- C++ -*-

#ifndef WC_SD_HH
#define WC_SD_HH

#include <G4VSensitiveDetector.hh>

#include "WCHit.hh"

class G4HCofThisEvent;
class G4Step;
class G4TouchableHistory;

//_____________________________________________________________________________
class WCSD : public G4VSensitiveDetector
{
public:
  static G4String ClassName();
  WCSD(const G4String& name);
  virtual ~WCSD();

private:
  G4THitsCollection<WCHit>* m_hits_collection;
  G4double                  m_refractive_index;

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
WCSD::ClassName()
{
  static G4String s_name("WCSD");
  return s_name;
}

#endif

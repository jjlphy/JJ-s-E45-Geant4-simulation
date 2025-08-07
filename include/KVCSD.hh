// -*- C++ -*-

#ifndef KVC_SD_HH
#define KVC_SD_HH

#include <G4VSensitiveDetector.hh>

#include "KVCHit.hh"

class G4HCofThisEvent;
class G4Step;
class G4TouchableHistory;

//_____________________________________________________________________________
class KVCSD : public G4VSensitiveDetector
{
public:
  static G4String ClassName();
  KVCSD(const G4String& name);
  virtual ~KVCSD();

private:
  G4THitsCollection<KVCHit>* m_hits_collection;
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
KVCSD::ClassName()
{
  static G4String s_name("KVCSD");
  return s_name;
}

#endif

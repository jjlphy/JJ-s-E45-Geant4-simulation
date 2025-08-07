#ifndef BVH_U_SD_HH
#define BVH_U_SD_HH

#include <G4VSensitiveDetector.hh>
#include <G4Step.hh>
#include <G4THitsCollection.hh>
#include <G4ThreeVector.hh>

#include "BVH_UHit.hh"

class BVH_U_SD : public G4VSensitiveDetector
{
public:
  BVH_U_SD(const G4String& name);
  virtual ~BVH_U_SD() = default;

  virtual void Initialize(G4HCofThisEvent* HCE) override;
  virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory*) override;

private:
  G4THitsCollection<BVH_UHit>* m_hits;
  G4int m_hitID;
};

#endif


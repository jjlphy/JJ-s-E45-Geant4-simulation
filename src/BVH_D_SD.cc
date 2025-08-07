#include "BVH_D_SD.hh"
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"

BVH_D_SD::BVH_D_SD(const G4String& name)
  : G4VSensitiveDetector(name), m_hits(nullptr), m_hitID(-1)
{
  collectionName.insert("hit");
}

void BVH_D_SD::Initialize(G4HCofThisEvent* HCE)
{
  m_hits = new G4THitsCollection<BVH_DHit>(SensitiveDetectorName, collectionName[0]);
  m_hitID = GetCollectionID(0);
  HCE->AddHitsCollection(m_hitID, m_hits);
}

G4bool BVH_D_SD::ProcessHits(G4Step* step, G4TouchableHistory*)
{
  auto hit = new BVH_DHit(SensitiveDetectorName, step);
  m_hits->insert(hit);
  return true;
}

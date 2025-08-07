#include "BVH_U_SD.hh"
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"

BVH_U_SD::BVH_U_SD(const G4String& name)
  : G4VSensitiveDetector(name), m_hits(nullptr), m_hitID(-1)
{
  collectionName.insert("hit");
}

void BVH_U_SD::Initialize(G4HCofThisEvent* HCE)
{
  m_hits = new G4THitsCollection<BVH_UHit>(SensitiveDetectorName, collectionName[0]);
  m_hitID = GetCollectionID(0);
  HCE->AddHitsCollection(m_hitID, m_hits);
}

G4bool BVH_U_SD::ProcessHits(G4Step* step, G4TouchableHistory*)
{
  auto hit = new BVH_UHit(SensitiveDetectorName, step);
  m_hits->insert(hit);
  return true;
}

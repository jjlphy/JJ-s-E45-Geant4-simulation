#include "T0SD.hh"

#include <G4Step.hh>
#include <G4TouchableHistory.hh>
#include <G4Track.hh>
#include <G4VPhysicalVolume.hh>

#include "T0Hit.hh"

T0SD::T0SD(const G4String& name)
  : G4VSensitiveDetector(name), m_hits_collection()
{
  collectionName.insert("hit");
}

T0SD::~T0SD() {}

void T0SD::Initialize(G4HCofThisEvent* HCTE)
{
  m_hits_collection = new G4THitsCollection<T0Hit>(SensitiveDetectorName, collectionName[0]);
  HCTE->AddHitsCollection(GetCollectionID(0), m_hits_collection);
}

G4bool T0SD::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
  const auto pre = aStep->GetPreStepPoint();
  const auto def = aStep->GetTrack()->GetDefinition();

  if (pre->GetStepStatus() != fGeomBoundary) return false;
  if (def->GetPDGCharge() == 0.) return false;

  m_hits_collection->insert(new T0Hit(SensitiveDetectorName, aStep));
  return true;
}

void T0SD::EndOfEvent(G4HCofThisEvent*) {}
void T0SD::DrawAll() {}
void T0SD::PrintAll() { m_hits_collection->PrintAllHits(); }

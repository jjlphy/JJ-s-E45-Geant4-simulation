// -*- C++ -*-

#include "HTOFSD.hh"

#include <G4Step.hh>
#include <G4TouchableHistory.hh>
#include <G4Track.hh>
#include <G4VPhysicalVolume.hh>
#include <G4VTouchable.hh>

#include "FuncName.hh"
#include "HTOFHit.hh"

//_____________________________________________________________________________
HTOFSD::HTOFSD(const G4String& name)
  : G4VSensitiveDetector(name),
    m_hits_collection()
{
  collectionName.insert("hit");
}

//_____________________________________________________________________________
HTOFSD::~HTOFSD()
{
}

//_____________________________________________________________________________
void
HTOFSD::Initialize(G4HCofThisEvent* HCTE)
{
  m_hits_collection = new G4THitsCollection<HTOFHit>(SensitiveDetectorName,
                                                     collectionName[0]);
  HCTE->AddHitsCollection(GetCollectionID(0), m_hits_collection);
}

//_____________________________________________________________________________
G4bool
HTOFSD::ProcessHits(G4Step* aStep, G4TouchableHistory* /* ROhist */)
{
  const auto preStepPoint = aStep->GetPreStepPoint();
  const auto aTrack = aStep->GetTrack();
  const auto Definition = aTrack->GetDefinition();
  const G4String particleName = Definition->GetParticleName();
  const G4String particleType = Definition->GetParticleType();

  if(preStepPoint->GetStepStatus() != fGeomBoundary)
    return false;
  if(Definition->GetPDGCharge() == 0.)
    return false;

  // if(particleName == "e-")
  //   return false;
  // if(particleName == "e+")
  //   return false;
  // if(particleName != "kaon+")
  //   return false;
  // if(particleName != "pi-" && particleName != "pi+")
  //   return false;
  // if(particleName != "pi+" && particleName != "pi-" &&
  //     particleName != "proton")
  //   return false;
  // if(particleType == "lepton")
  //   return false;

  m_hits_collection->insert(new HTOFHit(SensitiveDetectorName, aStep));

  return true;
}

//_____________________________________________________________________________
void
HTOFSD::EndOfEvent(G4HCofThisEvent* /* HCTE */)
{
}

//_____________________________________________________________________________
void
HTOFSD::DrawAll()
{
}

//_____________________________________________________________________________
void
HTOFSD::PrintAll()
{
  m_hits_collection->PrintAllHits();
}

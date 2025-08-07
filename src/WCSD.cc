// -*- C++ -*-

#include "WCSD.hh"

#include <G4Step.hh>
#include <G4TouchableHistory.hh>
#include <G4Track.hh>
#include <G4VPhysicalVolume.hh>
#include <G4VTouchable.hh>

#include "FuncName.hh"
#include "WCHit.hh"

//_____________________________________________________________________________
WCSD::WCSD(const G4String& name)
  : G4VSensitiveDetector(name),
    m_hits_collection(),
    m_refractive_index()
{
  collectionName.insert("hit");
}

//_____________________________________________________________________________
WCSD::~WCSD()
{
}

//_____________________________________________________________________________
void
WCSD::Initialize(G4HCofThisEvent* HCTE)
{
  m_hits_collection = new G4THitsCollection<WCHit>(SensitiveDetectorName,
                                                   collectionName[0]);
  HCTE->AddHitsCollection(GetCollectionID(0), m_hits_collection);
}

//_____________________________________________________________________________
G4bool
WCSD::ProcessHits(G4Step* aStep, G4TouchableHistory* /* ROhist */)
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

  // const G4double momentum_threshold =
  //   hit->GetMass() / std::sqrt(m_refractive_index*m_refractive_index - 1.);
  // if(hit->GetMomentum().mag() < momentum_threshold)
  //   return false;

  // if(hit->GetParticleName().contains("proton")){
  //   hit->Print();
  //   G4cout << momentum_threshold << G4endl;
  // }

  m_hits_collection->insert(new WCHit(SensitiveDetectorName, aStep));

  return true;
}

//_____________________________________________________________________________
void
WCSD::EndOfEvent(G4HCofThisEvent* /* HCTE */)
{
}

//_____________________________________________________________________________
void
WCSD::DrawAll()
{
}

//_____________________________________________________________________________
void
WCSD::PrintAll()
{
  m_hits_collection->PrintAllHits();
}

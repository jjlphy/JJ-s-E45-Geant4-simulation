// -*- C++ -*-

#include "TPCSD.hh"

#include <G4VPhysicalVolume.hh>
#include <G4Step.hh>
#include <G4Track.hh>
#include <G4VTouchable.hh>
#include <G4TouchableHistory.hh>
#include <G4PhysicalConstants.hh>
#include <G4SystemOfUnits.hh>
#include <G4DynamicParticle.hh>
#include <G4DecayProducts.hh>
#include <G4PhysicsLogVector.hh>
#include <G4ParticleChangeForDecay.hh>
#include <G4DecayProcessType.hh>
#include <Randomize.hh>

#include "ConfMan.hh"
#include "FuncName.hh"
#include "TPCHit.hh"
#include "padHelper.hh"

namespace
{
const auto& gConf = ConfMan::GetInstance();
}

//_____________________________________________________________________________
TPCSD::TPCSD(const G4String& name)
  : G4VSensitiveDetector(name),
    m_hits_collection()
{
  collectionName.insert("hit");
}

//_____________________________________________________________________________
TPCSD::~TPCSD()
{
}

//_____________________________________________________________________________
void
TPCSD::Initialize(G4HCofThisEvent* HCTE)
{
  m_hits_collection = new G4THitsCollection<TPCHit>(SensitiveDetectorName,
						    collectionName[0]);
  HCTE->AddHitsCollection(GetCollectionID(0), m_hits_collection);
}

//_____________________________________________________________________________
G4bool
TPCSD::ProcessHits(G4Step* aStep, G4TouchableHistory* /* ROhist */)
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

  
  

  m_hits_collection->insert(new TPCHit(SensitiveDetectorName, aStep));

  return true;
}

//_____________________________________________________________________________
void
TPCSD::EndOfEvent(G4HCofThisEvent* /* HCTE */)
{
}

//_____________________________________________________________________________
void
TPCSD::DrawAll()
{
}

//_____________________________________________________________________________
void
TPCSD::PrintAll()
{
  m_hits_collection->PrintAllHits();
}

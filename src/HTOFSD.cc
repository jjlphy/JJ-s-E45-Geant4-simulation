// -*- C++ -*-

#include "HTOFSD.hh"

#include <G4Step.hh>
#include <G4TouchableHistory.hh>
#include <G4Track.hh>
#include <G4VPhysicalVolume.hh>
#include <G4VTouchable.hh>

#include "FuncName.hh"
#include "HTOFHit.hh"
#include "TParticle.h"  // TParticle::SetStatusCode 사용

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

  // 히트 생성
  auto hit = new HTOFHit(SensitiveDetectorName, aStep);

  // ★ copy_no를 TParticle의 StatusCode에 저장
  if (const auto touch = preStepPoint->GetTouchable()) {
    const int copyNo = touch->GetCopyNumber(0);  // 현재 PV의 copy number
    if (auto tp = hit->GetParticle()) {          // GetParticle()이 TParticle* 반환
      tp->SetStatusCode(copyNo);
    }
  }

  m_hits_collection->insert(hit);
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

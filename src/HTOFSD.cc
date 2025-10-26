// -*- C++ -*-
#include "HTOFSD.hh"

#include <G4Step.hh>
#include <G4Track.hh>
#include <G4VTouchable.hh>
#include <G4VPhysicalVolume.hh>
#include <G4TouchableHistory.hh>
#include <G4ParticleDefinition.hh>

#include "HTOFHit.hh"
#include "FuncName.hh"
#include "TParticle.h"

// ★ 추가: OriginTag 읽기 위해
#include "TrackTag.hh"  // EOrigin, OriginTag

// prefix 체크 헬퍼는 그대로…
static inline bool StartsWith(const G4String& s, const char* prefix){
  return s.find(prefix) == 0;
}

HTOFSD::HTOFSD(const G4String& name)
  : G4VSensitiveDetector(name),
    m_hits_collection()
{
  collectionName.insert("hit");
}

HTOFSD::~HTOFSD() = default;

void HTOFSD::Initialize(G4HCofThisEvent* HCTE)
{
  m_hits_collection =
      new G4THitsCollection<HTOFHit>(SensitiveDetectorName, collectionName[0]);
  HCTE->AddHitsCollection(GetCollectionID(0), m_hits_collection);
}

G4bool HTOFSD::ProcessHits(G4Step* aStep, G4TouchableHistory* /* ROhist */)
{
  const auto pre = aStep->GetPreStepPoint();
  const auto trk = aStep->GetTrack();
  const auto def = trk->GetDefinition();

  if (!def || def->GetPDGCharge() == 0.) return false;

  const auto stat = pre->GetStepStatus();
  const bool accept =
      aStep->IsFirstStepInVolume() ||
      (stat == fGeomBoundary) ||
      (stat == fUndefined);
  if (!accept) return false;

  // --- copy-no 추출 ---
  int copyNo = -1;
  G4String pvName = "";
  if (const auto touch = pre->GetTouchable()){
    if (const auto pv = touch->GetVolume()){
      pvName = pv->GetName();
      copyNo = pv->GetCopyNo();
    }
  }
  if (copyNo < 0 || !StartsWith(pvName, "HtofPV")){
    if (const auto post = aStep->GetPostStepPoint()){
      if (const auto touch2 = post->GetTouchable()){
        if (const auto pv2 = touch2->GetVolume()){
          const G4String nm2 = pv2->GetName();
          const int      cp2 = pv2->GetCopyNo();
          if (cp2 >= 0 && StartsWith(nm2, "HtofPV")){
            copyNo = cp2;
            pvName = nm2;
          }
        }
      }
    }
  }
  if (copyNo < 0 || !StartsWith(pvName, "HtofPV")) return false;

  // --- 히트 생성 ---
  auto hit = new HTOFHit(SensitiveDetectorName, aStep);

  // ★★★ 핵심: copy-no는 StatusCode에, Origin은 UniqueID에 저장 ★★★
  if (auto tp = hit->GetParticle()){
    // 1) 타일 ID로 쓰기 위해 copy-no 유지
    tp->SetStatusCode(copyNo);

    // 2) 출처(origin) 태그를 UniqueID로 함께 보냄
    int origin_code = -1; // unknown
    if (auto ui = trk->GetUserInformation()){
      if (auto ot = dynamic_cast<const OriginTag*>(ui)){
        // EOrigin enum을 정수로
        origin_code = static_cast<int>(ot->origin);
        // 예: kPrimaryBeam=0, kForced2PiChild=1 … (TrackTag.hh 정의에 따름)
      }
    }
    tp->SetUniqueID(origin_code);
  }

  m_hits_collection->insert(hit);
  return true;
}

void HTOFSD::EndOfEvent(G4HCofThisEvent*) {}
void HTOFSD::DrawAll() {}
void HTOFSD::PrintAll(){ if (m_hits_collection) m_hits_collection->PrintAllHits(); }

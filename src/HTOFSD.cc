// -*- C++ -*-
// HTOFSD.cc
// 2025-10-26 (user modification by jaejin & assistant)
// - TParticle::StatusCode  : HTOF tile ID (copy-no) 유지
// - TParticle::UniqueID    : OriginTag(EOrigin) 값을 정수로 저장
// - TParticle::PdgCode     : 방어용으로 PDG 코드도 함께 저장
// - OriginTag 멤버명은 origin 사용 (TrackTag.hh 기준)
// - SimpleTag("UPSTREAM_HELPER")를 kPrimaryBeam으로 폴백 태깅

#include "HTOFSD.hh"

#include <G4Step.hh>
#include <G4Track.hh>
#include <G4VTouchable.hh>
#include <G4VPhysicalVolume.hh>
#include <G4TouchableHistory.hh>
#include <G4ParticleDefinition.hh>
#include <G4StepStatus.hh>

#include "HTOFHit.hh"
#include "FuncName.hh"
#include "TParticle.h"

// ★ OriginTag / SimpleTag
#include "TrackTag.hh"  // EOrigin, OriginTag, SimpleTag

// 작은 헬퍼: 문자열 prefix 체크
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

  // 중성 제외 (전하 0)
  if (!def || def->GetPDGCharge() == 0.) return false;

  // 엔트리 조건: 첫스텝/경계/Undefined
  const auto stat = pre->GetStepStatus();
  const bool accept =
      aStep->IsFirstStepInVolume() ||
      (stat == fGeomBoundary) ||
      (stat == fUndefined);
  if (!accept) return false;

  // --- copy-no 추출 (pre 우선, post 폴백) ---
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

  // --- TParticle에 메타데이터 주입 ---
  if (auto tp = hit->GetParticle()){
    // (1) 타일 ID = copy-no
    tp->SetStatusCode(copyNo);

    // (2) OriginTag → UniqueID (정수화)
    int origin_code = -1; // unknown
    if (auto ui = trk->GetUserInformation()){
      if (auto ot = dynamic_cast<const OriginTag*>(ui)){
        // TrackTag.hh 에서 멤버명이 origin 임
        origin_code = static_cast<int>(ot->origin);
      } else if (auto st = dynamic_cast<const SimpleTag*>(ui)){
        // 하위호환: SimpleTag("UPSTREAM_HELPER") → kPrimaryBeam 간주
        if (st->why_ && std::string(st->why_)=="UPSTREAM_HELPER"){
          origin_code = static_cast<int>(EOrigin::kPrimaryBeam);
        }
      }
    }
    tp->SetUniqueID(origin_code);

    // (3) PDG 코드도 기록(오프라인 방어/진단에 유용)
    tp->SetPdgCode(def->GetPDGEncoding());
  }

  m_hits_collection->insert(hit);
  return true;
}

void HTOFSD::EndOfEvent(G4HCofThisEvent* /* HCTE */) {}
void HTOFSD::DrawAll() {}
void HTOFSD::PrintAll()
{
  if (m_hits_collection) m_hits_collection->PrintAllHits();
}

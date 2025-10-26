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

// 작은 헬퍼: 문자열 prefix 체크 (G4String은 std::string과 호환)
static inline bool StartsWith(const G4String& s, const char* prefix){
  return s.find(prefix) == 0; // 앞에서부터 일치
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

  // 중성입자 제외
  if (!def || def->GetPDGCharge() == 0.) return false;

  // --- 엔트리 조건 ---
  //   • 경계 진입(fGeomBoundary) 이거나
  //   • Boolean/UnionSolid의 첫 스텝이 fUndefined 로 들어오는 특수 케이스,
  //   • 또는 볼륨 내부에서 시작한 트랙의 첫 스텝(IsFirstStepInVolume)
  // 이 세 경우만 히트로 인정 (불필요한 다중 히트 방지)
  const auto stat = pre->GetStepStatus();
  const bool accept =
      aStep->IsFirstStepInVolume() ||
      (stat == fGeomBoundary) ||
      (stat == fUndefined);

  if (!accept) return false;

  // --- copy-no 추출(견고한 폴백 포함) ---
  int copyNo = -1;
  G4String pvName = "";

  // 1) pre-touchable에서 우선 시도
  if (const auto touch = pre->GetTouchable()){
    if (const auto pv = touch->GetVolume()){
      pvName = pv->GetName();
      copyNo = pv->GetCopyNo();
    }
  }
  // 2) pre가 비어있거나(드문 경우) 이름이 HTOF PV가 아니면 post로 폴백
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

  // 3) 그래도 유효하지 않으면(또는 다른 서브디텍터) 히트 버림
  if (copyNo < 0 || !StartsWith(pvName, "HtofPV"))
    return false;

  // --- 히트 생성 및 copy-no 저장 ---
  auto hit = new HTOFHit(SensitiveDetectorName, aStep);

  if (auto tp = hit->GetParticle()){
    // 분석에서 copy-no를 tile ID로 쓰기 위해 StatusCode에 저장
    tp->SetStatusCode(copyNo);
  }

  m_hits_collection->insert(hit);
  return true;
}

void HTOFSD::EndOfEvent(G4HCofThisEvent* /* HCTE */)
{
}

void HTOFSD::DrawAll()
{
}

void HTOFSD::PrintAll()
{
  if (m_hits_collection) m_hits_collection->PrintAllHits();
}

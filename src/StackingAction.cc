// -*- C++ -*-
// StackingAction.cc
// 목적: "진짜 원빔"만 태깅(OriginTag:kPrimaryBeam)해서 이후 SteppingAction의
//       VP4/VP5/타겟-스폰 로직을 원빔에만 적용.
// 2025-10-26 (user modification): 
//   - π+ 원빔( PDG=+211 ) 태깅 지원 추가
//   - 신규 제너레이터 7227(π+ beam-through) 지원 추가
//   - PDG=±211 모두 빔 정렬(cosθ) 판정 적용
//
// jaejin 2025-10-13: cosθ 가드 안전성 향상(영벡터 방지), conf로 임계값 제어 추가.

#include "StackingAction.hh"

#include <G4Track.hh>
#include <G4ParticleDefinition.hh>

#include "AnaManager.hh"
#include "ConfMan.hh"     // cosThr를 conf로 노출
#include "TrackTag.hh"    // SimpleTag(하위호환), OriginTag(신규)

namespace {
  auto& gAnaMan = AnaManager::GetInstance();
  const auto& gConf = ConfMan::GetInstance();
}

// ✅ [필수] 빈 정의들 (vtable 보존)
StackingAction::StackingAction() : G4UserStackingAction() {}
StackingAction::~StackingAction() = default;
void StackingAction::NewStage() { /* no-op */ }
void StackingAction::PrepareNewEvent() { /* no-op */ }

//----------------------------------------------------------------------------
//  원빔 정의(일반화):
//   - ParentID == 0 (primary)
//   - PDG == -211 (pi-) 또는 +211 (pi+)
//   - 빔파일에서 온 방향 +p̂_beam 과 충분히 정렬 (cosθ >= cosThr)
//   - 적용 generator: 7217 / 7218 / 7219 / 7227
//  태깅 내용:
//   - OriginTag(EOrigin::kPrimaryBeam)  [신규]
//   - (하위호환) 이전 코드가 SimpleTag를 기대하는 위치가 있으나,
//     UserInformation은 1개 포인터만 보유할 수 있어 OriginTag만 달고,
//     SteppingAction에서는 OriginTag 우선 + SimpleTag fallback을 지원함.
//----------------------------------------------------------------------------
G4ClassificationOfNewTrack
StackingAction::ClassifyNewTrack(const G4Track* track)
{
  // (0) 1차(primary)만 대상
  if (track->GetParentID() != 0) return fUrgent;

  // (1) 현재 generator 확인 (E45/E72 빔 모드들만)
  const int gen = gAnaMan.GetNextGenerator();
  const bool isBeamGen =
      (gen==7217 || gen==7218 || gen==7219   // π− 계열
    || gen==7227);                           // π+ beam-through (신규)
  if (!isBeamGen) return fUrgent;

  // (2) 입자: π− 또는 π+
  const auto* P = track->GetParticleDefinition();
  const int pdg = P ? P->GetPDGEncoding() : 0;
  if (pdg!= -211 && pdg!= +211) return fUrgent;

  // (3) 방향: +p̂_beam 과 정렬되어 있을 때만 (진짜 원빔으로 간주)
  //     - gAnaMan.GetNextMom(): BeamMan에서 읽어온 primary 빔 3-모멘텀(G4ThreeVector)
  const auto Pbeam = gAnaMan.GetNextMom();
  const double P2  = Pbeam.mag2();
  if (P2 <= 0.0) return fUrgent; // 방어(uninitialized beam vector)

  const auto n_beam = Pbeam.unit();
  const auto u_trk  = track->GetMomentumDirection().unit();

  // conf 키 예시: BeamCosThr: 0.97 (기본 0.95)
  double cosThr = 0.95;
  try { cosThr = gConf.Get<double>("BeamCosThr"); } catch(...) {}

  const double cosang = n_beam.dot(u_trk);
  if (cosang < cosThr) return fUrgent; // 빔축과 어긋난 primary(π±)는 태그하지 않음

  // (4) 태그 부여 (중복 방지)
  if (!track->GetUserInformation()) {
    // 신규 표준 태그: 출처를 명시
    track->SetUserInformation(new OriginTag(EOrigin::kPrimaryBeam, -1,
      (pdg==+211? "beam_pi+" : "beam_pi-")));
    // (하위호환) SimpleTag는 생략. SteppingAction에서 fallback 처리 지원.
  }

  return fUrgent;
}

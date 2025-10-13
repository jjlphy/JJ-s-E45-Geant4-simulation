// -*- C++ -*-
// StackingAction.cc
// 목적: "진짜 원빔"만 태깅(OriginTag:kPrimaryBeam)해서 이후 SteppingAction의
//       VP4/VP5/타겟-스폰 로직을 원빔에만 적용.
// jaejin 2025-10-13: cosθ 가드 안전성 향상(영벡터 방지), conf로 임계값 제어 추가.

#include "StackingAction.hh"

#include <G4Track.hh>
#include <G4ParticleDefinition.hh>

#include "AnaManager.hh"
#include "ConfMan.hh"     // jaejin 2025-10-13: cosThr를 conf로 뺄 수 있게
#include "TrackTag.hh"    // SimpleTag(하위호환), OriginTag(신규)

namespace {
  auto& gAnaMan = AnaManager::GetInstance();
  const auto& gConf = ConfMan::GetInstance();
}

// ✅ [필수] 여기 4개 “정의” 추가
StackingAction::StackingAction() : G4UserStackingAction() {}                 // jaejin 2025-10-13: vtable 생성을 위해 빈 정의 추가
StackingAction::~StackingAction() = default;                                 // jaejin 2025-10-13: 선언만 있던 소멸자 정의
void StackingAction::NewStage() { /* no-op */ }                              // jaejin 2025-10-13: 선언만 있던 가상함수 정의
void StackingAction::PrepareNewEvent() { /* no-op */ }                       // jaejin 2025-10-13: 선언만 있던 가상함수 정의

//----------------------------------------------------------------------------
//  원빔 정의:
//   - ParentID == 0 (primary)
//   - PDG == -211 (pi-)
//   - 빔파일에서 온 방향 +p̂ 과 충분히 정렬 (cosθ >= cosThr)
//   - 적용 generator: 7217 / 7218 / 7219
//  태깅 내용:
//   - OriginTag(EOrigin::kPrimaryBeam)  [신규]
//   - (하위호환) 이전 코드가 기대하는 SimpleTag("UPSTREAM_HELPER")도 함께 부여
//----------------------------------------------------------------------------
G4ClassificationOfNewTrack
StackingAction::ClassifyNewTrack(const G4Track* track)
{
  // (0) 1차(primary)만 대상
  if (track->GetParentID() != 0) return fUrgent;

  // (1) 현재 generator 확인 (E45/E72 빔 모드들만)
  const int gen = gAnaMan.GetNextGenerator();
  if (gen != 7217 && gen != 7218 && gen != 7219) return fUrgent;

  // (2) 입자: pi-
  const auto* P = track->GetParticleDefinition();
  if (!P || P->GetPDGEncoding() != -211) return fUrgent;

  // (3) 방향: +p̂_beam 과 정렬되어 있을 때만 (진짜 원빔)
  //     - gAnaMan.GetNextMom(): BeamMan에서 읽어온 primary 빔 3-모멘텀(G4ThreeVector)
  const auto Pbeam = gAnaMan.GetNextMom();
  const double P2  = Pbeam.mag2();
  if (P2 <= 0.0) return fUrgent; // jaejin 2025-10-13: 방어 (unit() NaN 방지)

  const auto n_beam = Pbeam.unit();
  const auto u_trk  = track->GetMomentumDirection().unit();
  // conf 키 예시: BeamCosThr: 0.97
  double cosThr = 0.95;
  try { cosThr = gConf.Get<double>("BeamCosThr"); } catch(...) {}

  const double cosang = n_beam.dot(u_trk);
  if (cosang < cosThr) return fUrgent; // 빔축과 어긋난 primary(pi-)는 태그하지 않음

  // (4) 태그 부여 (중복 방지)
  if (!track->GetUserInformation()) {
    // 신규 표준 태그: 출처를 명시
    track->SetUserInformation(new OriginTag(EOrigin::kPrimaryBeam, -1, "beam_pi-"));
    // 하위호환: 기존 코드가 SimpleTag를 기대하는 위치가 있으므로 함께 부여
    // (주의: 둘 다 UserInformation에 넣을 수 없어 SimpleTag는 생략하고,
    //        SteppingAction에서 OriginTag 우선 사용 + SimpleTag fallback 지원)
    // -> 만약 꼭 SimpleTag가 필요하면 OriginTag 안에 why 문자열을 추가하거나,
    //    StackingAction에서 SimpleTag만 달고, SteppingAction 첫 스텝에서
    //    OriginTag를 보강 부여하는 방법도 가능.
  } else {
    // 이미 다른 태그가 있다면, 가능하면 OriginTag로 업그레이드(선택)
    // auto* old = track->GetUserInformation();
    // (여기서는 무해하게 통과)
  }

  return fUrgent;
}

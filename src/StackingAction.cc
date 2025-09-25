// -*- C++ -*-
#include "StackingAction.hh"

#include <G4Track.hh>
#include <G4ParticleDefinition.hh>

#include "AnaManager.hh"
#include "TrackTag.hh"  // SimpleTag

namespace {
  auto& gAnaMan = AnaManager::GetInstance();
}

//----------------------------------------------------------------------------
//  목적: "진짜 원빔"만 태깅해서 SteppingAction의 VP4/VP5/타겟-스폰 로직이
//        오직 원빔에만 적용되도록 한다.
//
//  원빔 정의:
//   - ParentID == 0 (primary)
//   - PDG == -211 (pi-)
//   - 빔파일에서 온 방향 +p̂ 과 충분히 정렬 (cosθ >= cosThr)
//   - 적용 generator: 7217 / 7218 / 7219
//----------------------------------------------------------------------------
G4ClassificationOfNewTrack
StackingAction::ClassifyNewTrack(const G4Track* track)
{
  // (0) 1차(primary)만 대상
  if (track->GetParentID() != 0) return fUrgent;

  // (1) 현재 generator 확인 (E45 빔 모드들만)
  const int gen = gAnaMan.GetNextGenerator();
  if (gen != 7217 && gen != 7218 && gen != 7219) return fUrgent;

  // (2) 입자: pi-
  const auto* P = track->GetParticleDefinition();
  if (!P || P->GetPDGEncoding() != -211) return fUrgent;

  // (3) 방향: +p̂_beam 과 정렬되어 있을 때만 (진짜 원빔)
  //     - gAnaMan.GetNextMom(): BeamMan에서 읽어온 primary 빔 3-모멘텀(G4ThreeVector)
  const auto n_beam = gAnaMan.GetNextMom().unit();
  const auto u_trk  = track->GetMomentumDirection().unit();
  const double cosang = n_beam.dot(u_trk);
  static constexpr double cosThr = 0.95; // 필요시 0.97~0.99로 조절 가능

  if (cosang < cosThr) return fUrgent;   // 빔축과 어긋난 primary(pi-)는 태그하지 않음

  // (4) 태그 부여 (중복 방지)
  if (!track->GetUserInformation()) {
    track->SetUserInformation(new SimpleTag("UPSTREAM_HELPER"));
  }
  return fUrgent;
}


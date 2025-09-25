#include "StackingAction.hh"

#include <G4Track.hh>
#include <G4ParticleDefinition.hh>

#include "AnaManager.hh"
#include "TrackTag.hh"  // SimpleTag

namespace {
  auto& gAnaMan = AnaManager::GetInstance();
}

G4ClassificationOfNewTrack
StackingAction::ClassifyNewTrack(const G4Track* track)
{
  // (1) 1차(primary)만 대상
  if (track->GetParentID() != 0) return fUrgent;

  // (2) PDG = -211 (pi-) : 진짜 빔만 태그
  const auto* P = track->GetParticleDefinition();
  if (!P || P->GetPDGEncoding() != -211) return fUrgent;

  // (3) 빔이면 무조건 helper 태그 (generator 번호/방향 조건 없음)
  if (!track->GetUserInformation()) {
    track->SetUserInformation(new SimpleTag("UPSTREAM_HELPER"));
  }
  return fUrgent;
}

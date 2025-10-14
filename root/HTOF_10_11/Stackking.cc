// -*- C++ -*-
// StackingAction.cc
// Goal: Tag *only true primary beam* (OriginTag:kPrimaryBeam) so that
//       VP4/VP5/target-spawn logic later applies *only to the beam track*.
// jaejin 2025-10-13: safer cosθ guard (avoid unit() on zero-vector),
//                    make the threshold configurable via conf.

#include "StackingAction.hh"

#include <G4Track.hh>
#include <G4ParticleDefinition.hh>

#include "AnaManager.hh"
#include "ConfMan.hh"     // jaejin 2025-10-13: to read BeamCosThr from conf
#include "TrackTag.hh"    // SimpleTag (legacy), OriginTag (new)

namespace {
  auto& gAnaMan = AnaManager::GetInstance();
  const auto& gConf = ConfMan::GetInstance();
}

// Provide out-of-line definitions to generate vtable (required by linker)
StackingAction::StackingAction() : G4UserStackingAction() {}                 // jaejin 2025-10-13
StackingAction::~StackingAction() = default;                                 // jaejin 2025-10-13
void StackingAction::NewStage() { /* no-op */ }                              // jaejin 2025-10-13
void StackingAction::PrepareNewEvent() { /* no-op */ }                       // jaejin 2025-10-13

//----------------------------------------------------------------------------
//  True primary beam definition (selection order):
//   - ParentID == 0 (primary)
//   - PDG == -211 (pi-)
//   - Direction aligned with beam file expectation (cosθ >= cosThr)
//   - Apply only to generators: 7217 / 7218 / 7219
//  Tagging:
//   - OriginTag(EOrigin::kPrimaryBeam)  [new]
//   - (Legacy) Older code that checks SimpleTag("UPSTREAM_HELPER") will still work
//----------------------------------------------------------------------------
G4ClassificationOfNewTrack
StackingAction::ClassifyNewTrack(const G4Track* track)
{
  // (0) Only consider primary tracks
  if (track->GetParentID() != 0) return fUrgent;

  // (1) Limit to specific generators (E45/E72 beam modes only)
  const int gen = gAnaMan.GetNextGenerator();
  if (gen != 7217 && gen != 7218 && gen != 7219) return fUrgent;

  // (2) Particle must be pi-
  const auto* P = track->GetParticleDefinition();
  if (!P || P->GetPDGEncoding() != -211) return fUrgent;

  // (3) Direction must be aligned with +p̂_beam (true beam only)
  //     - gAnaMan.GetNextMom(): primary beam 3-momentum (G4ThreeVector) read by BeamMan
  const auto Pbeam = gAnaMan.GetNextMom();
  const double P2  = Pbeam.mag2();
  if (P2 <= 0.0) return fUrgent; // jaejin 2025-10-13: guard (avoid unit() NaN on zero-length vector)

  const auto n_beam = Pbeam.unit();
  const auto u_trk  = track->GetMomentumDirection().unit();

  // Conf key example: BeamCosThr: 0.97
  double cosThr = 0.95;
  try { cosThr = gConf.Get<double>("BeamCosThr"); } catch(...) {}

  const double cosang = n_beam.dot(u_trk);
  if (cosang < cosThr) return fUrgent; // not well-aligned with beam axis → do not tag as primary beam

  // (4) Tag (no duplication)
  if (!track->GetUserInformation()) {
    // Standard origin tag
    track->SetUserInformation(new OriginTag(EOrigin::kPrimaryBeam, -1, "beam_pi-"));
    // Legacy note:
    //   We do not attach SimpleTag here because G4 keeps only one UserInformation pointer.
    //   SteppingAction supports legacy SimpleTag as a fallback when present.
  } else {
    // If some other tag exists, we leave it as-is (optional: could upgrade to OriginTag)
  }

  return fUrgent;
}

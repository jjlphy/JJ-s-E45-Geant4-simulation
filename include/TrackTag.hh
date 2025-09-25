//jaejin 25.9.25

// -*- C++ -*-
// TrackTag.hh
#ifndef TRACK_TAG_HH
#define TRACK_TAG_HH

#include <G4VUserTrackInformation.hh>

// StackingAction에서 달고, SteppingAction에서 읽을 공용 태그
struct SimpleTag : public G4VUserTrackInformation {
  explicit SimpleTag(const char* why) : why_(why) {}
  const char* why_;  // 예: "UPSTREAM_HELPER"
};

#endif // TRACK_TAG_HH

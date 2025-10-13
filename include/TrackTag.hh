// -*- C++ -*-
// TrackTag.hh
// jaejin 2025-10-13: 기존 SimpleTag(UPSTREAM_HELPER)를 보존하면서,
//                    원빔과 강제 2π 자식을 확실히 구분하기 위해 OriginTag를 도입.
//                    오프라인에서 origin, channel 기반으로 간단히 컷할 수 있도록 한다.

#ifndef TRACK_TAG_HH
#define TRACK_TAG_HH

#include <G4VUserTrackInformation.hh>

// ────────────────────────────────────────────────────────────────────────────
// [하위호환] 간단 태그 (기존 코드 유지)
// ────────────────────────────────────────────────────────────────────────────
struct SimpleTag : public G4VUserTrackInformation {
  explicit SimpleTag(const char* why) : why_(why) {}
  const char* why_;  // 예: "UPSTREAM_HELPER"
  // jaejin 2025-10-13: 기존 SteppingAction에서 쓰던 태그. 지우지 말고 유지.
};

// ────────────────────────────────────────────────────────────────────────────
// [신규] 출처 태그 (원빔/강제2π 자식 구분 + 채널 저장)
// ────────────────────────────────────────────────────────────────────────────
enum class EOrigin : int {
  kUnknown        = 0, // 기본값/미지정
  kPrimaryBeam    = 1, // 원빔(ParentID==0, pi-)
  kForced2PiChild = 2  // 타겟 첫진입에서 인위적으로 생성한 2π 자식들
};

struct OriginTag : public G4VUserTrackInformation {
  EOrigin   origin;   // 출처
  int       channel;  // 0:(π+ π− n), 1:(π− π0 p), -1: 없음
  const char* name;   // 선택: "pi-","pi+","pi0","neutron","proton"
  OriginTag(EOrigin o=EOrigin::kUnknown, int ch=-1, const char* nm=nullptr)
  : origin(o), channel(ch), name(nm) {}
  // jaejin 2025-10-13: 오프라인에서 origin==kForced2PiChild && pdg==-211 조건으로
  //                    "2π에서 나온 π−"만 간단히 선택 가능.
};

// ────────────────────────────────────────────────────────────────────────────
// [헬퍼] 안전하게 OriginTag 읽기/체크
// ────────────────────────────────────────────────────────────────────────────
inline OriginTag* GetOriginTag(const G4Track* trk){
  return trk ? dynamic_cast<OriginTag*>(trk->GetUserInformation()) : nullptr;
}
inline bool HasOrigin(const G4Track* trk, EOrigin o){
  if(!trk) return false;
  auto* t = GetOriginTag(trk);
  return (t && t->origin==o);
}

#endif // TRACK_TAG_HH

// -*- C++ -*-

#ifndef TPC_STEPPING_ACTION_HH
#define TPC_STEPPING_ACTION_HH

#include <G4UserSteppingAction.hh>
#include <G4ThreeVector.hh>      // 추가 *jaejin_25_10_25*
#include <unordered_map>         // 추가
#include <G4Types.hh>            // G4int 등

//_____________________________________________________________________________
class SteppingAction : public G4UserSteppingAction
{
public:
  SteppingAction();
  virtual ~SteppingAction();
  virtual void UserSteppingAction(const G4Step* theStep);

private:
  /// === σ 불요 트랙-균일: 트랙별 상태 (이 파일 내부에서만 사용)
struct UniAlongTrack {
  bool          active  = false;   // LH2 내부 경로 수집 중인가
  bool          hasCand = false;   // 후보 버텍스 존재 여부
  double        S       = 0.0;     // 누적 경로 길이 (mm)
  G4ThreeVector vtx;               // 후보 버텍스 (위치)
  G4ThreeVector pcand;             // ★ 후보 지점의 모멘텀(3-벡터)
  G4double      tcand = 0.0;       // (선택) 후보 지점의 global time
};

  std::unordered_map<G4int, UniAlongTrack> m_uat; // key: TrackID
};

#endif
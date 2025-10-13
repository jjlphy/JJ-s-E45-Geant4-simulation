// -*- C++ -*-
#ifndef STACKING_ACTION_HH
#define STACKING_ACTION_HH

#include <G4UserStackingAction.hh>
class G4Track;

class StackingAction : public G4UserStackingAction {
public:
  StackingAction();
  ~StackingAction() override;  // 반드시 override + 정의 제공

  G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* track) override;

  // 굳이 선언 안 해도 되지만, 선언했다면 .cc에 빈 정의 넣기
  void NewStage() override;
  void PrepareNewEvent() override;
};

#endif

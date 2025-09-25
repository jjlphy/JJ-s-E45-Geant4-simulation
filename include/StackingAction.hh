// StackingAction.hh
#ifndef STACKING_ACTION_HH
#define STACKING_ACTION_HH

#include <G4UserStackingAction.hh>

class StackingAction : public G4UserStackingAction {
public:
  StackingAction() = default;
  ~StackingAction() override = default;

  G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* track) override;
};

#endif

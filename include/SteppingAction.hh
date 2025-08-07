// -*- C++ -*-

#ifndef TPC_STEPPING_ACTION_HH
#define TPC_STEPPING_ACTION_HH

#include <G4UserSteppingAction.hh>

//_____________________________________________________________________________
class SteppingAction : public G4UserSteppingAction
{
public:
  SteppingAction();
  virtual ~SteppingAction();
  virtual void UserSteppingAction(const G4Step* theStep);
};

#endif

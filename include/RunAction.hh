// -*- C++ -*-

#ifndef RUN_ACTION_HH
#define RUN_ACTION_HH

#include <G4UserRunAction.hh>
#include <G4String.hh>

class G4Run;
class TPCAnaManager;

//_____________________________________________________________________________
class RunAction : public G4UserRunAction
{
public:
  static G4String ClassName();
  RunAction();
  virtual ~RunAction();
  virtual void BeginOfRunAction(const G4Run* aRun);
  virtual void EndOfRunAction(const G4Run* aRun);
};

//_____________________________________________________________________________
inline G4String
RunAction::ClassName()
{
  static G4String s_name("RunAction");
  return s_name;
}

#endif

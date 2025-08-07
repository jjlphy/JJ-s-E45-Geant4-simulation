// -*- C++ -*-

#include "ActionInitialization.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"

//_____________________________________________________________________________
ActionInitialization::ActionInitialization()
{
}

//_____________________________________________________________________________
ActionInitialization::~ActionInitialization()
{}

//_____________________________________________________________________________
void
ActionInitialization::BuildForMaster() const
{
  SetUserAction(new RunAction);
}

//_____________________________________________________________________________
void
ActionInitialization::Build() const
{
  SetUserAction(new PrimaryGeneratorAction);
  SetUserAction(new RunAction);
  SetUserAction(new EventAction);
  SetUserAction(new SteppingAction);
}

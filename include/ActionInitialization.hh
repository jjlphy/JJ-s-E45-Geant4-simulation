// -*- C++ -*-

#ifndef ACTION_INITIALIZATION_HH
#define ACTION_INITIALIZATION_HH

#include "G4VUserActionInitialization.hh"

//_____________________________________________________________________________
class ActionInitialization : public G4VUserActionInitialization
{
public:
  ActionInitialization();
  ~ActionInitialization() override;

  void BuildForMaster() const override;
  void Build() const override;
};

#endif

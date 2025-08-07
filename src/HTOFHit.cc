// -*- C++ -*-

#include "HTOFHit.hh"

#include <CLHEP/Units/PhysicalConstants.h>
#include <G4Step.hh>

G4Allocator<HTOFHit> HTOFHitAllocator;

//_____________________________________________________________________________
HTOFHit::HTOFHit(const G4String& name, G4Step* step)
  : G4VHit(),
    VHitInfo(name, step)
{
}

//_____________________________________________________________________________
HTOFHit::~HTOFHit()
{
}

//_____________________________________________________________________________
void
HTOFHit::Draw()
{
}

//_____________________________________________________________________________
void
HTOFHit::Print()
{
  VHitInfo::Print();
}

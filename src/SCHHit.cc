// -*- C++ -*-

#include "SCHHit.hh"

#include <CLHEP/Units/PhysicalConstants.h>
#include <G4Step.hh>

G4Allocator<SCHHit> SCHHitAllocator;

//_____________________________________________________________________________
SCHHit::SCHHit(const G4String& name, G4Step* step)
  : G4VHit(),
    VHitInfo(name, step)
{
}

//_____________________________________________________________________________
SCHHit::~SCHHit()
{
}

//_____________________________________________________________________________
void
SCHHit::Draw()
{
}

//_____________________________________________________________________________
void
SCHHit::Print()
{
  VHitInfo::Print();
}

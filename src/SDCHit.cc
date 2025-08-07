// -*- C++ -*-

#include "SDCHit.hh"

#include <CLHEP/Units/PhysicalConstants.h>
#include <G4Step.hh>

G4Allocator<SDCHit> SDCHitAllocator;

//_____________________________________________________________________________
SDCHit::SDCHit(const G4String& name, G4Step* step)
  : G4VHit(),
    VHitInfo(name, step)
{
}

//_____________________________________________________________________________
SDCHit::~SDCHit()
{
}

//_____________________________________________________________________________
void
SDCHit::Draw()
{
}

//_____________________________________________________________________________
void
SDCHit::Print()
{
  VHitInfo::Print();
}

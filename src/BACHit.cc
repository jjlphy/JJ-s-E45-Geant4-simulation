// -*- C++ -*-

#include "BACHit.hh"

#include <CLHEP/Units/PhysicalConstants.h>
#include <G4Step.hh>

G4Allocator<BACHit> BACHitAllocator;

//_____________________________________________________________________________
BACHit::BACHit(const G4String& name, G4Step* step)
  : G4VHit(),
    VHitInfo(name, step)
{
}

//_____________________________________________________________________________
BACHit::~BACHit()
{
}

//_____________________________________________________________________________
void
BACHit::Draw()
{
}

//_____________________________________________________________________________
void
BACHit::Print()
{
  VHitInfo::Print();
}

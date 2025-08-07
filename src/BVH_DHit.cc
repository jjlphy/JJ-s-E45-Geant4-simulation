// -*- C++ -*-

#include "BVH_DHit.hh"

#include <CLHEP/Units/PhysicalConstants.h>
#include <G4Step.hh>

G4Allocator<BVH_DHit> BVH_DHitAllocator;

//_____________________________________________________________________________
BVH_DHit::BVH_DHit(const G4String& name, G4Step* step)
  : G4VHit(),
    VHitInfo(name, step)
{
}

//_____________________________________________________________________________
BVH_DHit::~BVH_DHit()
{
}

//_____________________________________________________________________________
void
BVH_DHit::Draw()
{
}

//_____________________________________________________________________________
void
BVH_DHit::Print()
{
  VHitInfo::Print();
}

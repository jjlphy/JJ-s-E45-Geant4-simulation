// -*- C++ -*-

#include "BVH_UHit.hh"

#include <CLHEP/Units/PhysicalConstants.h>
#include <G4Step.hh>

G4Allocator<BVH_UHit> BVH_UHitAllocator;

//_____________________________________________________________________________
BVH_UHit::BVH_UHit(const G4String& name, G4Step* step)
  : G4VHit(),
    VHitInfo(name, step)
{
}

//_____________________________________________________________________________
BVH_UHit::~BVH_UHit()
{
}

//_____________________________________________________________________________
void
BVH_UHit::Draw()
{
}

//_____________________________________________________________________________
void
BVH_UHit::Print()
{
  VHitInfo::Print();
}

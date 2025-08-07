// -*- C++ -*-

#include "BH2Hit.hh"

#include <CLHEP/Units/PhysicalConstants.h>
#include <G4Step.hh>

G4Allocator<BH2Hit> BH2HitAllocator;

//_____________________________________________________________________________
BH2Hit::BH2Hit(const G4String& name, G4Step* step)
  : G4VHit(),
    VHitInfo(name, step)
{
}

//_____________________________________________________________________________
BH2Hit::~BH2Hit()
{
}

//_____________________________________________________________________________
void
BH2Hit::Draw()
{
}

//_____________________________________________________________________________
void
BH2Hit::Print()
{
  VHitInfo::Print();
}

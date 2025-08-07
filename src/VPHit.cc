// -*- C++ -*-

#include "VPHit.hh"

#include <CLHEP/Units/PhysicalConstants.h>
#include <G4Step.hh>

G4Allocator<VPHit> VPHitAllocator;

//_____________________________________________________________________________
VPHit::VPHit(const G4String& name, G4Step* step)
  : G4VHit(),
    VHitInfo(name, step)
{
}

//_____________________________________________________________________________
VPHit::~VPHit()
{
}

//_____________________________________________________________________________
void
VPHit::Draw()
{
}

//_____________________________________________________________________________
void
VPHit::Print()
{
  VHitInfo::Print();
}

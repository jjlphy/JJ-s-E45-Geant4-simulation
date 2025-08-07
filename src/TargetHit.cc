// -*- C++ -*-

#include "TargetHit.hh"

#include <CLHEP/Units/PhysicalConstants.h>
#include <G4Step.hh>

G4Allocator<TargetHit> TargetHitAllocator;

//_____________________________________________________________________________
TargetHit::TargetHit(const G4String& name, G4Step* step,
                     G4bool use_center_point)
  : G4VHit(),
    VHitInfo(name, step, use_center_point)
{
}

//_____________________________________________________________________________
TargetHit::~TargetHit()
{
}

//_____________________________________________________________________________
void
TargetHit::Draw()
{
}

//_____________________________________________________________________________
void
TargetHit::Print()
{
  VHitInfo::Print();
}

// -*- C++ -*-

#include "LACHit.hh"

#include <CLHEP/Units/PhysicalConstants.h>
#include <G4Step.hh>

G4Allocator<LACHit> LACHitAllocator;

//_____________________________________________________________________________
LACHit::LACHit(const G4String& name, G4Step* step)
  : G4VHit(),
    VHitInfo(name, step)
{
}

//_____________________________________________________________________________
LACHit::~LACHit()
{
}

//_____________________________________________________________________________
void
LACHit::Draw()
{
}

//_____________________________________________________________________________
void
LACHit::Print()
{
  VHitInfo::Print();
}

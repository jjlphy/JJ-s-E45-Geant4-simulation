// -*- C++ -*-

#include "TPCHit.hh"

#include <CLHEP/Units/PhysicalConstants.h>

G4Allocator<TPCHit> TPCHitAllocator;

//_____________________________________________________________________________
TPCHit::TPCHit(const G4String& name, G4Step* step)
  : G4VHit(),
    VHitInfo(name, step)
{
}

//_____________________________________________________________________________
TPCHit::~TPCHit()
{
}

//_____________________________________________________________________________
void
TPCHit::Draw()
{
}

//_____________________________________________________________________________
void
TPCHit::Print()
{
  VHitInfo::Print();
}

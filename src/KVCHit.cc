// -*- C++ -*-

#include "KVCHit.hh"

#include <CLHEP/Units/PhysicalConstants.h>
#include <G4Step.hh>

G4Allocator<KVCHit> KVCHitAllocator;

//_____________________________________________________________________________
KVCHit::KVCHit(const G4String& name, G4Step* step)
  : G4VHit(),
    VHitInfo(name, step)
{
}

//_____________________________________________________________________________
KVCHit::~KVCHit()
{
}

//_____________________________________________________________________________
void
KVCHit::Draw()
{
}

//_____________________________________________________________________________
void
KVCHit::Print()
{
  VHitInfo::Print();
}

// -*- C++ -*-

#include "WCHit.hh"

#include <CLHEP/Units/PhysicalConstants.h>
#include <G4Step.hh>

G4Allocator<WCHit> WCHitAllocator;

//_____________________________________________________________________________
WCHit::WCHit(const G4String& name, G4Step* step)
  : G4VHit(),
    VHitInfo(name, step)
{
}

//_____________________________________________________________________________
WCHit::~WCHit()
{
}

//_____________________________________________________________________________
void
WCHit::Draw()
{
}

//_____________________________________________________________________________
void
WCHit::Print()
{
  VHitInfo::Print();
}

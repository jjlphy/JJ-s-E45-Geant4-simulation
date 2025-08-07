// -*- C++ -*-

#include "FTOFHit.hh"

#include <CLHEP/Units/PhysicalConstants.h>
#include <G4Step.hh>

G4Allocator<FTOFHit> FTOFHitAllocator;

//_____________________________________________________________________________
FTOFHit::FTOFHit(const G4String& name, G4Step* step)
  : G4VHit(),
    VHitInfo(name, step)
{
}

//_____________________________________________________________________________
FTOFHit::~FTOFHit()
{
}

//_____________________________________________________________________________
void
FTOFHit::Draw()
{
}

//_____________________________________________________________________________
void
FTOFHit::Print()
{
  VHitInfo::Print();
}

#include "T0Hit.hh"

#include <G4Step.hh>
#include <G4VPhysicalVolume.hh>
#include <G4ThreeVector.hh>
#include <G4UnitsTable.hh>
#include <TParticle.h>

G4Allocator<T0Hit> T0HitAllocator;

T0Hit::T0Hit(const G4String& name, G4Step* step)
  : VHitInfo(name, step) {}

T0Hit::~T0Hit() {}

void T0Hit::Draw() {}

void T0Hit::Print()
{
  G4cout << "T0Hit: detector=" << GetDetectorName()
         << " particle=" << GetParticle()->GetName()
         << " pos=(" << GetParticle()->Vx() << ", "
                     << GetParticle()->Vy() << ", "
                     << GetParticle()->Vz() << ")" << G4endl;
}

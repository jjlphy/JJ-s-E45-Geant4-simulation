// -*- C++ -*-

#ifndef VP_HIT_HH
#define VP_HIT_HH

#include <G4Allocator.hh>
#include <G4String.hh>
#include <G4THitsCollection.hh>
#include <G4ThreeVector.hh>
#include <G4VHit.hh>

#include "VHitInfo.hh"

class G4Step;

//_____________________________________________________________________________
class VPHit : public G4VHit, public VHitInfo
{
public:
  VPHit(const G4String& name, G4Step* step);
  virtual ~VPHit();

  VPHit(const VPHit& right);
  const VPHit& operator=(const VPHit& right);

  void* operator new(size_t);
  void operator delete(void* aHit);

public:
  virtual void Draw();
  virtual void Print();
};

//_____________________________________________________________________________
inline
VPHit::VPHit(const VPHit& right)
  : G4VHit(right),
    VHitInfo(right)
{
}

//_____________________________________________________________________________
inline const VPHit&
VPHit::operator =(const VPHit& right)
{
  VHitInfo::operator =(right);
  return *this;
}

//_____________________________________________________________________________
// externally instanciated.
extern G4Allocator<VPHit> VPHitAllocator;

//_____________________________________________________________________________
inline void*
VPHit::operator new(size_t)
{
  return VPHitAllocator.MallocSingle();
}

//_____________________________________________________________________________
inline void
VPHit::operator delete(void* aHit)
{
  VPHitAllocator.FreeSingle(static_cast<VPHit*>(aHit));
}

#endif

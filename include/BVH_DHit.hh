// -*- C++ -*-

#ifndef BVH_D_HIT_HH
#define BVH_D_HIT_HH

#include <G4Allocator.hh>
#include <G4String.hh>
#include <G4THitsCollection.hh>
#include <G4ThreeVector.hh>
#include <G4VHit.hh>

#include "VHitInfo.hh"

class G4Step;

//_____________________________________________________________________________
class BVH_DHit : public G4VHit, public VHitInfo
{
public:
  BVH_DHit(const G4String& name, G4Step* step);
  virtual ~BVH_DHit();

  BVH_DHit(const BVH_DHit& right);
  const BVH_DHit& operator=(const BVH_DHit& right);

  void* operator new(size_t);
  void operator delete(void* aHit);

public:
  virtual void Draw();
  virtual void Print();
};

//_____________________________________________________________________________
inline
BVH_DHit::BVH_DHit(const BVH_DHit& right)
  : G4VHit(right),
    VHitInfo(right)
{
}

inline const BVH_DHit&
BVH_DHit::operator=(const BVH_DHit& right)
{
  VHitInfo::operator=(right);
  return *this;
}

// externally instantiated
extern G4Allocator<BVH_DHit> BVH_DHitAllocator;

inline void* BVH_DHit::operator new(size_t)
{
  return (void*) BVH_DHitAllocator.MallocSingle();
}

inline void BVH_DHit::operator delete(void* aHit)
{
  BVH_DHitAllocator.FreeSingle((BVH_DHit*) aHit);
}

#endif

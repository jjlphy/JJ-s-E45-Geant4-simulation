#ifndef BVH_U_HIT_HH
#define BVH_U_HIT_HH

#include <G4Allocator.hh>
#include <G4String.hh>
#include <G4THitsCollection.hh>
#include <G4ThreeVector.hh>
#include <G4VHit.hh>

#include "VHitInfo.hh"

class G4Step;

//_____________________________________________________________________________
class BVH_UHit : public G4VHit, public VHitInfo
{
public:
  BVH_UHit(const G4String& name, G4Step* step);
  virtual ~BVH_UHit();

  BVH_UHit(const BVH_UHit& right);
  const BVH_UHit& operator=(const BVH_UHit& right);

  void* operator new(size_t);
  void operator delete(void* aHit);

public:
  virtual void Draw();
  virtual void Print();
};

//_____________________________________________________________________________
inline
BVH_UHit::BVH_UHit(const BVH_UHit& right)
  : G4VHit(right),
    VHitInfo(right)
{
}

inline const BVH_UHit&
BVH_UHit::operator =(const BVH_UHit& right)
{
  VHitInfo::operator =(right);
  return *this;
}

// externally instantiated
extern G4Allocator<BVH_UHit> BVH_UHitAllocator;

inline void* BVH_UHit::operator new(size_t)
{
  return (void*) BVH_UHitAllocator.MallocSingle();
}

inline void BVH_UHit::operator delete(void* aHit)
{
  BVH_UHitAllocator.FreeSingle((BVH_UHit*) aHit);
}

#endif

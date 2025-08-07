// -*- C++ -*-

#ifndef BH2_HIT_HH
#define BH2_HIT_HH

#include <G4Allocator.hh>
#include <G4String.hh>
#include <G4THitsCollection.hh>
#include <G4ThreeVector.hh>
#include <G4VHit.hh>

#include "VHitInfo.hh"

class G4Step;

//_____________________________________________________________________________
class BH2Hit : public G4VHit, public VHitInfo
{
public:
  BH2Hit(const G4String& name, G4Step* step);
  virtual ~BH2Hit();

  BH2Hit(const BH2Hit& right);
  const BH2Hit& operator=(const BH2Hit& right);

  void* operator new(size_t);
  void operator delete(void* aHit);

public:
  virtual void Draw();
  virtual void Print();
};

//_____________________________________________________________________________
inline
BH2Hit::BH2Hit(const BH2Hit& right)
  : G4VHit(right),
    VHitInfo(right)
{
}

//_____________________________________________________________________________
inline const BH2Hit&
BH2Hit::operator =(const BH2Hit& right)
{
  VHitInfo::operator =(right);
  return *this;
}

//_____________________________________________________________________________
// externally instanciated.
extern G4Allocator<BH2Hit> BH2HitAllocator;

//_____________________________________________________________________________
inline void*
BH2Hit::operator new(size_t)
{
  return BH2HitAllocator.MallocSingle();
}

//_____________________________________________________________________________
inline void
BH2Hit::operator delete(void* aHit)
{
  BH2HitAllocator.FreeSingle(static_cast<BH2Hit*>(aHit));
}

#endif

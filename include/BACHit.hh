// -*- C++ -*-

#ifndef BAC_HIT_HH
#define BAC_HIT_HH

#include <G4Allocator.hh>
#include <G4String.hh>
#include <G4THitsCollection.hh>
#include <G4ThreeVector.hh>
#include <G4VHit.hh>

#include "VHitInfo.hh"

class G4Step;

//_____________________________________________________________________________
class BACHit : public G4VHit, public VHitInfo
{
public:
  BACHit(const G4String& name, G4Step* step);
  virtual ~BACHit();

  BACHit(const BACHit& right);
  const BACHit& operator=(const BACHit& right);

  void* operator new(size_t);
  void operator delete(void* aHit);

public:
  virtual void Draw();
  virtual void Print();
};

//_____________________________________________________________________________
inline
BACHit::BACHit(const BACHit& right)
  : G4VHit(right),
    VHitInfo(right)
{
}

//_____________________________________________________________________________
inline const BACHit&
BACHit::operator =(const BACHit& right)
{
  VHitInfo::operator =(right);
  return *this;
}

//_____________________________________________________________________________
// externally instanciated.
extern G4Allocator<BACHit> BACHitAllocator;

//_____________________________________________________________________________
inline void*
BACHit::operator new(size_t)
{
  void* aHit = (void*)BACHitAllocator.MallocSingle();
  return aHit;
}

//_____________________________________________________________________________
inline void
BACHit::operator delete(void* aHit)
{
  BACHitAllocator.FreeSingle((BACHit*) aHit);
}

#endif

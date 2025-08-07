// -*- C++ -*-

#ifndef SCH_HIT_HH
#define SCH_HIT_HH

#include <G4Allocator.hh>
#include <G4String.hh>
#include <G4THitsCollection.hh>
#include <G4ThreeVector.hh>
#include <G4VHit.hh>

#include "VHitInfo.hh"

class G4Step;

//_____________________________________________________________________________
class SCHHit : public G4VHit, public VHitInfo
{
public:
  SCHHit(const G4String& name, G4Step* step);
  virtual ~SCHHit();

  SCHHit(const SCHHit& right);
  const SCHHit& operator=(const SCHHit& right);

  void* operator new(size_t);
  void operator delete(void* aHit);

public:
  virtual void Draw();
  virtual void Print();
};

//_____________________________________________________________________________
inline
SCHHit::SCHHit(const SCHHit& right)
  : G4VHit(right),
    VHitInfo(right)
{
}

//_____________________________________________________________________________
inline const SCHHit&
SCHHit::operator =(const SCHHit& right)
{
  VHitInfo::operator =(right);
  return *this;
}

//_____________________________________________________________________________
// externally instanciated.
extern G4Allocator<SCHHit> SCHHitAllocator;

//_____________________________________________________________________________
inline void*
SCHHit::operator new(size_t)
{
  return SCHHitAllocator.MallocSingle();
}

//_____________________________________________________________________________
inline void
SCHHit::operator delete(void* aHit)
{
  SCHHitAllocator.FreeSingle(static_cast<SCHHit*>(aHit));
}

#endif

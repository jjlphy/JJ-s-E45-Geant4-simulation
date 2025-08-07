// -*- C++ -*-

#ifndef HTOF_HIT_HH
#define HTOF_HIT_HH

#include <G4Allocator.hh>
#include <G4String.hh>
#include <G4THitsCollection.hh>
#include <G4ThreeVector.hh>
#include <G4VHit.hh>

#include "VHitInfo.hh"

class G4Step;

//_____________________________________________________________________________
class HTOFHit : public G4VHit, public VHitInfo
{
public:
  HTOFHit(const G4String& name, G4Step* step);
  virtual ~HTOFHit();

  HTOFHit(const HTOFHit& right);
  const HTOFHit& operator=(const HTOFHit& right);

  void* operator new(size_t);
  void operator delete(void* aHit);

public:
  virtual void Draw();
  virtual void Print();
};

//_____________________________________________________________________________
inline
HTOFHit::HTOFHit(const HTOFHit& right)
  : G4VHit(right),
    VHitInfo(right)
{
}

//_____________________________________________________________________________
inline const HTOFHit&
HTOFHit::operator =(const HTOFHit& right)
{
  VHitInfo::operator =(right);
  return *this;
}

//_____________________________________________________________________________
// externally instanciated.
extern G4Allocator<HTOFHit> HTOFHitAllocator;

//_____________________________________________________________________________
inline void*
HTOFHit::operator new(size_t)
{
  return HTOFHitAllocator.MallocSingle();
}

//_____________________________________________________________________________
inline void
HTOFHit::operator delete(void* aHit)
{
  HTOFHitAllocator.FreeSingle(static_cast<HTOFHit*>(aHit));
}

#endif

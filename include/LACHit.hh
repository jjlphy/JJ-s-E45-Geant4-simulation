// -*- C++ -*-

#ifndef LAC_HIT_HH
#define LAC_HIT_HH

#include <G4Allocator.hh>
#include <G4String.hh>
#include <G4THitsCollection.hh>
#include <G4ThreeVector.hh>
#include <G4VHit.hh>

#include "VHitInfo.hh"

class G4Step;

//_____________________________________________________________________________
class LACHit : public G4VHit, public VHitInfo
{
public:
  LACHit(const G4String& name, G4Step* step);
  virtual ~LACHit();

  LACHit(const LACHit& right);
  const LACHit& operator=(const LACHit& right);

  void* operator new(size_t);
  void operator delete(void* aHit);

public:
  virtual void Draw();
  virtual void Print();
};

//_____________________________________________________________________________
inline
LACHit::LACHit(const LACHit& right)
  : G4VHit(right),
    VHitInfo(right)
{
}

//_____________________________________________________________________________
inline const LACHit&
LACHit::operator =(const LACHit& right)
{
  VHitInfo::operator =(right);
  return *this;
}

//_____________________________________________________________________________
// externally instanciated.
extern G4Allocator<LACHit> LACHitAllocator;

//_____________________________________________________________________________
inline void*
LACHit::operator new(size_t)
{
  void* aHit = (void*)LACHitAllocator.MallocSingle();
  return aHit;
}

//_____________________________________________________________________________
inline void
LACHit::operator delete(void* aHit)
{
  LACHitAllocator.FreeSingle((LACHit*) aHit);
}

#endif

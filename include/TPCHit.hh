// -*- C++ -*-

#ifndef TPC_HIT_HH
#define TPC_HIT_HH

#include <G4Allocator.hh>
#include <G4String.hh>
#include <G4VHit.hh>
#include <G4THitsCollection.hh>
#include <G4ThreeVector.hh>
#include <G4VHit.hh>

#include "VHitInfo.hh"

class G4Step;

//_____________________________________________________________________________
class TPCHit : public G4VHit, public VHitInfo
{
public:
  TPCHit(const G4String& name, G4Step* step);
  virtual ~TPCHit();

  TPCHit(const TPCHit& right);
  const TPCHit& operator=(const TPCHit& right);

  void* operator new(size_t);
  void operator delete(void* aHit);

public:
  virtual void Draw();
  virtual void Print();
};

//_____________________________________________________________________________
inline
TPCHit::TPCHit(const TPCHit& right)
  : G4VHit(right),
    VHitInfo(right)
{
}

//_____________________________________________________________________________
inline const TPCHit&
TPCHit::operator=(const TPCHit& right)
{
  VHitInfo::operator =(right);
  return *this;
}

//_____________________________________________________________________________
// externally instanciated.
extern G4Allocator<TPCHit> TPCHitAllocator;

//_____________________________________________________________________________
inline void* TPCHit::operator new(size_t)
{
  void* aHit = (void*)TPCHitAllocator.MallocSingle();
  return aHit;
}

//_____________________________________________________________________________
inline void TPCHit::operator delete(void* aHit)
{
  TPCHitAllocator.FreeSingle((TPCHit*) aHit);
}

#endif

// -*- C++ -*-

#ifndef WC_HIT_HH
#define WC_HIT_HH

#include <G4Allocator.hh>
#include <G4String.hh>
#include <G4THitsCollection.hh>
#include <G4ThreeVector.hh>
#include <G4VHit.hh>

#include "VHitInfo.hh"

class G4Step;

//_____________________________________________________________________________
class WCHit : public G4VHit, public VHitInfo
{
public:
  WCHit(const G4String& name, G4Step* step);
  virtual ~WCHit();

  WCHit(const WCHit& right);
  const WCHit& operator=(const WCHit& right);

  void* operator new(size_t);
  void operator delete(void* aHit);

public:
  virtual void Draw();
  virtual void Print();
};

//_____________________________________________________________________________
inline
WCHit::WCHit(const WCHit& right)
  : G4VHit(right),
    VHitInfo(right)
{
}

//_____________________________________________________________________________
inline const WCHit&
WCHit::operator =(const WCHit& right)
{
  VHitInfo::operator =(right);
  return *this;
}

//_____________________________________________________________________________
// externally instanciated.
extern G4Allocator<WCHit> WCHitAllocator;

//_____________________________________________________________________________
inline void*
WCHit::operator new(size_t)
{
  void* aHit = (void*)WCHitAllocator.MallocSingle();
  return aHit;
}

//_____________________________________________________________________________
inline void
WCHit::operator delete(void* aHit)
{
  WCHitAllocator.FreeSingle((WCHit*) aHit);
}

#endif

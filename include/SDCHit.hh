// -*- C++ -*-

#ifndef SDC_HIT_HH
#define SDC_HIT_HH

#include <G4Allocator.hh>
#include <G4String.hh>
#include <G4THitsCollection.hh>
#include <G4ThreeVector.hh>
#include <G4VHit.hh>

#include "VHitInfo.hh"

class G4Step;

//_____________________________________________________________________________
class SDCHit : public G4VHit, public VHitInfo
{
public:
  SDCHit(const G4String& name, G4Step* step);
  virtual ~SDCHit();

  SDCHit(const SDCHit& right);
  const SDCHit& operator=(const SDCHit& right);

  void* operator new(size_t);
  void operator delete(void* aHit);

public:
  virtual void Draw();
  virtual void Print();
};

//_____________________________________________________________________________
inline
SDCHit::SDCHit(const SDCHit& right)
  : G4VHit(right),
    VHitInfo(right)
{
}

//_____________________________________________________________________________
inline const SDCHit&
SDCHit::operator =(const SDCHit& right)
{
  VHitInfo::operator =(right);
  return *this;
}

//_____________________________________________________________________________
// externally instanciated.
extern G4Allocator<SDCHit> SDCHitAllocator;

//_____________________________________________________________________________
inline void*
SDCHit::operator new(size_t)
{
  return SDCHitAllocator.MallocSingle();
}

//_____________________________________________________________________________
inline void
SDCHit::operator delete(void* aHit)
{
  SDCHitAllocator.FreeSingle(static_cast<SDCHit*>(aHit));
}

#endif

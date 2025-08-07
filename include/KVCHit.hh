// -*- C++ -*-

#ifndef KVC_HIT_HH
#define KVC_HIT_HH

#include <G4Allocator.hh>
#include <G4String.hh>
#include <G4THitsCollection.hh>
#include <G4ThreeVector.hh>
#include <G4VHit.hh>

#include "VHitInfo.hh"

class G4Step;

//_____________________________________________________________________________
class KVCHit : public G4VHit, public VHitInfo
{
public:
  KVCHit(const G4String& name, G4Step* step);
  virtual ~KVCHit();

  KVCHit(const KVCHit& right);
  const KVCHit& operator=(const KVCHit& right);

  void* operator new(size_t);
  void operator delete(void* aHit);

public:
  virtual void Draw();
  virtual void Print();
};

//_____________________________________________________________________________
inline
KVCHit::KVCHit(const KVCHit& right)
  : G4VHit(right),
    VHitInfo(right)
{
}

//_____________________________________________________________________________
inline const KVCHit&
KVCHit::operator =(const KVCHit& right)
{
  VHitInfo::operator =(right);
  return *this;
}

//_____________________________________________________________________________
// externally instanciated.
extern G4Allocator<KVCHit> KVCHitAllocator;

//_____________________________________________________________________________
inline void*
KVCHit::operator new(size_t)
{
  void* aHit = (void*)KVCHitAllocator.MallocSingle();
  return aHit;
}

//_____________________________________________________________________________
inline void
KVCHit::operator delete(void* aHit)
{
  KVCHitAllocator.FreeSingle((KVCHit*) aHit);
}

#endif

// -*- C++ -*-

#ifndef TARGET_HIT_HH
#define TARGET_HIT_HH

#include <G4Allocator.hh>
#include <G4String.hh>
#include <G4THitsCollection.hh>
#include <G4ThreeVector.hh>
#include <G4VHit.hh>

#include "VHitInfo.hh"

class G4Step;

//_____________________________________________________________________________
class TargetHit : public G4VHit, public VHitInfo
{
public:
  TargetHit(const G4String& name, G4Step* step, G4bool use_center_point=false);
  virtual ~TargetHit();

  TargetHit(const TargetHit& right);
  const TargetHit& operator=(const TargetHit& right);

  void* operator new(size_t);
  void operator delete(void* aHit);

public:
  virtual void Draw();
  virtual void Print();
};

//_____________________________________________________________________________
inline
TargetHit::TargetHit(const TargetHit& right)
  : G4VHit(right),
    VHitInfo(right)
{
}

//_____________________________________________________________________________
inline const TargetHit&
TargetHit::operator =(const TargetHit& right)
{
  VHitInfo::operator =(right);
  return *this;
}

//_____________________________________________________________________________
// externally instanciated.
extern G4Allocator<TargetHit> TargetHitAllocator;

//_____________________________________________________________________________
inline void*
TargetHit::operator new(size_t)
{
  return TargetHitAllocator.MallocSingle();
}

//_____________________________________________________________________________
inline void
TargetHit::operator delete(void* aHit)
{
  TargetHitAllocator.FreeSingle(static_cast<TargetHit*>(aHit));
}

#endif

// -*- C++ -*-

#ifndef FTOF_HIT_HH
#define FTOF_HIT_HH

#include <G4Allocator.hh>
#include <G4String.hh>
#include <G4THitsCollection.hh>
#include <G4ThreeVector.hh>
#include <G4VHit.hh>

#include "VHitInfo.hh"

class G4Step;

//_____________________________________________________________________________
class FTOFHit : public G4VHit, public VHitInfo
{
public:
  FTOFHit(const G4String& name, G4Step* step);
  virtual ~FTOFHit();

  FTOFHit(const FTOFHit& right);
  const FTOFHit& operator=(const FTOFHit& right);

  void* operator new(size_t);
  void operator delete(void* aHit);

public:
  virtual void Draw();
  virtual void Print();
};

//_____________________________________________________________________________
inline
FTOFHit::FTOFHit(const FTOFHit& right)
  : G4VHit(right),
    VHitInfo(right)
{
}

//_____________________________________________________________________________
inline const FTOFHit&
FTOFHit::operator =(const FTOFHit& right)
{
  VHitInfo::operator =(right);
  return *this;
}

//_____________________________________________________________________________
// externally instanciated.
extern G4Allocator<FTOFHit> FTOFHitAllocator;

//_____________________________________________________________________________
inline void*
FTOFHit::operator new(size_t)
{
  return FTOFHitAllocator.MallocSingle();
}

//_____________________________________________________________________________
inline void
FTOFHit::operator delete(void* aHit)
{
  FTOFHitAllocator.FreeSingle(static_cast<FTOFHit*>(aHit));
}

#endif

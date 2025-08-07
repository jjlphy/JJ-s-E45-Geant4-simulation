// -*- C++ -*-
#ifndef T0_HIT_HH
#define T0_HIT_HH

#include <G4Allocator.hh>
#include <G4String.hh>
#include <G4THitsCollection.hh>
#include <G4ThreeVector.hh>
#include <G4VHit.hh>

#include "VHitInfo.hh"

class G4Step;

class T0Hit : public G4VHit, public VHitInfo
{
public:
  T0Hit(const G4String& name, G4Step* step);
  virtual ~T0Hit();

  T0Hit(const T0Hit& right);
  const T0Hit& operator=(const T0Hit& right);

  void* operator new(size_t);
  void operator delete(void* aHit);

  virtual void Draw();
  virtual void Print();
};

extern G4Allocator<T0Hit> T0HitAllocator;

inline T0Hit::T0Hit(const T0Hit& right) : G4VHit(right), VHitInfo(right) {}
inline const T0Hit& T0Hit::operator=(const T0Hit& right)
{ VHitInfo::operator=(right); return *this; }

inline void* T0Hit::operator new(size_t) { return T0HitAllocator.MallocSingle(); }
inline void T0Hit::operator delete(void* aHit) { T0HitAllocator.FreeSingle((T0Hit*)aHit); }

#endif

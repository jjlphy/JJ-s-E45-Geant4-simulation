// -*- C++ -*-

#ifndef PHYSICS_LIST_HH
#define PHYSICS_LIST_HH 1

#include <globals.hh>
#include <G4VModularPhysicsList.hh>

//_____________________________________________________________________________
class PhysicsList : public G4VModularPhysicsList
{
public:
  static G4String ClassName();
  PhysicsList(G4int verbose_level=1);
  virtual ~PhysicsList() = default;

  PhysicsList(const PhysicsList&) = delete;
  PhysicsList& operator =(const PhysicsList&) = delete;

protected:
  virtual void ConstructParticle();
  virtual void ConstructProcess();
  virtual void SetCuts();

private:
  void AddIonGasProcess();
};

//_____________________________________________________________________________
inline G4String
PhysicsList::ClassName()
{
  static G4String s_name("PhysicsList");
  return s_name;
}

#endif

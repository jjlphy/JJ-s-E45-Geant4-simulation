// -*- C++ -*-

#ifndef TPC_EDEP_SD_HH
#define TPC_EDEP_SD_HH

#include <G4VSensitiveDetector.hh>

#include "TPCEdepHit.hh"

class G4HCofThisEvent;
class G4Step;
class G4TouchableHistory;

//_____________________________________________________________________________
class TPCEdepSD : public G4VSensitiveDetector
{
public:
  static G4String ClassName( void );
  TPCEdepSD( const G4String& name );
  virtual ~TPCEdepSD( void );

private:
  G4THitsCollection<TPCEdepHit>* m_hits_collection;
  G4double select_plane;
  G4int num_plane;
  
  G4double DensityEffectCorrection(G4double betagamma, G4double* par);
  G4double TPCdEdx(G4double mass, G4double beta);
  G4double TPCdEdxSig(G4double mass, G4double mom);
  
public:
  G4int ntrk;
  virtual G4bool ProcessHits( G4Step* aStep, G4TouchableHistory* ROhist );
  virtual void   Initialize( G4HCofThisEvent* HCTE );
  virtual void   EndOfEvent( G4HCofThisEvent* HCTE );
  virtual void   DrawAll( void );
  virtual void   PrintAll( void );
};

//_____________________________________________________________________________
inline G4String
TPCEdepSD::ClassName( void )
{
  static G4String s_name("TPCEdepSD");
  return s_name;
}
#endif

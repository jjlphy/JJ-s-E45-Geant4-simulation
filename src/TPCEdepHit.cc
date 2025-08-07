// -*- C++ -*-

#include "TPCEdepHit.hh"

#include <CLHEP/Units/PhysicalConstants.h>

G4Allocator<TPCEdepHit> TPCEdepHitAllocator;

//_____________________________________________________________________________
TPCEdepHit::TPCEdepHit( G4ThreeVector& axyz, G4ThreeVector& apxyz, G4double t,
		      G4int tid, G4int pid,
		      G4int ilay, G4int irow, G4double Beta, G4double Ed,
		      G4int Parentid,
		      G4double tlength, G4double Mass, G4int Charge,
		      G4ThreeVector& Vtxpos, G4ThreeVector& Vtxmom,
		      G4double avtxene, G4double slength,
		      G4int Parentid_pid )
  : xyz(axyz), pxyz(apxyz), tof(t), trackID(tid), particleID(pid),
    iLay(ilay),iRow(irow), beta(Beta),edep(Ed), parentID(Parentid),
    Length(tlength), mass(Mass), charge(Charge), vtxmome(Vtxmom),
    vtxposi(Vtxpos), vtxene(avtxene), SLength(slength),
    parentID_pid(Parentid_pid)
{
}

//_____________________________________________________________________________
TPCEdepHit::~TPCEdepHit( void )
{
}

//_____________________________________________________________________________
void
TPCEdepHit::Draw( void )
{
}

//_____________________________________________________________________________
void
TPCEdepHit::Print( void )
{
  G4cout << "Hit in Counter:" << xyz*(1./CLHEP::cm) << " cm, "
	 << tof/CLHEP::ns << " ns" << G4endl;
}

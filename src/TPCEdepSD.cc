// -*- C++ -*-

#include "TPCEdepSD.hh"

#include <G4VPhysicalVolume.hh>
#include <G4Step.hh>
#include <G4Track.hh>
#include <G4VTouchable.hh>
#include <G4TouchableHistory.hh>
#include <G4PhysicalConstants.hh>
#include <G4SystemOfUnits.hh>
#include <G4DynamicParticle.hh>
#include <G4DecayProducts.hh>
#include <G4PhysicsLogVector.hh>
#include <G4ParticleChangeForDecay.hh>
#include <G4DecayProcessType.hh>
#include <Randomize.hh>

#include "ConfMan.hh"
#include "FuncName.hh"
#include "TPCEdepHit.hh"
#include "padHelper.hh"
#include "TPCPadHelper.hh"



namespace
{
  const auto& gConf = ConfMan::GetInstance();
}

//_____________________________________________________________________________
TPCEdepSD::TPCEdepSD( const G4String& name )
  : G4VSensitiveDetector( name ),
    m_hits_collection()
{
  collectionName.insert("hit");
}

//_____________________________________________________________________________
TPCEdepSD::~TPCEdepSD( void )
{
}

//_____________________________________________________________________________
void
TPCEdepSD::Initialize( G4HCofThisEvent* HCTE )
{
  ntrk=0;
  m_hits_collection = new G4THitsCollection<TPCEdepHit>( SensitiveDetectorName,
						     collectionName[0] );
  // push H.C. to "Hit Collection of This Event"
  G4int hcid = GetCollectionID(0);
  HCTE->AddHitsCollection( hcid, m_hits_collection );
}

//_____________________________________________________________________________
G4bool
TPCEdepSD::ProcessHits( G4Step* aStep, G4TouchableHistory* /* ROhist */ )
{

  const G4StepPoint* preStepPoint= aStep-> GetPreStepPoint();
  const G4StepPoint* postStepPoint= aStep-> GetPostStepPoint();
  const G4Track* aTrack = aStep->GetTrack();
  G4int copyNo = preStepPoint -> GetPhysicalVolume()->GetCopyNo();
  //std::cout<<"copyNo : "<<copyNo<<std::endl;
  G4int copyNo_post = postStepPoint ->GetPhysicalVolume()->GetCopyNo();
  //if(copyNo!=copyNo_post) return false;

  
  //if(preStepPoint-> GetStepStatus() != fGeomBoundary) return false;
  //  if(preStepPoint-> GetStepStatus() == fGeomBoundary){
  G4String particleName;
  
	if(aStep-> GetTrack()-> GetDefinition()-> GetPDGCharge() == 0.)
    return false;
	
  particleName = aStep-> GetTrack()-> GetDefinition()-> GetParticleName();

  G4String particleType;
  particleType = aTrack->GetDefinition()->GetParticleType();

  G4TouchableHistory* theTouchable
    = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
  G4ThreeVector pos= preStepPoint-> GetPosition();
  G4double hitx=pos.getX();
  G4double hitz=pos.getZ();

    G4ThreeVector VertexPosition = aTrack->GetVertexPosition();
  G4ThreeVector VertexMomentum = aTrack->GetVertexMomentumDirection();
  G4double VertexEnergy = aTrack -> GetVertexKineticEnergy(); // Ek = sqrt(p^2+m^2)-m
  G4VPhysicalVolume* physVol = theTouchable->GetVolume();

  G4ThreeVector mom= preStepPoint-> GetMomentum();
  G4double tof= preStepPoint-> GetGlobalTime();
  G4double beta= preStepPoint-> GetBeta();
  G4int tid =  aStep-> GetTrack()-> GetTrackID();
  G4int pid =  aStep-> GetTrack()-> GetDefinition() -> GetPDGEncoding();
  //  if(abs(pid)==11)return false; //remove e+, e-
  
  G4double mass = aStep -> GetTrack()->GetDynamicParticle()->GetMass();
  //  G4cout<<mass<<G4endl;
  G4int charge = aStep-> GetTrack()-> GetDefinition()-> GetPDGCharge();
  //  G4double edep = aStep->GetTotalEnergyDeposit();
  //  G4double edep = aStep->GetNonIonizingEnergyDeposit();
  //  G4double deltaep = aStep->GetDeltaEnergy();
  G4double tlength = aStep->GetTrack()-> GetTrackLength();
  G4double slength = aStep->GetTrack()-> GetStepLength();
  // G4double length = aStep-> GetStepLength();
  G4int parentID =  aStep-> GetTrack()-> GetParentID();

  G4int parentID_pid = aStep-> GetTrack()->GetDynamicParticle()->GetDefinition()->GetPDGEncoding();

  G4int iLay_copyNo=copyNo - 2000;
  //G4int iPad = padHelper::findPadID(hitz, hitx);
  //G4int iLay= padHelper::getLayerID(iPad);
  G4int iPad = TPCPadHelper::FindPadID(hitz, hitx);
  //G4int iLay = TPCPadHelper::GetLayerID(iPad);
  G4int iLay = iLay_copyNo;
  //G4int iRow= padHelper::getRowID(iPad);
  G4int iRow = TPCPadHelper::GetRowID(iPad);
  
	G4double PadLen = padHelper::getLength(iLay);
	G4ThreeVector PadPos(hitx,0,hitz + 143); 
	G4ThreeVector MomT(mom.x(),0,mom.z());	
	G4double alpha = PadPos.theta()-MomT.theta();
	G4double PathT = PadLen * 1./cos(alpha);
	G4double Pitch = mom.y()/MomT.mag();
	G4double Path = PathT * sqrt(1+Pitch*Pitch);
	slength = Path;

	/*
	G4double edepMean =TPCdEdx(mass,beta)*Path; 
	//	G4double IonEn = (0.9 * 188 + 0.1 * 41.7)*eV;
//	G4double edepSig = sqrt(edepMean / IonEn ) * IonEn;
	G4double edepSig =TPCdEdxSig(mass,mom.mag())*Path;
	if(edepSig / edepMean < 0.01 or edepSig/edepMean > 0.5)edepSig = 0.2*edepMean;
	G4double edep = G4RandGauss::shoot(edepMean,edepSig);
	if(edep <0.1* edepMean)edep = 0.1*edepMean;
	*/

	G4double conversion_factor = 11073.3; //MeV/cm -> ADC/mm 
	G4double cmTomm = 10;
	G4double edep = aStep->GetTotalEnergyDeposit() * conversion_factor * cmTomm; 


	//if(iLay_copyNo==4 && pos.getY()>-1 * 97./2.)return false;
	
#ifdef DEBUG
  
  //for test 
  G4double radius = sqrt( hitx*hitx + (hitz+143.)*(hitz+143.));
  TVector3 Point = padHelper::getPoint(iPad);
  //G4int iPad_re = padHelper::findPadID(Point.z(), Point.x());
  G4int iPad_re = iPad;
  /*
  G4cout<<"hitx = "<< hitx
   	<<", pointx "<< Point.x()
   	<<", hitz =" << hitz
  	<<", pointz "<< Point.z()
   	<<"Pre :  radius = "<< radius
   	<<", iPad ="<<iPad
	<<", iPad_re ="<<iPad_re
   	<<", iLay_copyNo = " <<iLay_copyNo
   	<< ", iLay = "<<iLay
        <<", Pid = "<<pid
	<<", Edep = "<<edep<<G4endl;
  */

  /*G4cout<<"dx ="<<hitx-Point.x()
    <<", dz ="<<hitz-Point.z()<<std::endl;*/

    //check post step point

  G4ThreeVector pos_post= postStepPoint-> GetPosition();
  G4double hitx_post=pos_post.getX();
  G4double hitz_post=pos_post.getZ();
  G4int iPad_post = TPCPadHelper::FindPadID(hitz_post,hitx_post);
  G4int iLay_post = TPCPadHelper::GetLayerID(iPad_post);
    /*
  G4cout<<"Post : radius = "<<sqrt( hitx_post*hitx_post + (hitz_post+143.)*(hitz_post+143.))
	<<", iPad ="<<iPad_post
   	<<", iLay_copyNo = " <<copyNo_post-2000
   	<< ", iLay = "<<iLay_post
	<<", Edep = "<<edep<<G4endl;
    */

  if(iLay<0 || iLay > 31 || iLay!=iLay_copyNo)G4cout<<"pre strange layer"<<G4endl;
  if(iLay_post<0 || iLay_post > 31)G4cout<<"post strange layer"<<G4endl;

#endif


  
  G4String name = physVol->GetName();
  
  



  if(name=="TPC_PV"){
    name="EdepPV-1";
  }

  sscanf(name,"EdepPV%d",&iLay_copyNo);
  TPCEdepHit* ahit= new TPCEdepHit(pos, mom, tof, tid, pid, iLay, iRow, beta, edep,parentID,tlength, mass, charge, VertexPosition, VertexMomentum, VertexEnergy,slength, parentID_pid );

  m_hits_collection-> insert(ahit);
  return true;

}

//_____________________________________________________________________________
void
TPCEdepSD::EndOfEvent(G4HCofThisEvent* /* HCTE */ )
{
}

//_____________________________________________________________________________
void
TPCEdepSD::DrawAll( void )
{
}

G4double
TPCEdepSD::TPCdEdx(Double_t mass/*MeV/c2*/, Double_t beta){

  Double_t rho=0.; //[g cm-3]
  Double_t ZoverA=0.; //[mol g-1]
  Double_t I=0.; //[eV]
  Double_t density_effect_par[6]={0.}; //Sternheimer’s parameterization
  //P10  
	rho = TMath::Power(10.,-3)*(0.9*1.662 + 0.1*0.6672);
	ZoverA = 17.2/37.6;
	I = 0.9*188.0 + 0.1*41.7;
	density_effect_par[0] = 0.9*0.19714 + 0.1*0.09253;
	density_effect_par[1] = 0.9*2.9618 + 0.1*3.6257;
	density_effect_par[2] = 0.9*1.7635 + 0.1*1.6263;
	density_effect_par[3] = 0.9*4.4855 + 0.1*3.9716;
	density_effect_par[4] = 0.9*11.9480 + 0.1*9.5243;
	density_effect_par[5] = 0.;

  Double_t Z = 1.;
  Double_t me = 0.5109989461; //[MeV]
  Double_t K = 0.307075; //[MeV cm2 mol-1]
  Double_t constant = rho*K*ZoverA; //[MeV cm-1]
  constant = constant;
  Double_t I2 = I*I; //Mean excitaion energy [eV]
  Double_t beta2 = beta*beta;
  Double_t gamma2 = 1./(1.-beta2);
  Double_t MeVToeV = TMath::Power(10.,6);
  Double_t Wmax = 2*me*beta2*gamma2/((me/mass+1.)*(me/mass+1.)+2*(me/mass)*(TMath::Sqrt(gamma2)-1));
  Double_t delta = DensityEffectCorrection(TMath::Sqrt(beta2*gamma2), density_effect_par);
  Double_t dedx = constant*Z*Z/beta2*(0.5*TMath::Log(2*me*beta2*gamma2*Wmax*MeVToeV*MeVToeV/I2) - beta2 - 0.5*delta);
  
	
  G4double conversion_factor = 11073.3;
  return conversion_factor*dedx;

}
double
TPCEdepSD::TPCdEdxSig(Double_t mass/*MeV/c2*/, Double_t mom){
	double par[3]={0,0,0} ;
	if(mass < 0.2){//pion;
		par[0] = 7.792;
		par[1] =	-8.704;
		par[2] = 4.477;
	}
	else if(mass > 0.7){//proton;
		par[0] = 33.92;
		par[1] = -26.24;
		par[2] = 6.259; 
	}
	double value = par[0]+par[1]*mom+par[2]*mom*mom;
	return value;
}
G4double
TPCEdepSD::DensityEffectCorrection(Double_t betagamma, Double_t *par){

    //reference : Sternheimer’s parameterizatio(PDG)
    //notation : par[0] : a, par[1] : k, par[2] : x0, par[3] : x1, par[4] : _C, par[5] : delta0
    Double_t constant = 2*TMath::Log(10);
    Double_t delta = 0.;
    Double_t X = log10(betagamma);
    if(X<=par[2]) delta = par[5]*TMath::Power(10., 2*(X - par[2]));
    else if(par[2]<X && X<par[3]) delta = constant*X - par[4] + par[0]*pow((par[3] - X), par[1]);
    else if(X>=par[3]) delta = constant*X - par[4];

  return delta;

}
//_____________________________________________________________________________
void
TPCEdepSD::PrintAll( void )
{
  m_hits_collection-> PrintAllHits();
}


// -*- C++ -*-

#include "PrimaryGeneratorAction.hh"

#include <G4Event.hh>
#include <G4IonTable.hh>
#include <G4ParticleGun.hh>
#include <G4ParticleTable.hh>
#include <G4ParticleDefinition.hh>
#include <G4UImanager.hh>
#include <G4IonConstructor.hh>
#include <Randomize.hh>

#include "AnaManager.hh"
#include "BeamMan.hh"
#include "ConfMan.hh"
#include "DCGeomMan.hh"
#include "DetSizeMan.hh"
#include "FuncName.hh"
#include "JamMan.hh"
#include "Kinema3Resonance.hh"
#include "KinemaHResonance.hh"
#include "Kinema3Body.hh"
#include "Kinema4Body.hh"
#include "KinemaHybrid.hh"
#include "KinemaHweak.hh"
#include "KinemaFermi.hh"
#include "KinemaKstar.hh"
#include "Kinematics.hh"

namespace
{
auto& gAnaMan = AnaManager::GetInstance();
const auto& gBeam = BeamMan::GetInstance();
const auto& gConf = ConfMan::GetInstance();
const auto& gGeom = DCGeomMan::GetInstance();
const auto& gSize = DetSizeMan::GetInstance();
const auto& gJam  = JamMan::GetInstance();
}

//_____________________________________________________________________________
void
PrimaryGeneratorAction::GenerateE07Study(G4Event* anEvent)
{

  double pbm[4]={-9999.9999}, pka[4]={-9999.9999}, vtx[3]={-9999.9999};
  // double pL1[4]={-9999.9999}, pL2[4]={-9999.9999};
  //,xL1[3]={-9999.9999}, xL2[3]={-9999.9999};
  //  double pp1[4]={-9999.9999}, ppi1[4]={-9999.9999}, pp2[4]={-9999.9999}, ppi2[4]={-9999.9999};

  char fnam[30] = "inc64cu.dat";

  FILE *fp;

  if ((fp = fopen(fnam,"r")) == NULL){
    fprintf(stderr, "inc64cu.dat:: Cannot open file: %s\n", fnam);
    exit(-1);
  }

  double data[100]={-9999.9999};
  int check=0.;
  // 22478 lines --> remove NaN
  int ran = G4RandFlat::shoot(1.,22478.);

  while(1){
    fscanf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n"
	   ,&data[0],&data[1],&data[2],&data[3],&data[4]
	   ,&data[5],&data[6],&data[7],&data[8],&data[9]
	   ,&data[10],&data[11],&data[12],&data[13],&data[14]
	   ,&data[15],&data[16],&data[17],&data[18],&data[19]
	   ,&data[20],&data[21],&data[22],&data[23],&data[24]
	   ,&data[25],&data[26],&data[27],&data[28],&data[29]
	   ,&data[30],&data[31],&data[32],&data[33],&data[34]
	   ,&data[35],&data[36],&data[37],&data[38],&data[39]
	   ,&data[40]);
    check=check+1.;
    if(check==ran) break;
    //    if(res==EOF) break;
  }

  for(int ii=0;ii<4;ii++){
    pbm[ii]=data[ii]/1000.;
    pka[ii]=data[ii+4]/1000.;
    // pL1[ii]=data[ii+11]/1000.;
    // pL2[ii]=data[ii+15]/1000.;
  }

  for(int ii=0;ii<3;ii++){
    vtx[ii]=data[ii+8];
  }

  fclose(fp);

  //  G4cout<<"m_target_size.z() :"<<atof(env_m_target_size.z().c_str())<<G4endl;
  vtx[2]= G4RandFlat::shoot(100.-m_target_size.z()/2,100.+m_target_size.z()/2);

  // G4double Energy_beam = pbm[3];
  G4ThreeVector momentumBeam(pbm[0],pbm[1],pbm[2]);
  //  G4ThreeVector vertexPos(vtx[0]*CLHEP::mm,vtx[1]*CLHEP::mm, vtx[2]*CLHEP::mm);
  G4ThreeVector vertexPos(0.*CLHEP::mm,0.*CLHEP::mm, vtx[2]*CLHEP::mm);

  G4double Energy_ka = pka[3]; //total energy sqrt(pp^2+rmk^2)??
  // G4double Energy_L1 = pL1[3]; //lambda momenta ?? --> 4 momomenta sqrt(p^2+m^2)
  // G4double Energy_L2 = pL2[3]; //lambda momenta ?? --> 4 momomenta sqrt(p^2+m^2)

  // G4double beta= pbm[2]/(2*0.93827203+pbm[3]);
  // G4double momcmk[4] = {};
  // G4double momk[4];
  // momk[0]=pka[0];
  // momk[1]=pka[1];
  // momk[2]=pka[2];
  // momk[3]=sqrt(pka[0]*pka[0]+pka[1]*pka[1]+pka[2]*pka[2]+0.493677*0.493677);
  // G4double momkpp=sqrt(pow(momk[0],2)+pow(momk[1],2)+pow(momk[2],2));
  // G4double test = lorentz(momk,beta,momcmk);
  // G4double momcmkpp=sqrt(pow(momcmk[0],2)+pow(momcmk[1],2)+pow(momcmk[2],2));

  // cmk->Fill(acos(momcmk[2]/momcmkpp)*180/3.141592654);
  // coscmk->Fill(momcmk[2]/momcmkpp);

  // labk->Fill(acos(momk[2]/momkpp)*180/3.141592654);
  // coslabk->Fill(momk[2]/momkpp);

  // G4double shphi=(atan2(momk[1],momk[0]))/3.141592654*180.;
  //phik->Fill(shphi);

  //  G4cout<<pka[1]<<G4endl;
  // ---- K+ -------------------
  G4ThreeVector momentumKaonPlus(pka[0], pka[1], pka[2]);
  m_particle_gun->SetParticleDefinition(m_KaonPlus);
  m_particle_gun->SetParticleMomentumDirection(momentumKaonPlus);
  m_particle_gun->SetParticleEnergy((Energy_ka - m_KaonPlus->GetPDGMass()/CLHEP::GeV)*CLHEP::GeV);
  m_particle_gun->SetParticlePosition(vertexPos);
  m_particle_gun->GeneratePrimaryVertex(anEvent);

  // ---- K-(beam)  -------------
  //  G4ThreeVector momentumLambda1(pbm[0], pbm[1], pbm[2]);
  //  m_particle_gun->SetParticleDefinition(m_KaonMinus);
  //  m_particle_gun->SetParticleMomentumDirection(momentumBeam);
  //  m_particle_gun->SetParticleEnergy((Energy_beam - m_KaonMinus->GetPDGMass()/CLHEP::GeV)*CLHEP::GeV);
  //  m_particle_gun->SetParticlePosition(vertexPos);
  //  m_particle_gun->GeneratePrimaryVertex(anEvent);

  // gAnaMan.SetPrimaryParticle(0,pka[0],pka[1],pka[2],m_KaonPlus->GetPDGMass()/CLHEP::GeV);
  //  gAnaMan.SetPrimaryParticle(1,pbm[0],pbm[1],pbm[2],m_KaonMinus->GetPDGMass()/CLHEP::GeV);

  // gAnaMan.SetPrimaryVertex(0,0.,0.,vtx[2]);
  //  gAnaMan.SetPrimaryVertex(1,0.,0.,vtx[2]);

  // G4double rmk=0.493677;
  // G4double misskp = miss1(pbm,rmk,momk);//pbeam, mass, momk

  //  missk->Fill(misskp);
  // G4double ptot=pow(pL1[0]+pL2[0],2)+pow(pL1[1]+pL2[1],2)+pow(pL1[2]+pL2[2],2);
  // G4double e1=Energy_L1;
  // G4double e2=Energy_L2;
  // G4double etot=pow(e1+e2,2);
  // G4double invm2=(etot-ptot);
  // if(invm2 > 0) invm=sqrt(invm2);
  // gen_im->Fill(invm);
}

//_____________________________________________________________________________
void
PrimaryGeneratorAction::GenerateE07StudyKp(G4Event* anEvent)
{

  double pbm[4]={-9999.9999}, pka[4]={-9999.9999}, vtx[3]={-9999.9999};
  // double pL1[4]={-9999.9999}, pL2[4]={-9999.9999};
  //,xL1[3]={-9999.9999}, xL2[3]={-9999.9999};
  //  double pp1[4]={-9999.9999}, ppi1[4]={-9999.9999}, pp2[4]={-9999.9999}, ppi2[4]={-9999.9999};

  char fnam[] = "inc64cu.dat";

  FILE *fp;

  if ((fp = fopen(fnam,"r")) == NULL){
    fprintf(stderr, "inc64cu.dat:: Cannot open file: %s\n", fnam);
    exit(-1);
  }

  double data[100]={-9999.9999};
  int check=0.;
  // 22478 lines --> remove NaN
  int ran = G4RandFlat::shoot(1.,22478.);

  while(1){
    fscanf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n"
	   ,&data[0],&data[1],&data[2],&data[3],&data[4]
	   ,&data[5],&data[6],&data[7],&data[8],&data[9]
	   ,&data[10],&data[11],&data[12],&data[13],&data[14]
	   ,&data[15],&data[16],&data[17],&data[18],&data[19]
	   ,&data[20],&data[21],&data[22],&data[23],&data[24]
	   ,&data[25],&data[26],&data[27],&data[28],&data[29]
	   ,&data[30],&data[31],&data[32],&data[33],&data[34]
	   ,&data[35],&data[36],&data[37],&data[38],&data[39]
	   ,&data[40]);
    check=check+1.;
    if(check==ran) {
      //      G4cout<<check<<G4endl;
      break;
    }
    //    if(res==EOF) break;
  }

  for(int ii=0;ii<4;ii++){
    pbm[ii]=data[ii]/1000.;
    pka[ii]=data[ii+4]/1000.;
    // pL1[ii]=data[ii+11]/1000.;
    // pL2[ii]=data[ii+15]/1000.;
  }

  for(int ii=0;ii<3;ii++){
    vtx[ii]=data[ii+8];
  }

  fclose(fp);

  G4double Energy_ka;
  //  G4cout<<"m_target_size.z() :"<<atof(env_m_target_size.z().c_str())<<G4endl;
  vtx[2]= G4RandFlat::shoot(100.-m_target_size.z()/2,100.+m_target_size.z()/2);
  // G4double Energy_beam = pbm[3];
  G4ThreeVector momentumBeam(pbm[0],pbm[1],pbm[2]);
  //  G4ThreeVector vertexPos(vtx[0]*CLHEP::mm,vtx[1]*CLHEP::mm, vtx[2]*CLHEP::mm);
  G4ThreeVector vertexPos(0.*CLHEP::mm,0.*CLHEP::mm, vtx[2]*CLHEP::mm);

  Energy_ka = pka[3]; //total energy sqrt(pp^2+rmk^2)??
  // ---- K+ -------------------
  G4ThreeVector momentumKaonPlus(pka[0], pka[1], pka[2]);
  m_particle_gun->SetParticleDefinition(m_KaonPlus);
  m_particle_gun->SetParticleMomentumDirection(momentumKaonPlus);
  m_particle_gun->SetParticleEnergy((Energy_ka - m_KaonPlus->GetPDGMass()/CLHEP::GeV)*CLHEP::GeV);
  m_particle_gun->SetParticlePosition(vertexPos);
  m_particle_gun->GeneratePrimaryVertex(anEvent);
  // gAnaMan.SetPrimaryParticle(0,pka[0],pka[1],pka[2],m_KaonPlus->GetPDGMass()/CLHEP::GeV);
  // gAnaMan.SetPrimaryVertex(0,0.,0.,vtx[2]);
}

//_____________________________________________________________________________
void
PrimaryGeneratorAction::GenerateE07StudyKpBeam(G4Event* anEvent)
{

  double pbm[4]={-9999.9999}, pka[4]={-9999.9999}, vtx[3]={-9999.9999};
  // double pL1[4]={-9999.9999}, pL2[4]={-9999.9999};
  //,xL1[3]={-9999.9999}, xL2[3]={-9999.9999};
  //  double pp1[4]={-9999.9999}, ppi1[4]={-9999.9999}, pp2[4]={-9999.9999}, ppi2[4]={-9999.9999};

  char fnam[] = "inc64cu.dat";

  FILE *fp;

  if ((fp = fopen(fnam,"r")) == NULL){
    fprintf(stderr, "inc64cu.dat:: Cannot open file: %s\n", fnam);
    exit(-1);
  }

  double data[100]={-9999.9999};
  int check=0.;
  // 22478 lines --> remove NaN
  int ran = G4RandFlat::shoot(1.,22478.);

  while(1){
    fscanf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n"
	   ,&data[0],&data[1],&data[2],&data[3],&data[4]
	   ,&data[5],&data[6],&data[7],&data[8],&data[9]
	   ,&data[10],&data[11],&data[12],&data[13],&data[14]
	   ,&data[15],&data[16],&data[17],&data[18],&data[19]
	   ,&data[20],&data[21],&data[22],&data[23],&data[24]
	   ,&data[25],&data[26],&data[27],&data[28],&data[29]
	   ,&data[30],&data[31],&data[32],&data[33],&data[34]
	   ,&data[35],&data[36],&data[37],&data[38],&data[39]
	   ,&data[40]);
    check=check+1.;
    if(check==ran) {
      //      G4cout<<check<<G4endl;
      break;
    }
    //    if(res==EOF) break;
  }

  for(int ii=0;ii<4;ii++){
    pbm[ii]=data[ii]/1000.;
    pka[ii]=data[ii+4]/1000.;
    // pL1[ii]=data[ii+11]/1000.;
    // pL2[ii]=data[ii+15]/1000.;
  }

  //  for(int ii=0;ii<3;ii++){
  //    vtx[ii]=data[ii+8];
  //  }
  fclose(fp);

  G4double Energy_ka;
  vtx[0] = G4RandFlat::shoot(-15.,15.)*CLHEP::mm;
  vtx[1] = G4RandFlat::shoot(-5.,5.)*CLHEP::mm;
  vtx[2]= G4RandFlat::shoot(100.-m_target_size.z()/2,100.+m_target_size.z()/2)*CLHEP::mm;

  // G4double Energy_beam = pbm[3];
  G4ThreeVector momentumBeam(pbm[0],pbm[1],pbm[2]);

  // gAnaMan.SetPrimaryBeam(pbm[0],pbm[1],pbm[2]);

  G4ThreeVector vertexPos(vtx[0],vtx[1],vtx[2]);

  Energy_ka = pka[3]; //total energy sqrt(pp^2+rmk^2)??
  // ---- K+ -------------------
  G4ThreeVector momentumKaonPlus(pka[0], pka[1], pka[2]);
  m_particle_gun->SetParticleDefinition(m_KaonPlus);
  m_particle_gun->SetParticleMomentumDirection(momentumKaonPlus);
  m_particle_gun->SetParticleEnergy((Energy_ka - m_KaonPlus->GetPDGMass()/CLHEP::GeV)*CLHEP::GeV);
  m_particle_gun->SetParticlePosition(vertexPos);
  m_particle_gun->GeneratePrimaryVertex(anEvent);
  // gAnaMan.SetPrimaryParticle(0,pka[0],pka[1],pka[2],m_KaonPlus->GetPDGMass()/CLHEP::GeV);
  // gAnaMan.SetPrimaryVertex(0,vtx[0],vtx[1],vtx[2]);

}

//_____________________________________________________________________________
void
PrimaryGeneratorAction::GenerateE07StudyAll(G4Event* anEvent)
{

  double pbm[4]={-9999.9999}, pka[4]={-9999.9999}, vtx[3]={-9999.9999};
  // double pL1[4]={-9999.9999}, pL2[4]={-9999.9999};
  //,xL1[3]={-9999.9999}, xL2[3]={-9999.9999};
  //  double pp1[4]={-9999.9999}, ppi1[4]={-9999.9999}, pp2[4]={-9999.9999}, ppi2[4]={-9999.9999};

  char fnam[] = "inc64cu.dat";

  FILE *fp;

  if ((fp = fopen(fnam,"r")) == NULL){
    fprintf(stderr, "inc64cu.dat:: Cannot open file: %s\n", fnam);
    exit(-1);
  }

  double data[100]={-9999.9999};
  int check=0.;
  // 22478 lines --> remove NaN
  int ran = G4RandFlat::shoot(1.,22478.);

  while(1){
    fscanf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n"
	   ,&data[0],&data[1],&data[2],&data[3],&data[4]
	   ,&data[5],&data[6],&data[7],&data[8],&data[9]
	   ,&data[10],&data[11],&data[12],&data[13],&data[14]
	   ,&data[15],&data[16],&data[17],&data[18],&data[19]
	   ,&data[20],&data[21],&data[22],&data[23],&data[24]
	   ,&data[25],&data[26],&data[27],&data[28],&data[29]
	   ,&data[30],&data[31],&data[32],&data[33],&data[34]
	   ,&data[35],&data[36],&data[37],&data[38],&data[39]
	   ,&data[40]);
    check=check+1.;
    if(check==ran) break;
    //    if(res==EOF) break;
  }

  for(int ii=0;ii<4;ii++){
    pbm[ii]=data[ii]/1000.;
    pka[ii]=data[ii+4]/1000.;
    // pL1[ii]=data[ii+11]/1000.;
    // pL2[ii]=data[ii+15]/1000.;
  }

  for(int ii=0;ii<3;ii++){
    vtx[ii]=data[ii+8];
  }

  fclose(fp);

  //  G4double Energy_p;
  G4double Energy_kp;
  //  G4double Energy_pip;
  //  G4double Energy_pin;
  //  G4double Energy_kn;

  //  G4cout<<"m_target_size.z() :"<<atof(env_m_target_size.z().c_str())<<G4endl;
  vtx[2]= G4RandFlat::shoot(100.-m_target_size.z()/2,100.+m_target_size.z()/2);

  // G4double Energy_beam = pbm[3];
  G4ThreeVector momentumBeam(pbm[0],pbm[1],pbm[2]);
  //  G4ThreeVector vertexPos(vtx[0]*CLHEP::mm,vtx[1]*CLHEP::mm, vtx[2]*CLHEP::mm);
  G4ThreeVector vertexPos(0.*CLHEP::mm,0.*CLHEP::mm, vtx[2]*CLHEP::mm);
  // G4double Energy_L1 = pL1[3]; //total energy sqrt(pp^2+rmk^2)??
  // G4double Energy_L2 = pL2[3]; //total energy sqrt(pp^2+rmk^2)??

  //  Energy_p = sqrt(pow(pka[0],2)+pow(pka[1],2)+pow(pka[2],2)+m_Proton->GetPDGMass()/CLHEP::GeV*m_Proton->GetPDGMass()/CLHEP::GeV); //total energy sqrt(pp^2+rmk^2)??
  Energy_kp = sqrt(pow(pka[0],2)+pow(pka[1],2)+pow(pka[2],2)+m_KaonPlus->GetPDGMass()/CLHEP::GeV*m_KaonPlus->GetPDGMass()/CLHEP::GeV); //total energy sqrt(pp^2+rmk^2)??
  //  Energy_pip = sqrt(pow(pka[0],2)+pow(pka[1],2)+pow(pka[2],2)+PionPlus->GetPDGMass()/CLHEP::GeV*PionPlus->GetPDGMass()/CLHEP::GeV); //total energy sqrt(pp^2+rmk^2)??
  //  Energy_pin = sqrt(pow(pka[0],2)+pow(pka[1],2)+pow(pka[2],2)+m_Proton->GetPDGMass()/CLHEP::GeV*m_Proton->GetPDGMass()/CLHEP::GeV); //total energy sqrt(pp^2+rmk^2)??
  //  Energy_kn = sqrt(pow(pka[0],2)+pow(pka[1],2)+pow(pka[2],2)+m_KaonMinus->GetPDGMass()/CLHEP::GeV*m_KaonMinus->GetPDGMass()/CLHEP::GeV); //total energy sqrt(pp^2+rmk^2)??

  // G4double beta= pbm[2]/(2*0.93827203+pbm[3]);
  // G4double momcmk[4]={0.};
  // G4double momk[4]={0.};
  // momk[0]=pka[0];
  // momk[1]=pka[1];
  // momk[2]=pka[2];
  // momk[3]=sqrt(pka[0]*pka[0]+pka[1]*pka[1]+pka[2]*pka[2]+0.493677*0.493677);

  // G4double momkpp = sqrt(pow(momk[0],2)+pow(momk[1],2)+pow(momk[2],2));

  // G4double test=lorentz(momk,beta,momcmk);
  // G4double momcmkpp=sqrt(pow(momcmk[0],2)+pow(momcmk[1],2)+pow(momcmk[2],2));

  // cmk->Fill(acos(momcmk[2]/momcmkpp)*180/3.141592654);
  // coscmk->Fill(momcmk[2]/momcmkpp);

  // labk->Fill(acos(momk[2]/momkpp)*180/3.141592654);
  // coslabk->Fill(momk[2]/momkpp);

  // G4double shphi=(atan2(momk[1],momk[0]))/3.141592654*180.;
  //phik->Fill(shphi);

  // ---- proton -------------------
  //  G4ThreeVector momentumKaonPlus(pka[0], pka[1], pka[2]);
  //  m_particle_gun->SetParticleDefinition(m_Proton);
  //  m_particle_gun->SetParticleMomentumDirection(momentumKaonPlus);
  //  m_particle_gun->SetParticleEnergy((Energy_p - m_Proton->GetPDGMass()/CLHEP::GeV)*CLHEP::GeV);
  //  m_particle_gun->SetParticlePosition(vertexPos);
  //  m_particle_gun->GeneratePrimaryVertex(anEvent);

  // ---- K+ -------------------
  G4ThreeVector momentumKaonPlus(pka[0], pka[1], pka[2]);
  m_particle_gun->SetParticleDefinition(m_KaonPlus);
  m_particle_gun->SetParticleMomentumDirection(momentumKaonPlus);
  m_particle_gun->SetParticleEnergy((Energy_kp - m_KaonPlus->GetPDGMass()/CLHEP::GeV)*CLHEP::GeV);
  m_particle_gun->SetParticlePosition(vertexPos);
  m_particle_gun->GeneratePrimaryVertex(anEvent);

  // ---- pi+ -------------------
  //  G4ThreeVector momentumKaonPlus(pka[0], pka[1], pka[2]);
  //  m_particle_gun->SetParticleDefinition(PionPlus);
  //  m_particle_gun->SetParticleMomentumDirection(momentumKaonPlus);
  //  m_particle_gun->SetParticleEnergy((Energy_pip - PionPlus->GetPDGMass()/CLHEP::GeV)*CLHEP::GeV);
  //  m_particle_gun->SetParticlePosition(vertexPos);
  //  m_particle_gun->GeneratePrimaryVertex(anEvent);
  /*
  // ---- pi- -------------------
  G4ThreeVector momentumKaonPlus(pka[0], pka[1], pka[2]);
  m_particle_gun->SetParticleDefinition(PionMinus);
  m_particle_gun->SetParticleMomentumDirection(momentumKaonPlus);
  m_particle_gun->SetParticleEnergy((Energy_pin - PionMinus->GetPDGMass()/CLHEP::GeV)*CLHEP::GeV);
  m_particle_gun->SetParticlePosition(vertexPos);
  m_particle_gun->GeneratePrimaryVertex(anEvent);

  // ---- K- -------------------
  G4ThreeVector momentumKaonPlus(pka[0], pka[1], pka[2]);
  m_particle_gun->SetParticleDefinition(m_KaonMinus);
  m_particle_gun->SetParticleMomentumDirection(momentumKaonPlus);
  m_particle_gun->SetParticleEnergy((Energy_kn - m_KaonMinus->GetPDGMass()/CLHEP::GeV)*CLHEP::GeV);
  m_particle_gun->SetParticlePosition(vertexPos);
  m_particle_gun->GeneratePrimaryVertex(anEvent);
  */

  //  // ---- K-(beam)  -------------
  //  G4ThreeVector momentumLambda1(pbm[0], pbm[1], pbm[2]);
  //  m_particle_gun->SetParticleDefinition(m_KaonMinus);
  //  m_particle_gun->SetParticleMomentumDirection(momentumBeam);
  //  m_particle_gun->SetParticleEnergy((Energy_beam - m_KaonMinus->GetPDGMass()/CLHEP::GeV)*CLHEP::GeV);
  //  m_particle_gun->SetParticlePosition(vertexPos);
  //  m_particle_gun->GeneratePrimaryVertex(anEvent);

  // gAnaMan.SetPrimaryParticle(0,pka[0],pka[1],pka[2],m_KaonPlus->GetPDGMass()/CLHEP::GeV);
  //  gAnaMan.SetPrimaryParticle(1,pka[0],pka[1],pka[2],m_KaonPlus->GetPDGMass()/CLHEP::GeV);
  //  gAnaMan.SetPrimaryParticle(2,pka[0],pka[1],pka[2],PionPlus->GetPDGMass()/CLHEP::GeV);

  // gAnaMan.SetPrimaryVertex(0,0.,0.,vtx[2]);
  //  gAnaMan.SetPrimaryVertex(1,0.,0.,vtx[2]);
  //  gAnaMan.SetPrimaryVertex(2,0.,0.,vtx[2]);

  // G4double rmk=0.493677;
  // G4double misskp = miss1(pbm,rmk,momk);//pbeam, mass, momk

  //missk->Fill(misskp);
  // G4double ptot=pow(pL1[0]+pL2[0],2)+pow(pL1[1]+pL2[1],2)+pow(pL1[2]+pL2[2],2);
  // G4double e1=Energy_L1;
  // G4double e2=Energy_L2;
  // G4double etot=pow(e1+e2,2);
  // G4double invm2=(etot-ptot);
  // if(invm2 > 0)
  //   invm=sqrt(invm2);
  //gen_im->Fill(invm);
}

//_____________________________________________________________________________
void
PrimaryGeneratorAction::GenerateE07StudyKnP(G4Event* anEvent)
{
  G4double mom[3];
  G4double Ebeam, pbeam;
  //  G4double rmk=0.493677;
  G4double pbm[4];
  G4double Energy_h, mom_h_x, mom_h_y, mom_h_z;
  G4double Energy_kp, mom_kp_x, mom_kp_y, mom_kp_z;

  //  G4double pimom=0.635+G4RandFlat::shoot()*(2.000-0.635);
  pbeam=1.8;
  // gAnaMan.SetPrimaryBeam(0,0,pbeam);
  Ebeam = sqrt(pbeam*pbeam+m_KaonMinus->GetPDGMass()/CLHEP::GeV*m_KaonMinus->GetPDGMass()/CLHEP::GeV);

  pbm[0]=0;
  pbm[1]=0;
  pbm[2]=pbeam;
  pbm[3]=Ebeam;
 up:
  KinemaHweak Hkinema(m_KaonMinus->GetPDGMass()/CLHEP::GeV,
                      m_Proton->GetPDGMass()/CLHEP::GeV,
                      m_KaonMinus->GetPDGMass()/CLHEP::GeV,
                      m_Proton->GetPDGMass()/CLHEP::GeV,
                      pbm[2], 0.0);

  Energy_kp=Hkinema.GetEnergy(3);
  // G4double momentum_kp = Hkinema.GetMomentum(3);
  Hkinema.GetMomentum(3,mom);
  mom_kp_x = mom[0];
  mom_kp_y = mom[1];
  mom_kp_z = mom[2];

  //  G4cout<<"kp mom:"<<mom[1]<<G4endl;

  Energy_h = Hkinema.GetEnergy(4);
  // G4double momentum_h = Hkinema.GetMomentum(4);
  Hkinema.GetMomentum(4,mom);
  mom_h_x = mom[0];
  mom_h_y = mom[1];
  mom_h_z = mom[2];
  //  cout<<mom_h_z<<endl;
  if((acos(mom[2]/sqrt(pow(mom[0],2)+pow(mom[1],2)+pow(mom[2],2)))*180/3.141592654)>20.) goto up;

  // G4double Theta_h = Hkinema.GetTheta(4);
  // G4double Phi_h = Hkinema.GetPhi(4);

  ////check kinematics
  // G4double momk[4];
  // momk[0]=mom_kp_x;
  // momk[1]=mom_kp_y;
  // momk[2]=mom_kp_z;
  // momk[3]=sqrt(mom_kp_x*mom_kp_x+mom_kp_y*mom_kp_y+mom_kp_z*mom_kp_z+0.493677*0.493677);
  //  G4double momkpp=sqrt(pow(mom_kp_x,2)+pow(mom_kp_y,2)+pow(mom_kp_z,2));
  // G4double momcmk[4]={0.};
  //  G4double beta= pbeam/(0.938272013+Ebeam);

  //  G4double test=lorentz(momk,beta,momcmk);
  // G4double momcmkpp=sqrt(pow(momcmk[0],2)+pow(momcmk[1],2)+pow(momcmk[2],2));
  //cmk->Fill(acos(momcmk[2]/momcmkpp)/3.141592654*180);
  //  G4double cmphik=atan2(momk[1],momk[0])*180/3.141592654;
  //coscmk->Fill(momcmk[2]/momcmkpp);
  //  G4double misskp = miss1(pbeam,rmk,momk);


  //  G4double vtx = 0.*CLHEP::mm;
  //  G4double vty = 0.*CLHEP::mm;

  G4double vtx = 0.*CLHEP::mm;
  G4double vty = 0.*CLHEP::mm;
  G4double vtz= G4RandFlat::shoot(env_target_pos_z-m_target_size.z()/2,env_target_pos_z+m_target_size.z()/2)*CLHEP::mm;

  //Kaon -
  m_particle_gun->SetParticleDefinition(m_KaonMinus);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_kp_x,mom_kp_y,mom_kp_z));
  m_particle_gun->SetParticleEnergy((Energy_kp - m_KaonMinus->GetPDGMass()/CLHEP::GeV)*CLHEP::GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex(anEvent);

  //proton
  m_particle_gun->SetParticleDefinition(m_Proton);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_h_x,mom_h_y,mom_h_z));
  m_particle_gun->SetParticleEnergy((Energy_h - m_Proton->GetPDGMass()/CLHEP::GeV)*CLHEP::GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex(anEvent);

  // gAnaMan.SetPrimaryParticle(0,mom_kp_x,mom_kp_y,mom_kp_z,m_KaonMinus->GetPDGMass()/CLHEP::GeV);
  // gAnaMan.SetPrimaryParticle(1,mom_h_x,mom_h_y,mom_h_z,m_Proton->GetPDGMass()/CLHEP::GeV);
  // gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  // gAnaMan.SetPrimaryVertex(1,vtx,vty,vtz);
}

//_____________________________________________________________________________
void
PrimaryGeneratorAction::GenerateE07StudyKnPBeam(G4Event* anEvent)
{
  G4double mom[3];
  G4double Ebeam, pbeam;
  //  G4double rmk=0.493677;
  G4double pbm[4];
  G4double Energy_h, mom_h_x, mom_h_y, mom_h_z;
  G4double Energy_kp, mom_kp_x, mom_kp_y, mom_kp_z;

  //  G4double pimom=0.635+G4RandFlat::shoot()*(2.000-0.635);
  pbeam=G4RandGauss::shoot(m_beam_p0,m_beam_p0*3.3*0.0001/2.3548);
  //  pbeam=G4RandGauss::shoot(m_beam_p0,m_beam_p0*3.3*0.0001/2.3548);
  // gAnaMan.SetPrimaryBeam(0,0,pbeam);
  Ebeam = sqrt(pbeam*pbeam+m_KaonMinus->GetPDGMass()/CLHEP::GeV*m_KaonMinus->GetPDGMass()/CLHEP::GeV);

  pbm[0]=0;
  pbm[1]=0;
  pbm[2]=pbeam;
  pbm[3]=Ebeam;
  KinemaHweak Hkinema(m_KaonMinus->GetPDGMass()/CLHEP::GeV,
                      m_Proton->GetPDGMass()/CLHEP::GeV,
                      m_KaonMinus->GetPDGMass()/CLHEP::GeV,
                      m_Proton->GetPDGMass()/CLHEP::GeV,
                      pbm[2], 0.0);

  Energy_kp=Hkinema.GetEnergy(3);
  // G4double momentum_kp = Hkinema.GetMomentum(3);
  Hkinema.GetMomentum(3,mom);
  mom_kp_x = mom[0];
  mom_kp_y = mom[1];
  mom_kp_z = mom[2];

  //  G4cout<<"kp mom:"<<mom[1]<<G4endl;

  Energy_h = Hkinema.GetEnergy(4);
  // G4double momentum_h = Hkinema.GetMomentum(4);
  Hkinema.GetMomentum(4,mom);
  mom_h_x = mom[0];
  mom_h_y = mom[1];
  mom_h_z = mom[2];
  //  cout<<mom_h_z<<endl;
  //    if((acos(mom[2]/sqrt(pow(mom[0],2)+pow(mom[1],2)+pow(mom[2],2)))*180/3.141592654)>20.) goto up;

  // G4double Theta_h = Hkinema.GetTheta(4);
  // G4double Phi_h = Hkinema.GetPhi(4);

  ////check kinematics
  // G4double momk[4];
  // momk[0]=mom_kp_x;
  // momk[1]=mom_kp_y;
  // momk[2]=mom_kp_z;
  // momk[3]=sqrt(mom_kp_x*mom_kp_x+mom_kp_y*mom_kp_y+mom_kp_z*mom_kp_z+0.493677*0.493677);
  //  G4double momkpp=sqrt(pow(mom_kp_x,2)+pow(mom_kp_y,2)+pow(mom_kp_z,2));
  // G4double momcmk[4]={0.};
  //  G4double beta= pbeam/(0.938272013+Ebeam);

  //  G4double test=lorentz(momk,beta,momcmk);
  // G4double momcmkpp=sqrt(pow(momcmk[0],2)+pow(momcmk[1],2)+pow(momcmk[2],2));
  //cmk->Fill(acos(momcmk[2]/momcmkpp)/3.141592654*180);
  //  G4double cmphik=atan2(momk[1],momk[0])*180/3.141592654;
  //coscmk->Fill(momcmk[2]/momcmkpp);
  //  G4double misskp = miss1(pbeam,rmk,momk);


  //  G4double vtx = 0.*CLHEP::mm;
  //  G4double vty = 0.*CLHEP::mm;

  G4double vtx = G4RandFlat::shoot(-15.,15.)*CLHEP::mm;
  G4double vty = G4RandFlat::shoot(-5.,5.)*CLHEP::mm;
  G4double vtz= G4RandFlat::shoot(env_target_pos_z-m_target_size.z()/2,env_target_pos_z+m_target_size.z()/2)*CLHEP::mm;

  //Kaon -
  m_particle_gun->SetParticleDefinition(m_KaonMinus);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_kp_x,mom_kp_y,mom_kp_z));
  m_particle_gun->SetParticleEnergy((Energy_kp - m_KaonMinus->GetPDGMass()/CLHEP::GeV)*CLHEP::GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex(anEvent);

  //proton
  m_particle_gun->SetParticleDefinition(m_Proton);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_h_x,mom_h_y,mom_h_z));
  m_particle_gun->SetParticleEnergy((Energy_h - m_Proton->GetPDGMass()/CLHEP::GeV)*CLHEP::GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex(anEvent);

  // gAnaMan.SetPrimaryParticle(0,mom_kp_x,mom_kp_y,mom_kp_z,m_KaonMinus->GetPDGMass()/CLHEP::GeV);
  // gAnaMan.SetPrimaryParticle(1,mom_h_x,mom_h_y,mom_h_z,m_Proton->GetPDGMass()/CLHEP::GeV);
  // gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  // gAnaMan.SetPrimaryVertex(1,vtx,vty,vtz);
}

//_____________________________________________________________________________
void
PrimaryGeneratorAction::GenerateE07StudyKpXiBeam(G4Event* anEvent)
{
  G4double mom[4];
  G4double Ebeam, pbeam;
  //  G4double rmk=0.493677;
  G4double pbm[4];
  G4double Energy_h, mom_h_x, mom_h_y, mom_h_z;
  G4double Energy_kp, mom_kp_x, mom_kp_y, mom_kp_z;

  //  G4double pimom=0.635+G4RandFlat::shoot()*(2.000-0.635);
  ////FWHM 3.3 x 10^-4, FWHM ~ 2.3548 sigma
  pbeam=G4RandGauss::shoot(m_beam_p0,m_beam_p0*3.3*0.0001/2.3548);
  //  pbeam=G4RandGauss::shoot(1.7,1.7*3.3*0.0001/2.3548);
  //  G4cout<<pbeam<<G4endl;

  // gAnaMan.SetPrimaryBeam(0,0,pbeam);
  Ebeam = sqrt(pbeam*pbeam+m_KaonMinus->GetPDGMass()/CLHEP::GeV*m_KaonMinus->GetPDGMass()/CLHEP::GeV);
  pbm[0]=0;
  pbm[1]=0;
  pbm[2]=pbeam;
  pbm[3]=Ebeam;
 up:
  KinemaHweak Hkinema(m_KaonMinus->GetPDGMass()/CLHEP::GeV,
                      m_Proton->GetPDGMass()/CLHEP::GeV,
                      m_XiMinus->GetPDGMass()/CLHEP::GeV,
                      m_KaonPlus->GetPDGMass()/CLHEP::GeV,
                      pbm[2], 0.0);

  Energy_kp=Hkinema.GetEnergy(4);
  // G4double momentum_kp = Hkinema.GetMomentum(4);
  Hkinema.GetMomentum(4,mom);
  mom_kp_x = mom[0];
  mom_kp_y = mom[1];
  mom_kp_z = mom[2];
  if((acos(mom[2]/sqrt(pow(mom[0],2)+pow(mom[1],2)+pow(mom[2],2)))*180/3.141592654)>20.)
    goto up;
  //  G4cout<<"kp mom:"<<mom[1]<<G4endl;

  Energy_h = Hkinema.GetEnergy(3);
  // G4double momentum_h = Hkinema.GetMomentum(3);
  Hkinema.GetMomentum(3,mom);
  mom_h_x = mom[0];
  mom_h_y = mom[1];
  mom_h_z = mom[2];
  //  cout<<mom_h_z<<endl;

  // G4double Theta_h = Hkinema.GetTheta(3);
  // G4double Phi_h = Hkinema.GetPhi(3);

  ////check kinematics
  // G4double momk[4];
  // momk[0]=mom_kp_x;
  // momk[1]=mom_kp_y;
  // momk[2]=mom_kp_z;
  // momk[3]=sqrt(mom_kp_x*mom_kp_x+mom_kp_y*mom_kp_y+mom_kp_z*mom_kp_z+0.493677*0.493677);
  //  G4double momkpp=sqrt(pow(mom_kp_x,2)+pow(mom_kp_y,2)+pow(mom_kp_z,2));
  // G4double momcmk[4]={0.};
  //  G4double beta= pbeam/(m_Proton->GetPDGMass()/CLHEP::GeV+Ebeam);

  //  G4double test=lorentz(momk,beta,momcmk);
  // G4double momcmkpp=sqrt(pow(momcmk[0],2)+pow(momcmk[1],2)+pow(momcmk[2],2));
  // cmk->Fill(acos(momcmk[2]/momcmkpp)/3.141592654*180);
  //  G4double cmphik=atan2(momk[1],momk[0])*180/3.141592654;
  //coscmk->Fill(momcmk[2]/momcmkpp);
  //  G4double misskp = miss1(pbeam,rmk,momk);


  //  G4double vtx = 0.*CLHEP::mm;
  //  G4double vty = 0.*CLHEP::mm;

  G4double vtx = G4RandFlat::shoot(-15.,15.)*CLHEP::mm;
  G4double vty = G4RandFlat::shoot(-5.,5.)*CLHEP::mm;
  G4double vtz= G4RandFlat::shoot(env_target_pos_z-m_target_size.z()/2,env_target_pos_z+m_target_size.z()/2)*CLHEP::mm;

  //Kaon -
  m_particle_gun->SetParticleDefinition(m_KaonPlus);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_kp_x,mom_kp_y,mom_kp_z));
  m_particle_gun->SetParticleEnergy((Energy_kp - m_KaonPlus->GetPDGMass()/CLHEP::GeV)*CLHEP::GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex(anEvent);

  //Proton
  m_particle_gun->SetParticleDefinition(m_XiMinus);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_h_x,mom_h_y,mom_h_z));
  m_particle_gun->SetParticleEnergy((Energy_h - m_XiMinus->GetPDGMass()/CLHEP::GeV)*CLHEP::GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex(anEvent);

  // gAnaMan.SetPrimaryParticle(0,mom_kp_x,mom_kp_y,mom_kp_z,m_KaonMinus->GetPDGMass()/CLHEP::GeV);
  // gAnaMan.SetPrimaryParticle(1,mom_h_x,mom_h_y,mom_h_z,m_XiMinus->GetPDGMass()/CLHEP::GeV);
  // gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  // gAnaMan.SetPrimaryVertex(1,vtx,vty,vtz);
}

//_____________________________________________________________________________
void
PrimaryGeneratorAction::GenerateE07StudyKpXiBeamOnlyKp(G4Event* anEvent)
{
  G4double mom[4];
  G4double Ebeam, pbeam;
  //  G4double rmk=0.493677;
  G4double pbm[4];
  G4double Energy_kp, mom_kp_x, mom_kp_y, mom_kp_z;

  //  G4double pimom=0.635+G4RandFlat::shoot()*(2.000-0.635);
  //  pbeam=1.7;
  //  pbeam=G4RandGauss::shoot(1.7,1.7*3.3*0.0001/2.3548);
  pbeam=G4RandGauss::shoot(m_beam_p0,m_beam_p0*3.3*0.0001/2.3548);
  // gAnaMan.SetPrimaryBeam(0,0,pbeam);
  Ebeam = sqrt(pbeam*pbeam+m_KaonMinus->GetPDGMass()/CLHEP::GeV*m_KaonMinus->GetPDGMass()/CLHEP::GeV);
  pbm[0]=0;
  pbm[1]=0;
  pbm[2]=pbeam;
  pbm[3]=Ebeam;
 up:
  KinemaHweak Hkinema(m_KaonMinus->GetPDGMass()/CLHEP::GeV,
                      m_Proton->GetPDGMass()/CLHEP::GeV,
                      m_XiMinus->GetPDGMass()/CLHEP::GeV,
                      m_KaonPlus->GetPDGMass()/CLHEP::GeV,
                      pbm[2], 0.0);

  Energy_kp=Hkinema.GetEnergy(4);
  // G4double momentum_kp = Hkinema.GetMomentum(4);
  Hkinema.GetMomentum(4,mom);
  mom_kp_x = mom[0];
  mom_kp_y = mom[1];
  mom_kp_z = mom[2];
  if((acos(mom[2]/sqrt(pow(mom[0],2)+pow(mom[1],2)+pow(mom[2],2)))*180/3.141592654)>20.)
    goto up;
  //  G4cout<<"kp mom:"<<mom[1]<<G4endl;

  // G4double Energy_h = Hkinema.GetEnergy(3);
  // G4double momentum_h = Hkinema.GetMomentum(3);
  Hkinema.GetMomentum(3,mom);
  // G4double mom_h_x = mom[0];
  // G4double mom_h_y = mom[1];
  // G4double mom_h_z = mom[2];
  //  cout<<mom_h_z<<endl;

  // G4double Theta_h = Hkinema.GetTheta(3);
  // G4double Phi_h = Hkinema.GetPhi(3);

  ////check kinematics
  // G4double momk[4];
  // momk[0]=mom_kp_x;
  // momk[1]=mom_kp_y;
  // momk[2]=mom_kp_z;
  // momk[3]=sqrt(mom_kp_x*mom_kp_x+mom_kp_y*mom_kp_y+mom_kp_z*mom_kp_z+0.493677*0.493677);
  //  G4double momkpp=sqrt(pow(mom_kp_x,2)+pow(mom_kp_y,2)+pow(mom_kp_z,2));
  // G4double momcmk[4]={0.};
  //  G4double beta= pbeam/(m_Proton->GetPDGMass()/CLHEP::GeV+Ebeam);

  //  G4double test=lorentz(momk,beta,momcmk);
  // G4double momcmkpp=sqrt(pow(momcmk[0],2)+pow(momcmk[1],2)+pow(momcmk[2],2));
  //cmk->Fill(acos(momcmk[2]/momcmkpp)/3.141592654*180);
  //  G4double cmphik=atan2(momk[1],momk[0])*180/3.141592654;
  //coscmk->Fill(momcmk[2]/momcmkpp);
  //  G4double misskp = miss1(pbeam,rmk,momk);


  //  G4double vtx = 0.*CLHEP::mm;
  //  G4double vty = 0.*CLHEP::mm;

  G4double vtx = G4RandFlat::shoot(-15.,15.)*CLHEP::mm;
  G4double vty = G4RandFlat::shoot(-5.,5.)*CLHEP::mm;
  G4double vtz= G4RandFlat::shoot(env_target_pos_z-m_target_size.z()/2,env_target_pos_z+m_target_size.z()/2)*CLHEP::mm;

  //Kaon -
  m_particle_gun->SetParticleDefinition(m_KaonPlus);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_kp_x,mom_kp_y,mom_kp_z));
  m_particle_gun->SetParticleEnergy((Energy_kp - m_KaonPlus->GetPDGMass()/CLHEP::GeV)*CLHEP::GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex(anEvent);

  //Proton
  /*  m_particle_gun->SetParticleDefinition(m_XiMinus);
      m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_h_x,mom_h_y,mom_h_z));
      m_particle_gun->SetParticleEnergy((Energy_h - m_XiMinus->GetPDGMass()/CLHEP::GeV)*CLHEP::GeV);
      m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
      m_particle_gun->GeneratePrimaryVertex(anEvent);
  */

  // gAnaMan.SetPrimaryParticle(0,mom_kp_x,mom_kp_y,mom_kp_z,m_KaonMinus->GetPDGMass()/CLHEP::GeV);
  // gAnaMan.SetPrimaryParticle(1,mom_h_x,mom_h_y,mom_h_z,m_XiMinus->GetPDGMass()/CLHEP::GeV);
  // gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  // gAnaMan.SetPrimaryVertex(1,vtx,vty,vtz);
}

//_____________________________________________________________________________
void
PrimaryGeneratorAction::GenerateE07StudyKpxi1530(G4Event* anEvent)
{
  G4double mom[4];
  G4double Ebeam, pbeam;
  //  G4double rmk=0.493677;
  G4double pbm[4];
  G4double Energy_h, mom_h_x, mom_h_y, mom_h_z;
  G4double Energy_kp, mom_kp_x, mom_kp_y, mom_kp_z;
  //  G4cout<<"test21"<<G4endl;
  //  G4cout<<m_Xi1530Minus->GetPDGMass()/CLHEP::GeV<<G4endl;
  //  G4double pimom=0.635+G4RandFlat::shoot()*(2.000-0.635);
  ////FWHM 3.3 x 10^-4, FWHM ~ 2.3548 sigma
  //  pbeam=G4RandGauss::shoot(1.7,1.7*3.3*0.0001/2.3548);
  pbeam=G4RandGauss::shoot(m_beam_p0,m_beam_p0*3.3*0.0001/2.3548);
  //  G4cout<<pbeam<<G4endl;

  // gAnaMan.SetPrimaryBeam(0,0,pbeam);
  Ebeam = sqrt(pbeam*pbeam+m_KaonMinus->GetPDGMass()/CLHEP::GeV*m_KaonMinus->GetPDGMass()/CLHEP::GeV);
  pbm[0]=0;
  pbm[1]=0;
  pbm[2]=pbeam;
  pbm[3]=Ebeam;
 up:
  KinemaHweak Hkinema(m_KaonMinus->GetPDGMass()/CLHEP::GeV,
                      m_Proton->GetPDGMass()/CLHEP::GeV,
                      m_Xi1530Minus->GetPDGMass()/CLHEP::GeV,
                      m_KaonPlus->GetPDGMass()/CLHEP::GeV,
                      pbm[2], 0.0);

  Energy_kp=Hkinema.GetEnergy(4);
  // G4double momentum_kp = Hkinema.GetMomentum(4);
  Hkinema.GetMomentum(4,mom);
  mom_kp_x = mom[0];
  mom_kp_y = mom[1];
  mom_kp_z = mom[2];

  if((atan2(mom_kp_y,mom_kp_z))/3.141592*180.<-18. ||
     (atan2(mom_kp_y,mom_kp_z))/3.141592*180.>20.)
    goto up;

  //    if((acos(mom[2]/sqrt(pow(mom[0],2)+pow(mom[1],2)+pow(mom[2],2)))*180/3.141592654)>20.) goto up;
  //  G4cout<<"kp mom:"<<mom[1]<<G4endl;

  Energy_h = Hkinema.GetEnergy(3);
  // G4double momentum_h = Hkinema.GetMomentum(3);
  Hkinema.GetMomentum(3,mom);
  mom_h_x = mom[0];
  mom_h_y = mom[1];
  mom_h_z = mom[2];
  //  cout<<mom_h_z<<endl;

  // G4double Theta_h = Hkinema.GetTheta(3);
  // G4double Phi_h = Hkinema.GetPhi(3);
  ////check kinematics
  // G4double momk[4];
  // momk[0]=mom_kp_x;
  // momk[1]=mom_kp_y;
  // momk[2]=mom_kp_z;
  // momk[3]=sqrt(mom_kp_x*mom_kp_x+mom_kp_y*mom_kp_y+mom_kp_z*mom_kp_z+0.493677*0.493677);
  //  G4double momkpp=sqrt(pow(mom_kp_x,2)+pow(mom_kp_y,2)+pow(mom_kp_z,2));
  // G4double momcmk[4]={0.};
  //  G4double beta= pbeam/(m_Proton->GetPDGMass()/CLHEP::GeV+Ebeam);

  //  G4double test=lorentz(momk,beta,momcmk);
  // G4double momcmkpp=sqrt(pow(momcmk[0],2)+pow(momcmk[1],2)+pow(momcmk[2],2));
  //cmk->Fill(acos(momcmk[2]/momcmkpp)/3.141592654*180.);
  //  G4double cmphik=atan2(momk[1],momk[0])*180/3.141592654;
  //coscmk->Fill(momcmk[2]/momcmkpp);
  //  G4double misskp = miss1(pbeam,rmk,momk);


  //  G4double vtx = 0.*CLHEP::mm;
  //  G4double vty = 0.*CLHEP::mm;

  G4double vtx = G4RandFlat::shoot(-15.,15.)*CLHEP::mm;
  G4double vty = G4RandFlat::shoot(-5.,5.)*CLHEP::mm;
  G4double vtz= G4RandFlat::shoot(env_target_pos_z-m_target_size.z()/2,env_target_pos_z+m_target_size.z()/2)*CLHEP::mm;

  //Kaon -
  m_particle_gun->SetParticleDefinition(m_KaonPlus);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_kp_x,mom_kp_y,mom_kp_z));
  m_particle_gun->SetParticleEnergy((Energy_kp - m_KaonPlus->GetPDGMass()/CLHEP::GeV)*CLHEP::GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex(anEvent);

  //  //Xi1530-
  m_particle_gun->SetParticleDefinition(m_Xi1530Minus);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_h_x,mom_h_y,mom_h_z));
  m_particle_gun->SetParticleEnergy((Energy_h - m_Xi1530Minus->GetPDGMass()/CLHEP::GeV)*CLHEP::GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex(anEvent);

  // gAnaMan.SetPrimaryParticle(0,mom_kp_x,mom_kp_y,mom_kp_z,m_KaonMinus->GetPDGMass()/CLHEP::GeV);
  // gAnaMan.SetPrimaryParticle(1,mom_h_x,mom_h_y,mom_h_z,m_Xi1530Minus->GetPDGMass()/CLHEP::GeV);
  // gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  // gAnaMan.SetPrimaryVertex(1,vtx,vty,vtz);
}

//_____________________________________________________________________________
void
PrimaryGeneratorAction::GenerateE07StudyTakahashi(G4Event* anEvent)
{
  G4double mom[3];
  //  G4double rmk=0.493677;
  G4double Energy_h, mom_h_x, mom_h_y, mom_h_z;
  G4double Energy_kp, mom_kp_x, mom_kp_y, mom_kp_z;

  G4int type = G4RandFlat::shoot() * 6. < 2. ? 0.:1.;
  //  G4cout<<"test111:"<<angular_mom<<G4endl;
  //  G4cout<<G4RandFlat::shoot()<<G4endl;
  auto p_proton = Kinematics::HarmonicFermiMomentum(type);
  G4double p_fermi = sqrt(pow(p_proton[0],2)+pow(p_proton[1],2)+
			  pow(p_proton[2],2));
  G4double ke = sqrt(p_fermi * p_fermi +
		     m_Proton->GetPDGMass()/CLHEP::GeV * m_Proton->GetPDGMass()/CLHEP::GeV);
  G4cout << ke << G4endl;
  G4double Al[13] = { 1., -1.22, 1.55, -1.08, 0.37,
		      -0.15, 0.16, -0.38, -0.18, 0.09,
		      0.05, -0.01, 0.20 };
  // from Dauber, et al., PR 179 5 1968 at K- beam energy 1.7 GeV/c

  G4double cross_section=0.;

  G4cout<<"p_proton:"<<p_proton[0]<<":"<<p_proton[1]<<":"<<p_proton[2]<<G4endl;
  //  G4double pimom=0.635+G4RandFlat::shoot()*(2.000-0.635);
  ////FWHM 3.3 x 10^-4, FWHM ~ 2.3548 sigma
  G4double pbeam=G4RandGauss::shoot(m_beam_p0,m_beam_p0*3.3*0.0001/2.3548);
  //  pbeam=G4RandGauss::shoot(1.7,1.7*3.3*0.0001/2.3548);
  //  G4cout<<pbeam<<G4endl;

  // gAnaMan.SetPrimaryBeam(0,0,pbeam);
  // G4double Ebeam = sqrt(pbeam*pbeam+
  // 	m_KaonMinus->GetPDGMass()/CLHEP::GeV*m_KaonMinus->GetPDGMass()/CLHEP::GeV);
  G4ThreeVector pbm(0., 0., pbeam);
  // up:
  G4double cosx=G4RandFlat::shoot(-1.,1.);
  KinemaFermi Hkinema(m_KaonMinus->GetPDGMass()/CLHEP::GeV,
		      ke,
		      m_XiMinus->GetPDGMass()/CLHEP::GeV,
		      m_KaonPlus->GetPDGMass()/CLHEP::GeV,
		      pbm, p_proton,cosx);
  for(G4int le=0;le<13;le++){
    cross_section=cross_section+Al[le];
  }

  Energy_kp=Hkinema.GetEnergy(4);
  // G4double momentum_kp = Hkinema.GetMomentum(4);
  Hkinema.GetMomentum(4,mom);
  mom_kp_x = mom[0];
  mom_kp_y = mom[1];
  mom_kp_z = mom[2];
  //    if((acos(mom[2]/sqrt(pow(mom[0],2)+pow(mom[1],2)+pow(mom[2],2)))*180/3.141592654)>20.) goto up;
  //  G4cout<<"kp mom:"<<mom[1]<<G4endl;

  Energy_h = Hkinema.GetEnergy(3);
  // G4double momentum_h = Hkinema.GetMomentum(3);
  Hkinema.GetMomentum(3,mom);
  mom_h_x = mom[0];
  mom_h_y = mom[1];
  mom_h_z = mom[2];

  //  cout<<mom_h_z<<endl;

  // G4double Theta_h = Hkinema.GetTheta(3);
  // G4double Phi_h = Hkinema.GetPhi(3);

  ////check kinematics
  // G4double momk[4];
  // momk[0]=mom_kp_x;
  // momk[1]=mom_kp_y;
  // momk[2]=mom_kp_z;
  // momk[3]=sqrt(mom_kp_x*mom_kp_x+mom_kp_y*mom_kp_y+mom_kp_z*mom_kp_z+0.493677*0.493677);
  //  G4double momkpp=sqrt(pow(mom_kp_x,2)+pow(mom_kp_y,2)+pow(mom_kp_z,2));
  // G4double momcmk[4]={0.};
  //  G4double beta= pbeam/(m_Proton->GetPDGMass()/CLHEP::GeV+Ebeam);

  //  G4double test=lorentz(momk,beta,momcmk);
  // G4double momcmkpp=sqrt(pow(momcmk[0],2)+pow(momcmk[1],2)+pow(momcmk[2],2));
  //cmk->Fill(acos(momcmk[2]/momcmkpp)/3.141592654*180);
  //  G4double cmphik=atan2(momk[1],momk[0])*180/3.141592654;
  //coscmk->Fill(momcmk[2]/momcmkpp);
  //  G4double misskp = miss1(pbeam,rmk,momk);


  //  G4double vtx = 0.*CLHEP::mm;
  //  G4double vty = 0.*CLHEP::mm;
  G4double vtx = G4RandFlat::shoot(-15.,15.)*CLHEP::mm;
  G4double vty = G4RandFlat::shoot(-5.,5.)*CLHEP::mm;
  G4double vtz= G4RandFlat::shoot(env_target_pos_z-m_target_size.z()/2,env_target_pos_z+m_target_size.z()/2)*CLHEP::mm;

  //Kaon -
  m_particle_gun->SetParticleDefinition(m_KaonPlus);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_kp_x,mom_kp_y,mom_kp_z));
  m_particle_gun->SetParticleEnergy((Energy_kp - m_KaonPlus->GetPDGMass()/CLHEP::GeV)*CLHEP::GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex(anEvent);

  //Proton
  m_particle_gun->SetParticleDefinition(m_XiMinus);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_h_x,mom_h_y,mom_h_z));
  m_particle_gun->SetParticleEnergy((Energy_h - m_XiMinus->GetPDGMass()/CLHEP::GeV)*CLHEP::GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex(anEvent);

  // gAnaMan.SetPrimaryParticle(0,mom_kp_x,mom_kp_y,mom_kp_z,m_KaonMinus->GetPDGMass()/CLHEP::GeV);
  // gAnaMan.SetPrimaryParticle(1,mom_h_x,mom_h_y,mom_h_z,m_XiMinus->GetPDGMass()/CLHEP::GeV);
  // gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  // gAnaMan.SetPrimaryVertex(1,vtx,vty,vtz);
}

//_____________________________________________________________________________
void
PrimaryGeneratorAction::GenerateE07StudyP08to20(G4Event* anEvent)
{
  G4double Energy_p, mom_p_x, mom_p_y, mom_p_z;
  //  G4double mom[3];
  //  G4double Ebeam,pg_x,pg_y,pg_z, pbeam;

  /* proton */
 up:
  G4double phi=G4RandFlat::shoot(-1.,1.)*3.141592654;
  G4double mom_p=G4RandFlat::shoot(0.8,2.0);
  //  G4double theta=acos(G4RandFlat::shoot(-1.,1.));
  ///  G4double theta=acos(G4RandFlat::shoot(0.939692646,1.));
  //  G4double theta=acos(G4RandFlat::shoot(0.866025458,1.));//30 deg
  G4double theta=acos(G4RandFlat::shoot(0.80,1.));//

  //  G4double theta=acos(G4RandFlat::shoot(0.962646,0.97));
  //  G4cout<<theta*180/3.141592<<G4endl;
  mom_p_x = mom_p*sin(theta)*cos(phi);
  mom_p_y = mom_p*sin(theta)*sin(phi);
  mom_p_z = mom_p*cos(theta);

  if((atan2(mom_p_y,mom_p_z))/3.141592*180.<-18. || (atan2(mom_p_y,mom_p_z))/3.141592*180.>18.)
    goto up;

  //  G4cout<<atan2(mom_p_y,mom_p_z)*180./3.141592<<G4endl;

  // labk->Fill(theta*180/3.141592);
  // coslabk->Fill(cos(theta));
  // phik->Fill(phi/3.141592*180);

  //  G4cout<<mom_p_x<<":"<<mom_p_y<<":"<<mom_p_z<<G4endl;

  Energy_p=sqrt(pow(mom_p,2)+pow(m_Proton->GetPDGMass()/CLHEP::GeV,2));
  /////vertex
  G4double vtx = G4RandFlat::shoot(-15.,15.)*CLHEP::mm;
  G4double vty = G4RandFlat::shoot(-5.,5.)*CLHEP::mm;
  G4double vtz = G4RandFlat::shoot(env_target_pos_z-m_target_size.z()/2,env_target_pos_z+m_target_size.z()/2)*CLHEP::mm;

  ///////Proton PG
  m_particle_gun->SetParticleDefinition(m_Proton);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_p_x,mom_p_y,mom_p_z));
  m_particle_gun->SetParticleEnergy((Energy_p - m_Proton->GetPDGMass()/CLHEP::GeV)*CLHEP::GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex(anEvent);

  // gAnaMan.SetPrimaryParticle(0,mom_p_x,mom_p_y,mom_p_z,m_Proton->GetPDGMass()/CLHEP::GeV);
  // gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
}

//_____________________________________________________________________________
void
PrimaryGeneratorAction::GenerateE07StudyKp04to15(G4Event* anEvent)
{
  //  G4cout<<"e07_generator"<<G4endl;
  G4double Energy_p, mom_p_x, mom_p_y, mom_p_z;
  //  G4double mom[3];
  //  G4double Ebeam,pg_x,pg_y,pg_z, pbeam;

  /* kaon+ */
 up:
  G4double phi=G4RandFlat::shoot(-1.,1.)*3.141592654;
  G4double mom_p=G4RandFlat::shoot(0.4,1.5);
  //  G4double theta=acos(G4RandFlat::shoot(-1.,1.));
  //  G4double theta=acos(G4RandFlat::shoot(0.866025458,1.));//30 deg

  G4double theta=acos(G4RandFlat::shoot(0.8,1.));
  //  G4double theta=acos(G4RandFlat::shoot(0.939692646,1.));//25
  //  G4double theta=acos(G4RandFlat::shoot(0.962646,0.97));


  //  G4cout<<theta*180/3.141592<<G4endl;
  mom_p_x = mom_p*sin(theta)*cos(phi);
  mom_p_y = mom_p*sin(theta)*sin(phi);
  mom_p_z = mom_p*cos(theta);

  if((atan2(mom_p_y,mom_p_z))/3.141592*180.<-18. || (atan2(mom_p_y,mom_p_z))/3.141592*180.>18.)
    goto up;

  // labk->Fill(theta*180/3.141592);
  // coslabk->Fill(cos(theta));
  // phik->Fill(phi/3.141592*180);

  Energy_p=sqrt(pow(mom_p,2)+pow(m_KaonPlus->GetPDGMass()/CLHEP::GeV,2));
  //  Energy_p=pow(mom_p,2)+m_KaonPlus->GetPDGMass()/CLHEP::GeV;
  /////vertex
  G4double vtx = G4RandFlat::shoot(-15.,15.)*CLHEP::mm;
  G4double vty = G4RandFlat::shoot(-5.,5.)*CLHEP::mm;
  G4double vtz = G4RandFlat::shoot(env_target_pos_z-m_target_size.z()/2,
                                   env_target_pos_z+m_target_size.z()/2)*CLHEP::mm;

  ///////Proton PG
  m_particle_gun->SetParticleDefinition(m_KaonPlus);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_p_x,mom_p_y,mom_p_z));
  m_particle_gun->SetParticleEnergy((Energy_p - m_KaonPlus->GetPDGMass()/CLHEP::GeV)*CLHEP::GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex(anEvent);

  // gAnaMan.SetPrimaryParticle(0,mom_p_x,mom_p_y,mom_p_z,m_KaonPlus->GetPDGMass()/CLHEP::GeV);
  // gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
}

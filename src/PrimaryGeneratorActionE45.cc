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
// E45 elastic scattering pip
void
PrimaryGeneratorAction::GenerateE45ElasticPionPlus(G4Event* anEvent)
{
  G4double mom[3];
  G4double Ebeam, pbeam=0.;

  G4double pbm[4];
  G4double Energy_h, mom_h_x, mom_h_y, mom_h_z;
  G4double Energy_kp, mom_kp_x, mom_kp_y, mom_kp_z;
  G4double protonMass=m_Proton->GetPDGMass()/CLHEP::GeV;//unit GeV
  G4double pipMass=m_PionPlus->GetPDGMass()/CLHEP::GeV;//unit GeV

  /*
 ///first w/o beam
 ///beam optics from simulation file.
 char fname[100] ;
 FILE *fp;
 sprintf(fname,"./beam_simulation/profile_ve07-5.dat.txt");
 if ((fp = fopen(fname,"r")) == NULL){
 fprintf(stderr, ": Cannot open file: %s\n", fname);
 exit(-1);
 }

 double data[100]={-9999.9999};
 int check;
 int ran;
 int res;
 up1:
 check=0.;
 ran=G4RandFlat::shoot(1,17005);


 while(1){
 /// file structure : x, u, y, v, p, PID, ???
 /// unit           : cm, mrad, cm, mrad, gev, PID, ???
 res=fscanf(fp,"%lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf",&data[0],&data[1],&data[2],&data[3],&data[4],&data[5], &data[6]);
 check=check+1.;
 if(check==ran) break;
 //    if(res==EOF) break;
 }

 G4double dxdz,dydz,pp;
 dxdz=atan(data[1]*0.001);
 dydz=atan(data[3]*0.001);
 pp=data[4]/1.8*m_beam_p0;

 //  G4cout<<"pp:"<<pp<<G4endl;
 //  G4cout<<"data 4:"<<data[4]<<G4endl;


 */
  // G4double x0;
  // G4double y0;
  // G4double z0;
  // G4double mom_beam;
  G4String particleName;
  // G4ParticleDefinition* particle;

  G4double vtx;
  G4double vty;
  G4double vtz;
  ///for E45
  G4double rn_vtx,rn_vtz;
  while(1){
    rn_vtx = G4RandFlat::shoot(-m_target_size.x(), m_target_size.x());
    rn_vtz = G4RandFlat::shoot(-m_target_size.x(), m_target_size.x());
    //    G4cout<<rn_vtx<<G4endl;
    if((rn_vtx*rn_vtx+rn_vtz*rn_vtz) < m_target_size.x()*m_target_size.x()) break;
  }
  vty = G4RandFlat::shoot(-m_target_size.z(), m_target_size.z());
  vtx=rn_vtx;
  vtz=rn_vtz+env_target_pos_z;
  /*   ///E42
       vtz= G4RandFlat::shoot(env_Target_pos_z-env_Target_width/2.,env_Target_pos_z+env_Target_width/2.)*CLHEP::mm;
       vtx = (data[0]*10.+dxdz*(vtz-env_Target_pos_z))*CLHEP::mm;
       vty = (data[2]*10.+dydz*(vtz-env_Target_pos_z))*CLHEP::mm;

       //  G4cout<<"vtx and vty ::"<<fabs(vtx)<<", "<<fabs(vty)<<G4endl;
       //  G4cout<<"target_x, target_y::"<<env_Target_x<<", "<<env_Target_y<<G4endl;
       if(fabs(vtx)>env_Target_x/2. || fabs(vty)>env_Target_y/2.){
       goto up1;
       }
  */
  //  G4cout<<"passed-->vtx and vty ::"<<fabs(vtx)<<", "<<fabs(vty)<<G4endl;

  ///close file, if not it will make an error.
  //  fclose(fp);

  G4double pbeam_x;
  G4double pbeam_y;
  G4double pbeam_z;
  G4double pp;
  pp=m_beam_p0+G4RandGauss::shoot(0.,0.01294*m_beam_p0);
  pbeam_z=pp;
  pbeam_x=0.;
  pbeam_y=0.;

  // gAnaMan.SetPrimaryBeam(pbeam_x,pbeam_y,pbeam_z);
  G4double cosx = G4RandFlat::shoot(-1.,1.);
  G4double p_proton[4]={0};
  //  gAnaMan.SetFermiMotion(p_proton);

  ///prepare real distribution
  Ebeam = sqrt(pbeam*pbeam+pipMass/CLHEP::GeV*pipMass/CLHEP::GeV);
  pbm[0]=pbeam_x;
  pbm[1]=pbeam_y;
  pbm[2]=pbeam_z;
  pbm[3]=Ebeam;
  ///first m_PionPlus
  KinemaFermi Hkinema(m_PionPlus->GetPDGMass()/CLHEP::GeV,
                      m_Proton->GetPDGMass()/CLHEP::GeV,
                      m_PionPlus->GetPDGMass()/CLHEP::GeV,
                      m_Proton->GetPDGMass()/CLHEP::GeV,
                      pbm, p_proton,cosx);

  Energy_kp=Hkinema.GetEnergy(3);
  // G4double momentum_kp = Hkinema.GetMomentum(3);
  Hkinema.GetMomentum(3,mom);
  mom_kp_x = mom[0];
  mom_kp_y = mom[1];
  mom_kp_z = mom[2];
  //  if((acos(mom[2]/sqrt(pow(mom[0],2)+pow(mom[1],2)+pow(mom[2],2)))*180/3.141592654)>30.) goto up_sh;
  //  G4cout<<"kp mom:"<<mom[2]<<G4endl;

  //  cross_section->Fill(cosx);
  //  gAnaMan.SetCrossSection(cosx);

  Energy_h = Hkinema.GetEnergy(4);
  // G4double momentum_h = Hkinema.GetMomentum(4);
  Hkinema.GetMomentum(4,mom);
  mom_h_x = mom[0];
  mom_h_y = mom[1];
  mom_h_z = mom[2];

  //  cout<<mom_h_z<<endl;
  // G4double Theta_h = Hkinema.GetTheta(4);
  // G4double Phi_h = Hkinema.GetPhi(4);
  ////check kinematics
  // G4double momk[4];
  // momk[0]=mom_kp_x;
  // momk[1]=mom_kp_y;
  // momk[2]=mom_kp_z;
  // momk[3]=sqrt(mom_kp_x*mom_kp_x+mom_kp_y*mom_kp_y+mom_kp_z*mom_kp_z+pipMass*pipMass);
  //  G4double momkpp=sqrt(pow(mom_kp_x,2)+pow(mom_kp_y,2)+pow(mom_kp_z,2));
  // G4double momcmk[4]={0};
  //  G4double beta= pbeam/(m_Proton->GetPDGMass()/CLHEP::GeV+Ebeam);
  //  G4double test=lorentz(momk,beta,momcmk);
  // G4double momcmkpp=sqrt(pow(momcmk[0],2)+pow(momcmk[1],2)+pow(momcmk[2],2));
  //cmk->Fill(acos(momcmk[2]/momcmkpp)/3.141592654*180);
  //  G4double cmphik=atan2(momk[1],momk[0])*180/3.141592654;
  //coscmk->Fill(momcmk[2]/momcmkpp);
  //  G4double misskp = miss1(pbeam,rmk,momk);

  //pi
  m_particle_gun->SetParticleDefinition(m_PionPlus);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_kp_x,mom_kp_y,mom_kp_z));
  m_particle_gun->SetParticleEnergy((Energy_kp - pipMass/CLHEP::GeV)*CLHEP::GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex(anEvent);

  //proton
  m_particle_gun->SetParticleDefinition(m_Proton);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_h_x,mom_h_y,mom_h_z));
  m_particle_gun->SetParticleEnergy((Energy_h - protonMass/CLHEP::GeV)*CLHEP::GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex(anEvent);

  // gAnaMan.SetPrimaryParticle(0,mom_kp_x,mom_kp_y,mom_kp_z,pipMass/CLHEP::GeV);///pip
  // gAnaMan.SetPrimaryParticle(1,mom_h_x,mom_h_y,mom_h_z,protonMass/CLHEP::GeV);///proton
  // gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  // gAnaMan.SetPrimaryVertex(1,vtx,vty,vtz);
}

//_____________________________________________________________________________
// E45 elastic scattering pin
void
PrimaryGeneratorAction::GenerateE45ElasticPionMinus(G4Event* anEvent)
{
  G4double mom[3];
  G4double Ebeam, pbeam=0.;
  G4double pbm[4];
  G4double Energy_h, mom_h_x, mom_h_y, mom_h_z;
  G4double Energy_kp, mom_kp_x, mom_kp_y, mom_kp_z;
  G4double protonMass=m_Proton->GetPDGMass()/CLHEP::GeV;//unit CLHEP::GeV
  G4double pinMass=m_PionMinus->GetPDGMass()/CLHEP::GeV;//unit CLHEP::GeV

  /*
 ///first w/o beam
 ///beam optics from simulation file.
 char fname[100] ;
 FILE *fp;
 sprintf(fname,"./beam_simulation/profile_ve07-5.dat.txt");
 if ((fp = fopen(fname,"r")) == NULL){
 fprintf(stderr, ": Cannot open file: %s\n", fname);
 exit(-1);
 }

 double data[100]={-9999.9999};
 int check;
 int ran;
 int res;
 up1:
 check=0.;
 ran=G4RandFlat::shoot(1,17005);


 while(1){
 /// file structure : x, u, y, v, p, PID, ???
 /// unit           : cm, mrad, cm, mrad, gev, PID, ???
 res=fscanf(fp,"%lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf",&data[0],&data[1],&data[2],&data[3],&data[4],&data[5], &data[6]);
 check=check+1.;
 if(check==ran) break;
 //    if(res==EOF) break;
 }

 G4double dxdz,dydz,pp;
 dxdz=atan(data[1]*0.001);
 dydz=atan(data[3]*0.001);
 pp=data[4]/1.8*m_beam_p0;

 //  G4cout<<"pp:"<<pp<<G4endl;
 //  G4cout<<"data 4:"<<data[4]<<G4endl;


 */
  // G4double x0;
  // G4double y0;
  // G4double z0;
  // G4double mom_beam;
  G4String particleName;
  // G4ParticleDefinition* particle;

  G4double vtx;
  G4double vty;
  G4double vtz;
  ///for E45
  G4double rn_vtx,rn_vtz;
  while(1){
    rn_vtx = G4RandFlat::shoot(-m_target_size.x(), m_target_size.x());
    rn_vtz = G4RandFlat::shoot(-m_target_size.x(), m_target_size.x());
    //    G4cout<<rn_vtx<<G4endl;
    if((rn_vtx*rn_vtx+rn_vtz*rn_vtz) < m_target_size.x()*m_target_size.x()) break;
  }
  vty = G4RandFlat::shoot(-m_target_size.z(), m_target_size.z());
  vtx=rn_vtx;
  vtz=rn_vtz+env_target_pos_z;
  /*   ///E42
       vtz= G4RandFlat::shoot(env_Target_pos_z-env_Target_width/2.,env_Target_pos_z+env_Target_width/2.)*CLHEP::mm;
       vtx = (data[0]*10.+dxdz*(vtz-env_Target_pos_z))*CLHEP::mm;
       vty = (data[2]*10.+dydz*(vtz-env_Target_pos_z))*CLHEP::mm;

       //  G4cout<<"vtx and vty ::"<<fabs(vtx)<<", "<<fabs(vty)<<G4endl;
       //  G4cout<<"target_x, target_y::"<<env_Target_x<<", "<<env_Target_y<<G4endl;
       if(fabs(vtx)>env_Target_x/2. || fabs(vty)>env_Target_y/2.){
       goto up1;
       }
  */
  //  G4cout<<"passed-->vtx and vty ::"<<fabs(vtx)<<", "<<fabs(vty)<<G4endl;

  ///close file, if not it will make an error.
  //  fclose(fp);

  G4double pbeam_x;
  G4double pbeam_y;
  G4double pbeam_z;
  G4double pp;
  pp=m_beam_p0+G4RandGauss::shoot(0.,0.01294*m_beam_p0);
  pbeam_z=pp;
  pbeam_x=0.;
  pbeam_y=0.;

  // gAnaMan.SetPrimaryBeam(pbeam_x,pbeam_y,pbeam_z);
  G4double cosx = G4RandFlat::shoot(-1.,1.);
  G4double p_proton[4]={0};
  //  gAnaMan.SetFermiMotion(p_proton);

  ///prepare real distribution
  Ebeam = sqrt(pbeam*pbeam+pinMass/CLHEP::GeV*pinMass/CLHEP::GeV);
  pbm[0]=pbeam_x;
  pbm[1]=pbeam_y;
  pbm[2]=pbeam_z;
  pbm[3]=Ebeam;
  ///first pin
  KinemaFermi Hkinema(m_PionMinus->GetPDGMass()/CLHEP::GeV,
		      m_Proton->GetPDGMass()/CLHEP::GeV,
		      m_PionMinus->GetPDGMass()/CLHEP::GeV,
		      m_Proton->GetPDGMass()/CLHEP::GeV,
		      pbm, p_proton,cosx);

  Energy_kp=Hkinema.GetEnergy(3);
  // G4double momentum_kp = Hkinema.GetMomentum(3);
  Hkinema.GetMomentum(3,mom);
  mom_kp_x = mom[0];
  mom_kp_y = mom[1];
  mom_kp_z = mom[2];
  //  if((acos(mom[2]/sqrt(pow(mom[0],2)+pow(mom[1],2)+pow(mom[2],2)))*180/3.141592654)>30.) goto up_sh;
  //  G4cout<<"kp mom:"<<mom[2]<<G4endl;

  //  cross_section->Fill(cosx);
  //  gAnaMan.SetCrossSection(cosx);

  Energy_h = Hkinema.GetEnergy(4);
  // G4double   momentum_h = Hkinema.GetMomentum(4);
  Hkinema.GetMomentum(4,mom);
  mom_h_x = mom[0];
  mom_h_y = mom[1];
  mom_h_z = mom[2];

  //  cout<<mom_h_z<<endl;
  // G4double   Theta_h = Hkinema.GetTheta(4);
  // G4double   Phi_h = Hkinema.GetPhi(4);
  ////check kinematics
  // G4double momk[4];
  // momk[0]=mom_kp_x;
  // momk[1]=mom_kp_y;
  // momk[2]=mom_kp_z;
  // momk[3]=sqrt(mom_kp_x*mom_kp_x+mom_kp_y*mom_kp_y+mom_kp_z*mom_kp_z+pinMass*pinMass);
  //  G4double momkpp=sqrt(pow(mom_kp_x,2)+pow(mom_kp_y,2)+pow(mom_kp_z,2));
  // G4double momcmk[4]={0};
  //  G4double beta= pbeam/(m_Proton->GetPDGMass()/CLHEP::GeV+Ebeam);
  //  G4double test=lorentz(momk,beta,momcmk);
  // G4double momcmkpp=sqrt(pow(momcmk[0],2)+pow(momcmk[1],2)+pow(momcmk[2],2));
  //cmk->Fill(acos(momcmk[2]/momcmkpp)/3.141592654*180);
  //  G4double cmphik=atan2(momk[1],momk[0])*180/3.141592654;
  //coscmk->Fill(momcmk[2]/momcmkpp);
  //  G4double misskp = miss1(pbeam,rmk,momk);

  //pi
  m_particle_gun->SetParticleDefinition(m_PionMinus);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_kp_x,mom_kp_y,mom_kp_z));
  m_particle_gun->SetParticleEnergy((Energy_kp - pinMass/CLHEP::GeV)*CLHEP::GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex(anEvent);

  //proton
  m_particle_gun->SetParticleDefinition(m_Proton);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_h_x,mom_h_y,mom_h_z));
  m_particle_gun->SetParticleEnergy((Energy_h - protonMass/CLHEP::GeV)*CLHEP::GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex(anEvent);

  // gAnaMan.SetPrimaryParticle(0,mom_kp_x,mom_kp_y,mom_kp_z,pinMass/CLHEP::GeV);///pin
  // gAnaMan.SetPrimaryParticle(1,mom_h_x,mom_h_y,mom_h_z,protonMass/CLHEP::GeV);///proton
  // gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  // gAnaMan.SetPrimaryVertex(1,vtx,vty,vtz);
}

// -*- C++ -*-

#include "PrimaryGeneratorAction.hh"

#include <G4Event.hh>
#include <G4IonTable.hh>
#include <G4ParticleGun.hh>
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
PrimaryGeneratorAction::GenerateE42Hdibaryon2(G4Event* anEvent)
{
  G4double mass_hdibaryon, width_hdibaryon;
  //  G4double Energy_h, momentum_h, mom_h_x, mom_h_y, mom_h_z;

  G4double Energy_L1, mom_L1_x, mom_L1_y, mom_L1_z;
  G4double Energy_L2, mom_L2_x, mom_L2_y, mom_L2_z;
  G4double Energy_kp, mom_kp_x, mom_kp_y, mom_kp_z;
  G4double mom[3];
  G4double pg_x,pg_y,pg_z;
  // G4double rmk=0.493677;
  //  Ebeam = 1.9;       // Incident gamma energy (GeV)
  pg_x = 0.0;
  pg_y = 0.0;
  //  G4double pbeam=G4RandGauss::shoot(m_beam_p0,m_beam_p0*3.3*0.0001/2.3548)*CLHEP::GeV;
  //  pg_z = m_beam_p0;
  pg_z = G4RandGauss::shoot(m_beam_p0, 0.01294*m_beam_p0);
  //  G4cout<<"Ebeam:"<<Ebeam<<G4endl;
  // G4double pbeam=sqrt(pow(pg_x,2)+pow(pg_y,2)+pow(pg_z,2));
  // G4double Ebeam = sqrt(pbeam*pbeam+m_KaonMinus->GetPDGMass()/CLHEP::GeV*m_KaonMinus->GetPDGMass()/CLHEP::GeV);       // Incident gamma energy (GeV)

  // gAnaMan.SetPrimaryBeam(pg_x,pg_y,pg_z);
  mass_hdibaryon = gConf.Get<G4double>("HdibaryonMass");
  width_hdibaryon = gConf.Get<G4double>("HdibaryonWidth");

 up:
  Kinema3Resonance Hkinema(m_KaonMinus->GetPDGMass()/CLHEP::GeV,
                           ((m_Proton->GetPDGMass()/CLHEP::GeV)*2),
                           m_Lambda->GetPDGMass()/CLHEP::GeV,
                           m_Lambda->GetPDGMass()/CLHEP::GeV,
                           m_KaonPlus->GetPDGMass()/CLHEP::GeV,
                           mass_hdibaryon, width_hdibaryon, pg_z, 0.0);

  Energy_kp = Hkinema.GetEnergy(5);
  // G4double momentum_kp = Hkinema.GetMomentum(5);
  Hkinema.GetMomentum(5,mom);

  //  if(atan((mom[1])/mom[2])*180/3.141592654 > 15.) goto up;
  if(fabs(atan2(mom[1],mom[2])*180/3.141592654) > 15. || fabs(atan2(mom[0],mom[2])*180/3.141592654) > 20.) goto up;
  //  if((acos(mom[2]/sqrt(pow(mom[0],2)+pow(mom[1],2)+pow(mom[2],2)))*180/3.141592654)>15.) goto up;

  /* L1 */
  Energy_L1 = Hkinema.GetEnergy(3);
  // G4double momentum_L1 = Hkinema.GetMomentum(3);
  Hkinema.GetMomentum(3,mom);

  mom_L1_x = mom[0];
  mom_L1_y = mom[1];
  mom_L1_z = mom[2];

  // G4double ThetaL1 = Hkinema.GetTheta(3);
  // G4double PhiL1 = Hkinema.GetPhi(3);

  /* L2 */
  Energy_L2 = Hkinema.GetEnergy(4);
  // G4double momentum_L2 = Hkinema.GetMomentum(4);
  Hkinema.GetMomentum(4,mom);

  mom_L2_x = mom[0];
  mom_L2_y = mom[1];
  mom_L2_z = mom[2];

  // G4double ThetaL2 = Hkinema.GetTheta(4);
  // G4double PhiL2 = Hkinema.GetPhi(4);

  /* Kp */
  Energy_kp = Hkinema.GetEnergy(5);
  // G4double momentum_kp = Hkinema.GetMomentum(5);
  Hkinema.GetMomentum(5,mom);

  mom_kp_x = mom[0];
  mom_kp_y = mom[1];
  mom_kp_z = mom[2];

  //  Thetakp = Hkinema.GetThetaCM(5);
  // G4double Thetakp = Hkinema.GetThetaCM(1);
  //  Phikp = Hkinema.GetPhi(5);
  //  G4cout<<Phikp<<G4endl;
  // G4double shphi=(atan2(mom_kp_y,mom_kp_x))/3.141592654*180.;
  //phik->Fill(shphi);

  //  G4cout<<"phi test: "<<Phikp<<":"<<shphi<<G4endl;
  //  G4cout<<"test: "<<Hkinema.GetPhiCM(1)<<":"<<Hkinema.GetPhi(5)<<G4endl;
  // Vertex (Reaction) point
  //  G4cout<<Energy_kp<<G4endl;
  //  //  G4double beta= sqrt(1.8*1.8-0.493677*0.493677)/(2*0.93827203+Energy_kp);
  //  G4cout<<(m_Proton->GetPDGMass()/CLHEP::GeV)<<G4endl;
  // G4double beta= pbeam/(2.*0.938272013+Ebeam);
  //  G4cout<<"Pri:"<<beta<<G4endl;
  //  G4double beta= sqrt(pow(mom_kp_x,2)+pow(mom_kp_y,2)+pow(mom_kp_z,2))/(2*0.93827203+1.8);
  // G4double test;
  // G4double momk[4]={0.};
  // G4double momkpp;
  // G4double momcmk[4]={0.};
  // G4double momcmkpp;

  // momk[0]=mom_kp_x;
  // momk[1]=mom_kp_y;
  // momk[2]=mom_kp_z;
  // momk[3]=sqrt(mom_kp_x*mom_kp_x+mom_kp_y*mom_kp_y+mom_kp_z*mom_kp_z+0.493677*0.493677);
  // momkpp=sqrt(pow(mom_kp_x,2)+pow(mom_kp_y,2)+pow(mom_kp_z,2));
  //  G4cout<<"beta:"<<beta<<G4endl;
  // momcmkpp=sqrt(pow(momcmk[0],2)+pow(momcmk[1],2)+pow(momcmk[2],2));
  //  G4cout<<"mom comp-----------------------"<<G4endl;
  //  G4cout<<"mom comp:"<<momk[2]<<":"<<momcmk[2]<<G4endl;
  //  G4cout<<"mom comp:"<<momk[3]<<":"<<momcmk[3]<<G4endl;
  //cmk->Fill(acos(momcmk[2]/momcmkpp)/3.141592654*180);

  //  cmk->Fill(Thetakp);
  //  G4cout<<"Theta LAB :"<<":"<<acos(momk[2]/momkpp)*180/3.141592654<<G4endl;
  //  G4cout<<"Theta CM gen:"<<Thetakp<<G4endl;
  //  G4cout<<"Theta CM cal:"<<acos(momcmk[2]/momcmkpp)*180/3.141592654<<G4endl;
  // G4double cmphik=atan2(momk[1],momk[0])*180/3.141592654;

  //  G4cout<<"sh phi:"<<shphi<<G4endl;
  //  G4cout<<"cm phi:"<<cmphik<<G4endl;

  //coscmk->Fill((momcmk[2]/momcmkpp));

  // labk->Fill(acos(momk[2]/momkpp)*180/3.141592654);
  // coslabk->Fill(momk[2]/momkpp);


  ////check angle in the CM frame
  G4double hlab[4];
  hlab[0]=mom_L1_x+mom_L2_x;
  hlab[1]=mom_L1_y+mom_L2_y;
  hlab[2]=mom_L1_z+mom_L2_z;
  //  hlab[3]=sqrt(pow(hlab[0],2)+pow(hlab[1],2)+pow(hlab[2],2)+pow(m_Lambda->GetPDGMass()/CLHEP::GeV,2)*2);
  //  mass_hdibaryon = 2.250;    // Mass for H
  hlab[3]=sqrt(hlab[0]*hlab[0]+hlab[1]*hlab[1]+hlab[2]*hlab[2]+mass_hdibaryon*mass_hdibaryon);
  // G4double hcm[4];

  // G4double hcmpp=sqrt(pow(hcm[0],2)+pow(hcm[1],2)+pow(hcm[2],2));
  // cmh->Fill(acos(hcm[2]/hcmpp)/3.141592654*180);

  // coscmh->Fill(hcm[2]/hcmpp);

  //  G4double checkphi;
  // G4double hphi=atan2(hcm[1],hcm[0])/3.141592654*180;
  //  G4double hphilab=atan2(hlab[1],hlab[0])/3.141592654*180;

  // phih->Fill(hphi);
  // thetadiff->Fill(acos(hcm[2]/hcmpp)/3.141592654*180+acos(momcmk[2]/momcmkpp)/3.141592654*180);

  // phidiff->Fill(hphi-cmphik);
  //  G4cout<<"theta diff H:"<<acos(hcm[2]/hcmpp)/3.141592654*180<<G4endl;
  //  G4cout<<"theta diff K+:"<<acos(momcmk[2]/momcmkpp)/3.141592654*180<<G4endl;
  // G4double momsum[3]={0};
  // momsum[0]=mom_kp_x+mom_L1_x+mom_L2_x;
  // momsum[1]=mom_kp_y+mom_L1_y+mom_L2_y;
  // momsum[2]=mom_kp_z+mom_L1_z+mom_L2_z;
  //  G4cout<<"mom sum:"<<momsum[0]<<":"<<momsum[1]<<":"<<momsum[2]<<G4endl;
  // G4double momH[3]={0};
  // momH[0]=mom_L1_x+mom_L2_x;
  // momH[1]=mom_L1_y+mom_L2_y;
  // momH[2]=mom_L1_z+mom_L2_z;
  //  G4cout<<"mom sum X H:kp:"<<momH[0]<<":"<<mom_kp_x<<G4endl;
  //  G4cout<<"mom sum Y H:kp:"<<momH[1]<<":"<<mom_kp_y<<G4endl;
  //  G4cout<<"mom sum Z H:kp:"<<momH[2]<<":"<<mom_kp_z<<G4endl;

  // G4double e1,e2,etot,invm2,ptot,invm;
  // invm=-1.;

  // ptot=pow(mom_L1_x+mom_L2_x,2)+pow(mom_L1_y+mom_L2_y,2)+pow(mom_L1_z+mom_L2_z,2);
  // e1=Energy_L1;
  // e2=Energy_L2;
  // etot=pow(e1+e2,2);
  // invm2=etot-ptot;
  // if(invm2 > 0) invm=sqrt(invm2);
  //gen_im->Fill(invm);

  //  G4cout<<"IM(LL)"<<invm<<G4endl;
  //  G4cout<<"-----------end-------------"<<G4endl;
  //  G4double vtr = 20.0*CLHEP::mm*((double) G4RandFlat::shoot());
  //  G4double vtphi = 2.0*pi*((double) G4RandFlat::shoot());
  //  G4double vtx = vtr*cos(vtphi);
  //  G4double vty = vtr*sin(vtphi);

  G4double vtx = 0.*CLHEP::mm;
  G4double vty = 0.*CLHEP::mm;
  G4double vtz= G4RandFlat::shoot(m_target_pos.z()-m_target_size.z()/2,
				  m_target_pos.z()+m_target_size.z()/2)*CLHEP::mm;

  ///////kp PG
  m_particle_gun->SetParticleDefinition(m_KaonPlus);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_kp_x,mom_kp_y,mom_kp_z));
  m_particle_gun->SetParticleEnergy((Energy_kp - m_KaonPlus->GetPDGMass()/CLHEP::GeV)*CLHEP::GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex(anEvent);
  ///////L1 PG
  m_particle_gun->SetParticleDefinition(m_Lambda);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_L1_x,mom_L1_y,mom_L1_z));
  m_particle_gun->SetParticleEnergy((Energy_L1 - m_Lambda->GetPDGMass()/CLHEP::GeV)*CLHEP::GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex(anEvent);
  ///////L2 PG
  m_particle_gun->SetParticleDefinition(m_Lambda);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_L2_x,mom_L2_y,mom_L2_z));
  m_particle_gun->SetParticleEnergy((Energy_L2 - m_Lambda->GetPDGMass()/CLHEP::GeV)*CLHEP::GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex(anEvent);

  // gAnaMan.SetPrimaryParticle(0,mom_kp_x,mom_kp_y,mom_kp_z,m_KaonPlus->GetPDGMass()/CLHEP::GeV);
  // gAnaMan.SetPrimaryParticle(1,mom_L1_x,mom_L1_y,mom_L1_z,m_Lambda->GetPDGMass()/CLHEP::GeV);
  // gAnaMan.SetPrimaryParticle(2,mom_L2_x,mom_L2_y,mom_L2_z,m_Lambda->GetPDGMass()/CLHEP::GeV);
  // gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  // gAnaMan.SetPrimaryVertex(1,vtx,vty,vtz);
  // gAnaMan.SetPrimaryVertex(2,vtx,vty,vtz);
}

//_____________________________________________________________________________
void
PrimaryGeneratorAction::GenerateE42HdibaryonPHSG(G4Event* anEvent)
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
  //  Ebeam = sqrt(pimom*pimom+m_PionMinus->GetPDGMass()/CLHEP::GeV*m_PionMinus->GetPDGMass()/CLHEP::GeV);
  Ebeam = sqrt(pbeam*pbeam+m_KaonMinus->GetPDGMass()/CLHEP::GeV*m_KaonMinus->GetPDGMass()/CLHEP::GeV);

  pbm[0]=0;
  pbm[1]=0;
  pbm[2]=pbeam;
  pbm[3]=Ebeam;
 up:
  KinemaHweak Hkinema(m_KaonMinus->GetPDGMass()/CLHEP::GeV,
                      (m_Proton->GetPDGMass()/CLHEP::GeV)*2,
                      m_Hdibaryon->GetPDGMass()/CLHEP::GeV,
                      m_KaonPlus->GetPDGMass()/CLHEP::GeV,
                      pbm[2], 0.0);

  Energy_kp=Hkinema.GetEnergy(4);
  // G4double momentum_kp = Hkinema.GetMomentum(4);
  Hkinema.GetMomentum(4,mom);
  mom_kp_x = mom[0];
  mom_kp_y = mom[1];
  mom_kp_z = mom[2];
  if((acos(mom[2]/sqrt(pow(mom[0],2)+pow(mom[1],2)+pow(mom[2],2)))*180/3.141592654)>15.) goto up;
  //  G4cout<<"kp mom:"<<mom[1]<<G4endl;

  Energy_h = Hkinema.GetEnergy(3);
  // G4double momentum_h = Hkinema.GetMomentum(3);
  Hkinema.GetMomentum(3,mom);

  //  G4cout<<"mom kp :"<<mom[0]<<":"<<mom[1]<<":"<<mom[2]<<G4endl;
  //  G4cout<<"mom comp:"<<momentum_h<<":"<<sqrt(mom[0]*mom[0]+mom[1]*mom[1]+mom[2]*mom[2])<<G4endl;
  mom_h_x = mom[0];
  mom_h_y = mom[1];
  mom_h_z = mom[2];

  //  G4cout<<"h mom:"<<mom[1]<<G4endl;

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
  //  G4double beta= pbeam/(2.*0.938272013+Ebeam);

  //  G4cout<<"beta:"<<beta<<G4endl;
  // G4double momcmkpp=sqrt(pow(momcmk[0],2)+pow(momcmk[1],2)+pow(momcmk[2],2));
  //  G4cout<<"mom comp-----------------------"<<G4endl;
  //  G4cout<<"mom comp:"<<momk[2]<<":"<<momcmk[2]<<G4endl;
  //  G4cout<<"mom comp:"<<momk[3]<<":"<<momcmk[3]<<G4endl;
  //cmk->Fill(acos(momcmk[2]/momcmkpp)/3.141592654*180);

  //  G4double cmphik=atan2(momk[1],momk[0])*180/3.141592654;

  //  G4cout<<"comp phi:"<<cmphik<<":"<<Phi_h<<G4endl;

  //coscmk->Fill(momcmk[2]/momcmkpp);


  G4double vtx = 0.*CLHEP::mm;
  G4double vty = 0.*CLHEP::mm;
  G4double vtz = G4RandFlat::shoot(m_target_pos.z() - m_target_size.z()/2,
                                   m_target_pos.z() + m_target_size.z()/2)*CLHEP::mm;
  //  G4double vtz= G4RandFlat::shoot(-150.-env_target_width/2,-150.+env_target_width/2);

  //Kaon +
  m_particle_gun->SetParticleDefinition(m_KaonPlus);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_kp_x,mom_kp_y,mom_kp_z));
  m_particle_gun->SetParticleEnergy((Energy_kp - m_KaonPlus->GetPDGMass()/CLHEP::GeV)*CLHEP::GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex(anEvent);

  //Hdibaryon
  m_particle_gun->SetParticleDefinition(m_Hdibaryon);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_h_x,mom_h_y,mom_h_z));
  m_particle_gun->SetParticleEnergy((Energy_h - m_Hdibaryon->GetPDGMass()/CLHEP::GeV)*CLHEP::GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex(anEvent);

  // gAnaMan.SetPrimaryParticle(0,mom_kp_x,mom_kp_y,mom_kp_z,m_KaonPlus->GetPDGMass()/CLHEP::GeV);
  // gAnaMan.SetPrimaryParticle(1,mom_h_x,mom_h_y,mom_h_z,m_Hdibaryon->GetPDGMass()/CLHEP::GeV);
  // gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  // gAnaMan.SetPrimaryVertex(1,vtx,vty,vtz);
}

//_____________________________________________________________________________
void
PrimaryGeneratorAction::GenerateE42HdibaryonPHSGS(G4Event* anEvent)
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
  //  Ebeam = sqrt(pimom*pimom+m_PionMinus->GetPDGMass()/CLHEP::GeV*m_PionMinus->GetPDGMass()/CLHEP::GeV);
  Ebeam = sqrt(pbeam*pbeam+m_KaonMinus->GetPDGMass()/CLHEP::GeV*m_KaonMinus->GetPDGMass()/CLHEP::GeV);

  pbm[0]=0;
  pbm[1]=0;
  pbm[2]=pbeam;
  pbm[3]=Ebeam;
  //  G4cout<<"m_Hdibaryon mass"<<m_Hdibaryon->GetPDGMass()/CLHEP::GeV<<G4endl;
 up:
  KinemaHweak Hkinema(m_KaonMinus->GetPDGMass()/CLHEP::GeV,
		      (m_Proton->GetPDGMass()/CLHEP::GeV)*2,
		      m_HdibaryonS->GetPDGMass()/CLHEP::GeV,
		      m_KaonPlus->GetPDGMass()/CLHEP::GeV,
		      pbm[2], 0.0);

  Energy_kp=Hkinema.GetEnergy(4);
  // G4double momentum_kp = Hkinema.GetMomentum(4);
  Hkinema.GetMomentum(4,mom);
  mom_kp_x = mom[0];
  mom_kp_y = mom[1];
  mom_kp_z = mom[2];

  if((acos(mom[2]/sqrt(pow(mom[0],2)+pow(mom[1],2)+pow(mom[2],2)))*180/3.141592654)>15.) goto up;
  //  G4cout<<"kp mom:"<<mom[1]<<G4endl;

  Energy_h = Hkinema.GetEnergy(3);
  // G4double momentum_h = Hkinema.GetMomentum(3);
  Hkinema.GetMomentum(3,mom);

  //  G4cout<<"mom kp :"<<mom[0]<<":"<<mom[1]<<":"<<mom[2]<<G4endl;
  //  G4cout<<"mom comp:"<<momentum_h<<":"<<sqrt(mom[0]*mom[0]+mom[1]*mom[1]+mom[2]*mom[2])<<G4endl;

  mom_h_x = mom[0];
  mom_h_y = mom[1];
  mom_h_z = mom[2];

  //  G4cout<<"h mom:"<<mom[1]<<G4endl;

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
  //  G4double beta= pbeam/(2.*0.938272013+Ebeam);

  //  G4cout<<"beta:"<<beta<<G4endl;
  // G4double momcmkpp=sqrt(pow(momcmk[0],2)+pow(momcmk[1],2)+pow(momcmk[2],2));
  //  G4cout<<"mom comp-----------------------"<<G4endl;
  //  G4cout<<"mom comp:"<<momk[2]<<":"<<momcmk[2]<<G4endl;
  //  G4cout<<"mom comp:"<<momk[3]<<":"<<momcmk[3]<<G4endl;
  //cmk->Fill(acos(momcmk[2]/momcmkpp)/3.141592654*180);

  //  G4double cmphik=atan2(momk[1],momk[0])*180/3.141592654;

  //  G4cout<<"comp phi:"<<cmphik<<":"<<Phi_h<<G4endl;

  //coscmk->Fill(momcmk[2]/momcmkpp);


  G4double vtx = 0.*CLHEP::mm;
  G4double vty = 0.*CLHEP::mm;
  //  G4cout<<"env_target_width :"<<atof(env_env_target_width.c_str())<<G4endl;
  //  G4double vtz= G4RandFlat::shoot(-150.-env_target_width/2,-150.+env_target_width/2);
  G4double vtz= G4RandFlat::shoot(m_target_pos.z() - m_target_size.z()/2,
                                  m_target_pos.z() + m_target_size.z()/2)*CLHEP::mm;

  //Kaon +
  m_particle_gun->SetParticleDefinition(m_KaonPlus);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_kp_x,mom_kp_y,mom_kp_z));
  m_particle_gun->SetParticleEnergy((Energy_kp - m_KaonPlus->GetPDGMass()/CLHEP::GeV)*CLHEP::GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex(anEvent);

  //Hdibaryon
  m_particle_gun->SetParticleDefinition(m_HdibaryonS);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_h_x,mom_h_y,mom_h_z));
  m_particle_gun->SetParticleEnergy((Energy_h - m_HdibaryonS->GetPDGMass()/CLHEP::GeV)*CLHEP::GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex(anEvent);

  // gAnaMan.SetPrimaryParticle(0,mom_kp_x,mom_kp_y,mom_kp_z,m_KaonPlus->GetPDGMass()/CLHEP::GeV);
  // gAnaMan.SetPrimaryParticle(1,mom_h_x,mom_h_y,mom_h_z,m_HdibaryonS->GetPDGMass()/CLHEP::GeV);
  // gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  // gAnaMan.SetPrimaryVertex(1,vtx,vty,vtz);
}

//_____________________________________________________________________________
void
PrimaryGeneratorAction::GenerateE42HdibaryonPHSGLL(G4Event* anEvent)
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
  //  Ebeam = sqrt(pimom*pimom+m_PionMinus->GetPDGMass()/CLHEP::GeV*m_PionMinus->GetPDGMass()/CLHEP::GeV);
  Ebeam = sqrt(pbeam*pbeam+m_KaonMinus->GetPDGMass()/CLHEP::GeV*m_KaonMinus->GetPDGMass()/CLHEP::GeV);

  pbm[0]=0;
  pbm[1]=0;
  pbm[2]=pbeam;
  pbm[3]=Ebeam;
 up:
  KinemaHweak Hkinema(m_KaonMinus->GetPDGMass()/CLHEP::GeV,
                      (m_Proton->GetPDGMass()/CLHEP::GeV)*2,
                      m_HdibaryonLL->GetPDGMass()/CLHEP::GeV,
                      m_KaonPlus->GetPDGMass()/CLHEP::GeV,
                      pbm[2], 0.0);

  Energy_kp=Hkinema.GetEnergy(4);
  // G4double momentum_kp = Hkinema.GetMomentum(4);
  Hkinema.GetMomentum(4,mom);
  mom_kp_x = mom[0];
  mom_kp_y = mom[1];
  mom_kp_z = mom[2];


  if((acos(mom[2]/sqrt(pow(mom[0],2)+pow(mom[1],2)+pow(mom[2],2)))*180/3.141592654)>15. ) goto up;
  //  G4cout<<"kp mom:"<<mom[1]<<G4endl;

  Energy_h = Hkinema.GetEnergy(3);
  // G4double momentum_h = Hkinema.GetMomentum(3);
  Hkinema.GetMomentum(3,mom);

  //  G4cout<<"mom kp :"<<mom[0]<<":"<<mom[1]<<":"<<mom[2]<<G4endl;
  //  G4cout<<"mom comp:"<<momentum_h<<":"<<sqrt(mom[0]*mom[0]+mom[1]*mom[1]+mom[2]*mom[2])<<G4endl;

  mom_h_x = mom[0];
  mom_h_y = mom[1];
  mom_h_z = mom[2];

  //  G4cout<<"h mom:"<<mom[1]<<G4endl;

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
  //  G4double beta= pbeam/(2.*0.938272013+Ebeam);

  //  G4cout<<"beta:"<<beta<<G4endl;
  // G4double momcmkpp=sqrt(pow(momcmk[0],2)+pow(momcmk[1],2)+pow(momcmk[2],2));
  //  G4cout<<"mom comp-----------------------"<<G4endl;
  //  G4cout<<"mom comp:"<<momk[2]<<":"<<momcmk[2]<<G4endl;
  //  G4cout<<"mom comp:"<<momk[3]<<":"<<momcmk[3]<<G4endl;
  //cmk->Fill(acos(momcmk[2]/momcmkpp)/3.141592654*180);

  //  G4double cmphik=atan2(momk[1],momk[0])*180/3.141592654;

  //  G4cout<<"comp phi:"<<cmphik<<":"<<Phi_h<<G4endl;

  //coscmk->Fill(momcmk[2]/momcmkpp);

  G4double vtx = 0.*CLHEP::mm;
  G4double vty = 0.*CLHEP::mm;
  G4double vtz = G4RandFlat::shoot(m_target_pos.z() - m_target_size.z()/2,
                                   m_target_pos.z() + m_target_size.z()/2)*CLHEP::mm;

  //Kaon +
  m_particle_gun->SetParticleDefinition(m_KaonPlus);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_kp_x,mom_kp_y,mom_kp_z));
  m_particle_gun->SetParticleEnergy((Energy_kp - m_KaonPlus->GetPDGMass()/CLHEP::GeV)*CLHEP::GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex(anEvent);

  //Hdibaryon
  m_particle_gun->SetParticleDefinition(m_HdibaryonLL);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_h_x,mom_h_y,mom_h_z));
  m_particle_gun->SetParticleEnergy((Energy_h - m_HdibaryonLL->GetPDGMass()/CLHEP::GeV)*CLHEP::GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex(anEvent);

  // gAnaMan.SetPrimaryParticle(0,mom_kp_x,mom_kp_y,mom_kp_z,m_KaonPlus->GetPDGMass()/CLHEP::GeV);
  // gAnaMan.SetPrimaryParticle(1,mom_h_x,mom_h_y,mom_h_z,m_HdibaryonLL->GetPDGMass()/CLHEP::GeV);
  // gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  // gAnaMan.SetPrimaryVertex(1,vtx,vty,vtz);
}

//_____________________________________________________________________________
void
PrimaryGeneratorAction::GenerateE42HdibaryonNonReso(G4Event* anEvent)
{
  G4double  momk[4]={0.}, mom[3]={0.};
  G4double Ebeam, pbeam;
  //  G4double rmk=0.493677;
  G4double pbm[4];
  G4double Energy_h, mom_h_x, mom_h_y, mom_h_z;
  G4double Energy_kp, mom_kp_x, mom_kp_y, mom_kp_z;

  //  G4double pimom=0.635+G4RandFlat::shoot()*(2.000-0.635);
  pbeam=1.8;
  // gAnaMan.SetPrimaryBeam(0,0,pbeam);
  //  Ebeam = sqrt(pimom*pimom+m_PionMinus->GetPDGMass()/CLHEP::GeV*m_PionMinus->GetPDGMass()/CLHEP::GeV);
  Ebeam = sqrt(pbeam*pbeam+m_KaonMinus->GetPDGMass()/CLHEP::GeV*m_KaonMinus->GetPDGMass()/CLHEP::GeV);

  pbm[0]=0;
  pbm[1]=0;
  pbm[2]=pbeam;
  pbm[3]=Ebeam;
 up:
  KinemaHweak Hkinema(m_KaonMinus->GetPDGMass()/CLHEP::GeV,
                      (m_Proton->GetPDGMass()/CLHEP::GeV)*2,
                      m_Hdibaryon->GetPDGMass()/CLHEP::GeV,
                      m_KaonPlus->GetPDGMass()/CLHEP::GeV,
                      pbm[2], 0.0);

  Energy_kp=Hkinema.GetEnergy(4);
  // G4double momentum_kp = Hkinema.GetMomentum(4);
  Hkinema.GetMomentum(4,mom);
  mom_kp_x = mom[0];
  mom_kp_y = mom[1];
  mom_kp_z = mom[2];


  if((acos(momk[2]/sqrt(pow(momk[0],2)+pow(momk[1],2)+pow(momk[2],2)))*180/3.141592654)>15.) goto up;
  //  G4cout<<"kp mom:"<<mom[1]<<G4endl;

  Energy_h = Hkinema.GetEnergy(3);
  // G4double momentum_h = Hkinema.GetMomentum(3);
  Hkinema.GetMomentum(3,mom);

  //  G4cout<<"mom kp :"<<mom[0]<<":"<<mom[1]<<":"<<mom[2]<<G4endl;
  //  G4cout<<"mom comp:"<<momentum_h<<":"<<sqrt(mom[0]*mom[0]+mom[1]*mom[1]+mom[2]*mom[2])<<G4endl;

  mom_h_x = mom[0];
  mom_h_y = mom[1];
  mom_h_z = mom[2];

  //  G4cout<<"h mom:"<<mom[1]<<G4endl;

  // G4double Theta_h = Hkinema.GetTheta(3);
  // G4double Phi_h = Hkinema.GetPhi(3);

  ////check kinematics
  momk[0]=mom_kp_x;
  momk[1]=mom_kp_y;
  momk[2]=mom_kp_z;
  momk[3]=sqrt(mom_kp_x*mom_kp_x+mom_kp_y*mom_kp_y+mom_kp_z*mom_kp_z+0.493677*0.493677);
  //  G4double momkpp=sqrt(pow(mom_kp_x,2)+pow(mom_kp_y,2)+pow(mom_kp_z,2));
  // G4double momcmk[4]={0.};
  //  G4double beta= pbeam/(2.*0.938272013+Ebeam);

  //  G4cout<<"beta:"<<beta<<G4endl;
  // G4double momcmkpp=sqrt(pow(momcmk[0],2)+pow(momcmk[1],2)+pow(momcmk[2],2));
  //  G4cout<<"mom comp-----------------------"<<G4endl;
  //  G4cout<<"mom comp:"<<momk[2]<<":"<<momcmk[2]<<G4endl;
  //  G4cout<<"mom comp:"<<momk[3]<<":"<<momcmk[3]<<G4endl;
  //cmk->Fill(acos(momcmk[2]/momcmkpp)/3.141592654*180);

  //  G4double cmphik=atan2(momk[1],momk[0])*180/3.141592654;

  //  G4cout<<"comp phi:"<<cmphik<<":"<<Phi_h<<G4endl;

  //coscmk->Fill(momcmk[2]/momcmkpp);

  G4double vtx = 0.*CLHEP::mm;
  G4double vty = 0.*CLHEP::mm;
  G4double vtz = G4RandFlat::shoot(m_target_pos.z() - m_target_size.z()/2,
                                   m_target_pos.z() + m_target_size.z()/2)*CLHEP::mm;

  //Kaon +
  m_particle_gun->SetParticleDefinition(m_KaonPlus);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_kp_x,mom_kp_y,mom_kp_z));
  m_particle_gun->SetParticleEnergy((Energy_kp - m_KaonPlus->GetPDGMass()/CLHEP::GeV)*CLHEP::GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex(anEvent);

  //Hdibaryon
  m_particle_gun->SetParticleDefinition(m_Hdibaryon);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_h_x,mom_h_y,mom_h_z));
  m_particle_gun->SetParticleEnergy((Energy_h - m_Hdibaryon->GetPDGMass()/CLHEP::GeV)*CLHEP::GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex(anEvent);

  // gAnaMan.SetPrimaryParticle(0,mom_kp_x,mom_kp_y,mom_kp_z,m_KaonPlus->GetPDGMass()/CLHEP::GeV);
  // gAnaMan.SetPrimaryParticle(1,mom_h_x,mom_h_y,mom_h_z,m_Hdibaryon->GetPDGMass()/CLHEP::GeV);
  // gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  // gAnaMan.SetPrimaryVertex(1,vtx,vty,vtz);
}

//_____________________________________________________________________________
void
PrimaryGeneratorAction::GenerateE42Hdibaryon1(G4Event* anEvent)
{
  G4double mass_hdibaryon, width_hdibaryon;

  G4double Energy_L1, mom_L1_x, mom_L1_y, mom_L1_z;
  G4double Energy_L2, mom_L2_x, mom_L2_y, mom_L2_z;
  G4double Energy_kp, mom_kp_x, mom_kp_y, mom_kp_z;

  G4double mom[3];
  G4double pg_x,pg_y,pg_z;
  // G4double rmk=0.493677;

  // G4double Ebeam = sqrt(1.8*1.8+m_KaonMinus->GetPDGMass()/CLHEP::GeV*m_KaonMinus->GetPDGMass()/CLHEP::GeV);       // Incident gamma energy (CLHEP::GeV)

  pg_x = 0.0;
  pg_y = 0.0;
  pg_z = 1.8;
  // G4double pbeam=1.8;

  // gAnaMan.SetPrimaryBeam(pg_x,pg_y,pg_z);

  mass_hdibaryon = gConf.Get<G4double>("HdibaryonMass");
  width_hdibaryon = gConf.Get<G4double>("HdibaryonWidth");

  Kinema3Resonance Hkinema(m_KaonMinus->GetPDGMass()/CLHEP::GeV,
                           ((m_Proton->GetPDGMass()/CLHEP::GeV)*2),
                           m_Lambda->GetPDGMass()/CLHEP::GeV,
                           m_Lambda->GetPDGMass()/CLHEP::GeV,
                           m_KaonPlus->GetPDGMass()/CLHEP::GeV,
                           mass_hdibaryon, width_hdibaryon, pg_z, 0.0);

  //  Hkinema.Dump();

  Energy_kp = Hkinema.GetEnergy(5);
  // G4double   momentum_kp = Hkinema.GetMomentum(5);
  Hkinema.GetMomentum(5,mom);

  /* L1 */
  Energy_L1 = Hkinema.GetEnergy(3);
  // G4double   momentum_L1 = Hkinema.GetMomentum(3);
  Hkinema.GetMomentum(3,mom);

  mom_L1_x = mom[0];
  mom_L1_y = mom[1];
  mom_L1_z = mom[2];

  // G4double ThetaL1 = Hkinema.GetTheta(3);
  // G4double PhiL1 = Hkinema.GetPhi(3);

  /* L2 */
  Energy_L2 = Hkinema.GetEnergy(4);
  // G4double   momentum_L2 = Hkinema.GetMomentum(4);
  Hkinema.GetMomentum(4,mom);

  mom_L2_x = mom[0];
  mom_L2_y = mom[1];
  mom_L2_z = mom[2];

  // G4double   ThetaL2 = Hkinema.GetTheta(4);
  // G4double   PhiL2 = Hkinema.GetPhi(4);

  /* Kp */
  Energy_kp = Hkinema.GetEnergy(5);
  // G4double   momentum_kp = Hkinema.GetMomentum(5);
  Hkinema.GetMomentum(5,mom);

  mom_kp_x = mom[0];
  mom_kp_y = mom[1];
  mom_kp_z = mom[2];

  // G4double Thetakp = Hkinema.GetThetaCM(1);
  // G4double shphi=(atan2(mom_kp_y,mom_kp_x))/3.141592654*180.;
  //phik->Fill(shphi);

  // G4double beta= sqrt(pbeam*pbeam)/(2.*0.938272013+Ebeam);
  // G4double momk[4]={};
  // momk[0]=mom_kp_x;
  // momk[1]=mom_kp_y;
  // momk[2]=mom_kp_z;
  // momk[3]=sqrt(mom_kp_x*mom_kp_x+mom_kp_y*mom_kp_y+mom_kp_z*mom_kp_z+0.493677*0.493677);
  // G4double momkpp=sqrt(pow(mom_kp_x,2)+pow(mom_kp_y,2)+pow(mom_kp_z,2));
  // G4double momcmkpp=sqrt(pow(momcmk[0],2)+pow(momcmk[1],2)+pow(momcmk[2],2));
  //cmk->Fill(acos(momcmk[2]/momcmkpp)/3.141592654*180);

  // coscmk->Fill((momcmk[2]/momcmkpp));

  // labk->Fill(acos(momk[2]/momkpp)*180/3.141592654);
  // coslabk->Fill(momk[2]/momkpp);

  //cal for invariant mass
  // G4double ptot=pow(mom_L1_x+mom_L2_x,2)+pow(mom_L1_y+mom_L2_y,2)+pow(mom_L1_z+mom_L2_z,2);
  // G4double e1=Energy_L1;
  // G4double e2=Energy_L2;
  // G4double etot=pow(e1+e2,2);
  // G4double invm2=etot-ptot;
  // if(invm2 > 0) invm=sqrt(invm2);
  //gen_im->Fill(invm);
  double vtx = 0.*CLHEP::mm;
  double vty = 0.*CLHEP::mm;
  /////shhwang hdibaryon1

  //  G4cout<<"env_target_width :"<<atof(env_env_target_width.c_str())<<G4endl;
  //  G4double vtz= G4RandFlat::shoot(-150.-env_target_width/2,-150.+env_target_width/2);
  G4double vtz= G4RandFlat::shoot(m_target_pos.z() - m_target_size.z()/2,
                                  m_target_pos.z()+m_target_size.z()/2)*CLHEP::mm;

  ///////kp PG
  m_particle_gun->SetParticleDefinition(m_KaonPlus);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_kp_x,mom_kp_y,mom_kp_z));
  m_particle_gun->SetParticleEnergy((Energy_kp - m_KaonPlus->GetPDGMass()/CLHEP::GeV)*CLHEP::GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex(anEvent);
  ///////L1 PG
  m_particle_gun->SetParticleDefinition(m_Lambda);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_L1_x,mom_L1_y,mom_L1_z));
  m_particle_gun->SetParticleEnergy((Energy_L1 - m_Lambda->GetPDGMass()/CLHEP::GeV)*CLHEP::GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex(anEvent);
  ///////L2 PG
  m_particle_gun->SetParticleDefinition(m_Lambda);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_L2_x,mom_L2_y,mom_L2_z));
  m_particle_gun->SetParticleEnergy((Energy_L2 - m_Lambda->GetPDGMass()/CLHEP::GeV)*CLHEP::GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex(anEvent);

  // gAnaMan.SetPrimaryParticle(0,mom_kp_x,mom_kp_y,mom_kp_z, m_KaonPlus->GetPDGMass()/CLHEP::GeV);
  // gAnaMan.SetPrimaryParticle(1,mom_L1_x,mom_L1_y,mom_L1_z, m_Lambda->GetPDGMass()/CLHEP::GeV);
  // gAnaMan.SetPrimaryParticle(2,mom_L2_x,mom_L2_y,mom_L2_z, m_Lambda->GetPDGMass()/CLHEP::GeV);
  // gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
  // gAnaMan.SetPrimaryVertex(1,vtx,vty,vtz);
  // gAnaMan.SetPrimaryVertex(2,vtx,vty,vtz);
}

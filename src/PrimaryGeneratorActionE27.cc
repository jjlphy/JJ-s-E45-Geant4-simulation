// -*- C++ -*-

#include "PrimaryGeneratorAction.hh"

#include <G4Event.hh>
#include <G4IonConstructor.hh>
#include <G4IonTable.hh>
#include <G4KaonPlus.hh>
#include <G4LorentzVector.hh>
#include <G4ParticleGun.hh>
#include <G4ParticleGun.hh>
#include <G4ParticleDefinition.hh>
#include <G4ParticleTypes.hh>
#include <G4UImanager.hh>
#include <Randomize.hh>

#include "AnaManager.hh"
#include "AngDisGenerator.hh"
#include "ConfMan.hh"
#include "DCGeomMan.hh"
#include "DetSizeMan.hh"
#include "FuncName.hh"
#include "GeneratorHelper.hh"
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
const int MaxTry = 1000;
const double AtomicMassUnit = 0.9314932;
const auto& gConf = ConfMan::GetInstance();
const auto& gGeom = DCGeomMan::GetInstance();
const auto& gSize = DetSizeMan::GetInstance();
const auto particleTable = G4ParticleTable::GetParticleTable();
}

//_____________________________________________________________________________
// reactio No #2701 pi+ beam through
void
PrimaryGeneratorAction::GenerateE27BeamThrough(G4Event* anEvent)
{
  //  G4double  momk[3], mom[3],momkn[3];
  //  G4double rmk=0.493677;
  //  G4double pbm[4];
  G4double Energy_pip,  mom_pip_x, mom_pip_y, mom_pip_z;

  auto pionPlus = G4PionPlus::Definition();
  G4double pbeam = CLHEP::RandGauss::shoot(gConf.Get<G4double>("BeamMom"),
                                           0.01294*gConf.Get<G4double>("BeamMom"));
  //  pbeam=CLHEP::RandGauss::shoot(env_Beam_mom,env_Beam_mom*3.3*0.0001/2.3548);
  //  pbeam=1.8;
  mom_pip_x=0;
  mom_pip_y=0;
  mom_pip_z=pbeam;
  Energy_pip=sqrt(pionPlus->GetPDGMass()/CLHEP::GeV*pionPlus->GetPDGMass()/CLHEP::GeV+pbeam*pbeam);

  // gAnaMan.SetPrimaryBeam(0,0,pbeam);

  // G4double Ebeam = sqrt(pbeam*pbeam+pionPlus->GetPDGMass()/CLHEP::GeV*pionPlus->GetPDGMass()/CLHEP::GeV);

  //  G4double vtx = CLHEP::RandFlat::shoot(-15.,15.)*CLHEP::mm;
  //  G4double vty = CLHEP::RandFlat::shoot(-5.,5.)*CLHEP::mm;
  //  G4double vtz= CLHEP::RandFlat::shoot(env_target_pos_z-env_target_width/2,env_target_pos_z+env_target_width/2)*CLHEP::mm;

  G4double vtx = CLHEP::RandGauss::shoot(0,10.)*CLHEP::mm;
  G4double vty = CLHEP::RandFlat::shoot(0.,3.2)*CLHEP::mm;
  G4double vtz = gSize.Get("Target", ThreeVector::Z);
  //G4double vtz= CLHEP::RandFlat::shoot(m_particle_gun->Get_env_target_pos_z()-gSize.Get("Target", ThreeVector::Z)/2,m_particle_gun->Get_env_target_pos_z()+gSize.Get("Target", ThreeVector::Z)/2)*CLHEP::mm-250.*CLHEP::mm;
  std::cout<<"pbeam = "<<pbeam<<std::endl;
  //getchar();
  //beam pi+
  m_particle_gun->SetParticleDefinition(pionPlus);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_pip_x,mom_pip_y,mom_pip_z));
  m_particle_gun->SetParticleEnergy((Energy_pip - pionPlus->GetPDGMass()/CLHEP::GeV)*CLHEP::GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex(anEvent);
  //  gAnaMan->SetPrimaryParticle(0,mom_pip_x,mom_pip_y,mom_pip_z,pionPlus->GetPDGMass()/CLHEP::GeV);
  //  gAnaMan->SetPrimaryVertex(0,vtx,vty,vtz);
}

//_____________________________________________________________________________
//reactio No #2702 K+ gun for test
void
PrimaryGeneratorAction::GenerateE27Kptest(G4Event* anEvent)
{
  //  G4double  momk[3], mom[3],momkn[3];
  //  G4double rmk=0.493677;
  //  G4double pbm[4];
  G4double Energy_Kp,  mom_Kp_x, mom_Kp_y, mom_Kp_z;

  auto KaonPlus = G4KaonPlus::Definition();
  //  kaonMinus = particleTable->FindParticle("kaon-");
  //kaonMinus = particleTable->FindParticle("pi-");
  G4double pbeam=CLHEP::RandGauss::shoot(gConf.Get<G4double>("BeamMom"),0.01294*gConf.Get<G4double>("BeamMom"));
  //  pbeam=CLHEP::RandGauss::shoot(env_Beam_mom,env_Beam_mom*3.3*0.0001/2.3548);
  //  pbeam=1.8;
  mom_Kp_x=0;
  mom_Kp_y=0;
  mom_Kp_z=pbeam;
  Energy_Kp=sqrt(KaonPlus->GetPDGMass()/CLHEP::GeV*KaonPlus->GetPDGMass()/CLHEP::GeV+pbeam*pbeam);

  // gAnaMan.SetPrimaryBeam(0,0,pbeam);

  // G4double Ebeam = sqrt(pbeam*pbeam+KaonPlus->GetPDGMass()/CLHEP::GeV*KaonPlus->GetPDGMass()/CLHEP::GeV);

  //  G4double vtx = CLHEP::RandFlat::shoot(-15.,15.)*CLHEP::mm;
  //  G4double vty = CLHEP::RandFlat::shoot(-5.,5.)*CLHEP::mm;
  //  G4double vtz= CLHEP::RandFlat::shoot(env_target_pos_z-env_target_width/2,env_target_pos_z+env_target_width/2)*CLHEP::mm;

  // G4double vtx = CLHEP::RandGauss::shoot(0,10.)*CLHEP::mm;
  // G4double vty = CLHEP::RandFlat::shoot(0.,3.2)*CLHEP::mm;

  G4double vtx = 0.*CLHEP::mm;
  G4double vty = 0.*CLHEP::mm;
  //  G4double vtz= CLHEP::RandFlat::shoot(gGeom.GetGlobalPosition("Target").z()-gSize.Get("Target", ThreeVector::Z)/2,gGeom.GetGlobalPosition("Target").z()+gSize.Get("Target", ThreeVector::Z)/2)*CLHEP::mm-250.*CLHEP::mm;
  G4double vtz= gGeom.GetGlobalPosition("Target").z();
  std::cout<<"pbeam = "<<pbeam<<std::endl;
  // getchar();
  //scat K+
  m_particle_gun->SetParticleDefinition(KaonPlus);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_Kp_x,mom_Kp_y,mom_Kp_z));
  m_particle_gun->SetParticleEnergy((Energy_Kp - KaonPlus->GetPDGMass()/CLHEP::GeV)*CLHEP::GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex(anEvent);
  //  gAnaMan->SetPrimaryParticle(0,mom_Kp_x,mom_Kp_y,mom_Kp_z,pionPlus->GetPDGMass()/CLHEP::GeV);
  //  gAnaMan->SetPrimaryVertex(0,vtx,vty,vtz);
}

//_____________________________________________________________________________
//reactio No #2703 pi+ d -> K+ K-pp, K-pp -> Lambda p reaction
void
PrimaryGeneratorAction::GenerateE27KppFLambdaP(G4Event* anEvent)
{
  G4double Mi1=G4PionPlus::Definition()->GetPDGMass();
  G4double Mi2=G4Deuteron::Definition()->GetPDGMass();
  G4double Mf1=G4KaonPlus::Definition()->GetPDGMass();
  G4double Mf2=G4Lambda::Definition()->GetPDGMass();
  G4double Mf3=G4Proton::Definition()->GetPDGMass();

  G4double Mp=G4Proton::Definition()->GetPDGMass();
  // G4double Mpim=G4PionMinus::Definition()->GetPDGMass();
  // G4double Mn=G4Neutron::Definition()->GetPDGMass();
  // G4double Mpiz=G4PionZero::Definition()->GetPDGMass();

  G4ThreeVector LPos = GaussPosition_LqTarg(gConf.Get<G4double>("BeamX0"),
                                            gConf.Get<G4double>("BeamY0"),
                                            gGeom.GetGlobalPosition("Target").z(),
                                            gConf.Get<G4double>("BeamDX"),
                                            gConf.Get<G4double>("BeamDY"),
                                            gSize.Get("Target", ThreeVector::X),
                                            gSize.Get("Target", ThreeVector::Z));
  //Note!! env_target_width = Target_Size_z (height of target)

  G4ThreeVector LBeamDir =  GaussDirectionInUV(gConf.Get<G4double>("BeamU0"),
                                               gConf.Get<G4double>("BeamV0"),
                                               gConf.Get<G4double>("BeamDU"),
                                               gConf.Get<G4double>("BeamDV"));

  G4double pb = gConf.Get<G4double>("BeamMom")*CLHEP::GeV;
  G4double dpb = 0.;
  if(gConf.Get<G4double>("BeamMom")!=0.)
    dpb = G4RandGauss::shoot(0.,gConf.Get<G4double>("BeamWidth"))*CLHEP::GeV;
  pb += dpb;

  G4double Mm1 = BreitWigner(2.275, 0.162)*CLHEP::GeV; //E27

  G4double thetaK, theta_CM;
  G4ThreeVector LPm1, LPf1, LPf2, LPf3;

  AGUniform gen1;
  //AGLambda1405 gen1;
  AGUniform gen2;
  G4ThreeVector LPini2=G4ThreeVector(0.,0.,0.);

  bool status = false;
  bool status2 = false;

  G4int n=0;
  while(1){
    if(++n>MaxTry){
      G4Exception(FUNC_NAME,
                  "Production under threshold",
                  RunMustBeAborted,
                  "Production under Threshold!!");
    }

    status=Scattering2Body_theta(Mi1, Mi2, Mf1, Mm1,
                                 pb*LBeamDir,LPini2,
                                 LPf1, LPm1,theta_CM, gen1);
    theta_CM = theta_CM*(180./(acos(-1.)));

    if(status ==true){
      thetaK = LPf1.theta()*(180./(acos(-1.)));
      G4cout<<"thetaK= " <<thetaK <<G4endl;
      if(thetaK<20.){
        status2=Decay2Body(Mm1, Mf2, Mf3, LPm1, LPf2, LPf3, gen2);
        if(status2 == true)
          break;
        else
          std::cout<<"Mm1="<<Mm1<<std::endl;
      }
    }
    pb = gConf.Get<G4double>("BeamMom")*CLHEP::GeV;
    dpb = 0.;
    if(gConf.Get<G4double>("BeamMom")!=0.)
      dpb = G4RandGauss::shoot(0.,gConf.Get<G4double>("BeamWidth"))*CLHEP::GeV;
    pb += dpb;
    Mm1 = BreitWigner(2.275, 0.162)*CLHEP::GeV; //E27
  }



  G4ThreeVector beam_mom = pb*LBeamDir;

  G4LorentzVector Lv_beam, Lv_targ_D, Lv_targ_P, Lv_K;
  Lv_beam.setVectM(pb*LBeamDir, Mi1);
  Lv_targ_D.setVectM(G4ThreeVector(0.,0.,0.), Mi2);
  Lv_targ_P.setVectM(G4ThreeVector(0.,0.,0.), Mp);
  Lv_K.setVectM(LPf1, Mf1);

  // gAnaMan.SetPrimaryBeam(beam_mom.x(),beam_mom.y(),beam_mom.z());
  auto KaonPlus = G4KaonPlus::Definition();
  m_particle_gun->SetParticleDefinition(KaonPlus);
  m_particle_gun->SetParticlePosition(LPos);
  m_particle_gun->SetParticleMomentum(LPf1);
  m_particle_gun->GeneratePrimaryVertex(anEvent);

  auto Lambda= G4Lambda::Definition();
  m_particle_gun->SetParticleDefinition(Lambda);
  m_particle_gun->SetParticlePosition(LPos);
  m_particle_gun->SetParticleMomentum(LPf2);
  m_particle_gun->GeneratePrimaryVertex(anEvent);

  auto Proton = G4Proton::Definition();
  m_particle_gun->SetParticleDefinition(Proton);
  m_particle_gun->SetParticlePosition(LPos);
  m_particle_gun->SetParticleMomentum(LPf3);
  m_particle_gun->GeneratePrimaryVertex(anEvent);
  // double mm_d = (Lv_beam + Lv_targ_D + (-1.)*Lv_K).mag();
  // double mm_p = (Lv_beam + Lv_targ_P + (-1.)*Lv_K).mag();
  // double cos_theta_scat = (beam_mom.x()*LPf1.x()+beam_mom.y()*LPf1.y()+beam_mom.z()*LPf1.z())/(beam_mom.mag()*LPf1.mag());
  // double theta_scat = acos(cos_theta_scat)*(180./acos(-1.));
  // std::cout<<"theta_scat="<<theta_scat<<std::endl;
  // gAnaMan.SetPrimaryInfo(mm_d, mm_p, thetaK, theta_scat, theta_CM);
  // gAnaMan.SetPrimaryParticle(0,LPf1.x(),LPf1.y(),LPf1.z(),KaonPlus->GetPDGMass()/CLHEP::GeV);
  // gAnaMan.SetPrimaryParticle(1,LPf2.x(),LPf2.y(),LPf2.z(),Lambda->GetPDGMass()/CLHEP::GeV);
  // gAnaMan.SetPrimaryParticle(2,LPf3.x(),LPf3.y(),LPf3.z(),Proton->GetPDGMass()/CLHEP::GeV);
  // gAnaMan.SetPrimaryVertex(0,LPos.x(),LPos.y(),LPos.z());
  // gAnaMan.SetPrimaryVertex(1,LPos.x(),LPos.y(),LPos.z());
  // gAnaMan.SetPrimaryVertex(2,LPos.x(),LPos.y(),LPos.z());
}

//_____________________________________________________________________________
// reactio No #2704 pi+ d -> K+ K-pp, K-pp -> SigmaZ p reaction
void
PrimaryGeneratorAction::GenerateE27KppFSigmaZP(G4Event* anEvent)
{
  G4double Mi1=G4PionPlus::Definition()->GetPDGMass();
  G4double Mi2=G4Deuteron::Definition()->GetPDGMass();
  G4double Mf1=G4KaonPlus::Definition()->GetPDGMass();
  G4double Mf2=G4SigmaZero::Definition()->GetPDGMass();
  G4double Mf3=G4Proton::Definition()->GetPDGMass();

  G4double Mp=G4Proton::Definition()->GetPDGMass();
  // G4double Mpim=G4PionMinus::Definition()->GetPDGMass();
  // G4double Mn=G4Neutron::Definition()->GetPDGMass();
  // G4double Mpiz=G4PionZero::Definition()->GetPDGMass();

  G4ThreeVector LPos = GaussPosition_LqTarg(gConf.Get<G4double>("BeamX0"),
                                            gConf.Get<G4double>("BeamY0"),
                                            gGeom.GetGlobalPosition("Target").z(),
                                            gConf.Get<G4double>("BeamDX"),
                                            gConf.Get<G4double>("BeamDY"),
                                            gSize.Get("Target", ThreeVector::X),
                                            gSize.Get("Target", ThreeVector::Z));
  //Note!! env_target_width = Target_Size_z (height of target)

  G4ThreeVector LBeamDir =  GaussDirectionInUV(gConf.Get<G4double>("BeamU0"),
                                               gConf.Get<G4double>("BeamV0"),
                                               gConf.Get<G4double>("BeamDU"),
                                               gConf.Get<G4double>("BeamDV"));

  G4double pb = gConf.Get<G4double>("BeamMom")*CLHEP::GeV;
  G4double dpb = 0.;
  if(gConf.Get<G4double>("BeamMom")!=0.)
    dpb = G4RandGauss::shoot(0.,gConf.Get<G4double>("BeamWidth"))*CLHEP::GeV;
  pb += dpb;

  G4double Mm1 = BreitWigner(2.275, 0.162)*CLHEP::GeV; //E27

  G4double thetaK, theta_CM;
  G4ThreeVector LPm1, LPf1, LPf2, LPf3;

  AGUniform gen1;
  //AGLambda1405 gen1;
  AGUniform gen2;
  G4ThreeVector LPini2=G4ThreeVector(0.,0.,0.);

  bool status = false;
  bool status2 = false;

  G4int n=0;
  while(1){
    if(++n>MaxTry){
      G4Exception(FUNC_NAME,
                  "Production under threshold",
                  RunMustBeAborted,
                  "Production under Threshold!!");
    }

    status=Scattering2Body_theta(Mi1, Mi2, Mf1, Mm1,
                                 pb*LBeamDir,LPini2,
                                 LPf1, LPm1,theta_CM, gen1);
    theta_CM = theta_CM*(180./(acos(-1.)));

    if(status ==true){
      thetaK = LPf1.theta()*(180./(acos(-1.)));
      G4cout<<"thetaK= " <<thetaK <<G4endl;
      if(thetaK<20.){

        status2=Decay2Body(Mm1, Mf2, Mf3, LPm1, LPf2, LPf3, gen2);
        if(status2 == true)
          break;
        else
          std::cout<<"Mm1="<<Mm1<<std::endl;
      }
    }
    pb = gConf.Get<G4double>("BeamMom")*CLHEP::GeV;
    dpb = 0.;
    if(gConf.Get<G4double>("BeamMom")!=0.)
      dpb = G4RandGauss::shoot(0.,gConf.Get<G4double>("BeamWidth"))*CLHEP::GeV;
    pb += dpb;
    Mm1 = BreitWigner(2.275, 0.162)*CLHEP::GeV; //E27
  }

  G4ThreeVector beam_mom = pb*LBeamDir;
  G4LorentzVector Lv_beam, Lv_targ_D, Lv_targ_P, Lv_K;
  Lv_beam.setVectM(pb*LBeamDir, Mi1);
  Lv_targ_D.setVectM(G4ThreeVector(0.,0.,0.), Mi2);
  Lv_targ_P.setVectM(G4ThreeVector(0.,0.,0.), Mp);
  Lv_K.setVectM(LPf1, Mf1);

  // gAnaMan.SetPrimaryBeam(beam_mom.x(),beam_mom.y(),beam_mom.z());
  auto KaonPlus = G4KaonPlus::Definition();
  m_particle_gun->SetParticleDefinition(KaonPlus);
  m_particle_gun->SetParticlePosition(LPos);
  m_particle_gun->SetParticleMomentum(LPf1);
  m_particle_gun->GeneratePrimaryVertex(anEvent);

  auto SigmaZ= G4SigmaZero::Definition();
  m_particle_gun->SetParticleDefinition(SigmaZ);
  m_particle_gun->SetParticlePosition(LPos);
  m_particle_gun->SetParticleMomentum(LPf2);
  m_particle_gun->GeneratePrimaryVertex(anEvent);

  auto Proton= G4Proton::Definition();
  m_particle_gun->SetParticleDefinition(Proton);
  m_particle_gun->SetParticlePosition(LPos);
  m_particle_gun->SetParticleMomentum(LPf3);
  m_particle_gun->GeneratePrimaryVertex(anEvent);

  // double mm_d = (Lv_beam + Lv_targ_D + (-1.)*Lv_K).mag();
  // double mm_p = (Lv_beam + Lv_targ_P + (-1.)*Lv_K).mag();
  // double cos_theta_scat = (beam_mom.x()*LPf1.x()+beam_mom.y()*LPf1.y()+beam_mom.z()*LPf1.z())/(beam_mom.mag()*LPf1.mag());
  // double theta_scat = acos(cos_theta_scat)*(180./acos(-1.));
  // gAnaMan.SetPrimaryInfo(mm_d, mm_p, thetaK, theta_scat, theta_CM);
  // gAnaMan.SetPrimaryParticle(0,LPf1.x(),LPf1.y(),LPf1.z(),KaonPlus->GetPDGMass()/CLHEP::GeV);
  // gAnaMan.SetPrimaryParticle(1,LPf2.x(),LPf2.y(),LPf2.z(),SigmaZ->GetPDGMass()/CLHEP::GeV);
  // gAnaMan.SetPrimaryParticle(2,LPf3.x(),LPf3.y(),LPf3.z(),Proton->GetPDGMass()/CLHEP::GeV);
  // gAnaMan.SetPrimaryVertex(0,LPos.x(),LPos.y(),LPos.z());
  // gAnaMan.SetPrimaryVertex(1,LPos.x(),LPos.y(),LPos.z());
  // gAnaMan.SetPrimaryVertex(2,LPos.x(),LPos.y(),LPos.z());
}

//_____________________________________________________________________________
// reactio No #2705 pi+ d -> K+ K-pp, K-pp -> Lambda piz p reaction
void
PrimaryGeneratorAction::GenerateE27KppFLambdaPizP(G4Event* anEvent)
{
  G4double Mi1=G4PionPlus::Definition()->GetPDGMass();
  G4double Mi2=G4Deuteron::Definition()->GetPDGMass();
  G4double Mf1=G4KaonPlus::Definition()->GetPDGMass();
  G4double Mf2=G4Lambda::Definition()->GetPDGMass();
  G4double Mf3=G4Proton::Definition()->GetPDGMass();
  G4double Mf4=G4PionZero::Definition()->GetPDGMass();

  G4double Mp=G4Proton::Definition()->GetPDGMass();
  // G4double Mpim=G4PionMinus::Definition()->GetPDGMass();
  // G4double Mn=G4Neutron::Definition()->GetPDGMass();
  // G4double Mpiz=G4PionZero::Definition()->GetPDGMass();

  G4ThreeVector LPos = GaussPosition_LqTarg(gConf.Get<G4double>("BeamX0"),
                                            gConf.Get<G4double>("BeamY0"),
                                            gGeom.GetGlobalPosition("Target").z(),
                                            gConf.Get<G4double>("BeamDX"),
                                            gConf.Get<G4double>("BeamDY"),
                                            gSize.Get("Target", ThreeVector::X),
                                            gSize.Get("Target", ThreeVector::Z));
  //Note!! env_target_width = Target_Size_z (height of target)

  G4ThreeVector LBeamDir =  GaussDirectionInUV(gConf.Get<G4double>("BeamU0"),
                                               gConf.Get<G4double>("BeamV0"),
                                               gConf.Get<G4double>("BeamDU"),
                                               gConf.Get<G4double>("BeamDV"));

  G4double pb = gConf.Get<G4double>("BeamMom")*CLHEP::GeV;
  G4double dpb = 0.;
  if(gConf.Get<G4double>("BeamMom")!=0.)
    dpb = G4RandGauss::shoot(0.,gConf.Get<G4double>("BeamWidth"))*CLHEP::GeV;
  pb += dpb;

  G4double Mm1 = BreitWigner(2.275, 0.162)*CLHEP::GeV; //E27

  G4double thetaK, theta_CM;
  G4ThreeVector LPm1, LPf1, LPf2, LPf3, LPf4;

  AGUniform gen1;
  //AGLambda1405 gen1;
  AGUniform gen2;
  G4ThreeVector LPini2=G4ThreeVector(0.,0.,0.);

  bool status = false;
  bool status2 = false;

  G4int n=0;
  while(1){
    if(++n>MaxTry){
      G4Exception(FUNC_NAME,
                  "Production under threshold",
                  RunMustBeAborted,
                  "Production under Threshold!!");
    }

    status=Scattering2Body_theta(Mi1, Mi2, Mf1, Mm1,
                                 pb*LBeamDir,LPini2,
                                 LPf1, LPm1,theta_CM, gen1);
    theta_CM = theta_CM*(180./(acos(-1.)));

    if(status ==true){
      thetaK = LPf1.theta()*(180./(acos(-1.)));
      G4cout<<"thetaK= " <<thetaK <<G4endl;
      if(thetaK<20.){

        status2=Decay3BodyPhaseSpace(Mm1, Mf2, Mf3, Mf4, LPm1, LPf2, LPf3, LPf4);
        if(status2 == true)
          break;
        else
          std::cout<<"Mm1="<<Mm1<<std::endl;
      }
    }
    pb = gConf.Get<G4double>("BeamMom")*CLHEP::GeV;
    dpb = 0.;
    if(gConf.Get<G4double>("BeamMom")!=0.)
      dpb = G4RandGauss::shoot(0.,gConf.Get<G4double>("BeamWidth"))*CLHEP::GeV;
    pb += dpb;
    Mm1 = BreitWigner(2.275, 0.162)*CLHEP::GeV; //E27
  }

  G4ThreeVector beam_mom = pb*LBeamDir;
  G4LorentzVector Lv_beam, Lv_targ_D, Lv_targ_P, Lv_K;
  Lv_beam.setVectM(pb*LBeamDir, Mi1);
  Lv_targ_D.setVectM(G4ThreeVector(0.,0.,0.), Mi2);
  Lv_targ_P.setVectM(G4ThreeVector(0.,0.,0.), Mp);
  Lv_K.setVectM(LPf1, Mf1);

  // gAnaMan.SetPrimaryBeam(beam_mom.x(),beam_mom.y(),beam_mom.z());
  auto KaonPlus = G4KaonPlus::Definition();
  m_particle_gun->SetParticleDefinition(KaonPlus);
  m_particle_gun->SetParticlePosition(LPos);
  m_particle_gun->SetParticleMomentum(LPf1);
  m_particle_gun->GeneratePrimaryVertex(anEvent);

  auto Lambda = G4Lambda::Definition();
  m_particle_gun->SetParticleDefinition(Lambda);
  m_particle_gun->SetParticlePosition(LPos);
  m_particle_gun->SetParticleMomentum(LPf2);
  m_particle_gun->GeneratePrimaryVertex(anEvent);

  auto Proton = G4Proton::Definition();
  m_particle_gun->SetParticleDefinition(Proton);
  m_particle_gun->SetParticlePosition(LPos);
  m_particle_gun->SetParticleMomentum(LPf3);
  m_particle_gun->GeneratePrimaryVertex(anEvent);

  auto PiZero = G4PionZero::Definition();
  m_particle_gun->SetParticleDefinition(PiZero);
  m_particle_gun->SetParticlePosition(LPos);
  m_particle_gun->SetParticleMomentum(LPf4);
  m_particle_gun->GeneratePrimaryVertex(anEvent);

  // double mm_d = (Lv_beam + Lv_targ_D + (-1.)*Lv_K).mag();
  // double mm_p = (Lv_beam + Lv_targ_P + (-1.)*Lv_K).mag();
  // double cos_theta_scat = (beam_mom.x()*LPf1.x()+beam_mom.y()*LPf1.y()+beam_mom.z()*LPf1.z())/(beam_mom.mag()*LPf1.mag());
  // double theta_scat = acos(cos_theta_scat)*(180./acos(-1.));
  // gAnaMan.SetPrimaryInfo(mm_d, mm_p, thetaK, theta_scat, theta_CM);
  // gAnaMan.SetPrimaryParticle(0,LPf1.x(),LPf1.y(),LPf1.z(),KaonPlus->GetPDGMass()/CLHEP::GeV);
  // gAnaMan.SetPrimaryParticle(1,LPf2.x(),LPf2.y(),LPf2.z(),Lambda->GetPDGMass()/CLHEP::GeV);
  // gAnaMan.SetPrimaryParticle(2,LPf3.x(),LPf3.y(),LPf3.z(),Proton->GetPDGMass()/CLHEP::GeV);
  // gAnaMan.SetPrimaryParticle(3,LPf4.x(),LPf4.y(),LPf4.z(),PiZero->GetPDGMass()/CLHEP::GeV);
  // gAnaMan.SetPrimaryVertex(0,LPos.x(),LPos.y(),LPos.z());
  // gAnaMan.SetPrimaryVertex(1,LPos.x(),LPos.y(),LPos.z());
  // gAnaMan.SetPrimaryVertex(2,LPos.x(),LPos.y(),LPos.z());
  // gAnaMan.SetPrimaryVertex(3,LPos.x(),LPos.y(),LPos.z());
}

//_____________________________________________________________________________
// reactio No #2706 pi+ d -> K+ K-pp, K-pp -> SigmaZ piz p reaction
void
PrimaryGeneratorAction::GenerateE27KppFSigmaZPizP(G4Event* anEvent)
{
  G4double Mi1=G4PionPlus::Definition()->GetPDGMass();
  G4double Mi2=G4Deuteron::Definition()->GetPDGMass();
  G4double Mf1=G4KaonPlus::Definition()->GetPDGMass();
  G4double Mf2=G4SigmaZero::Definition()->GetPDGMass();
  G4double Mf3=G4Proton::Definition()->GetPDGMass();
  G4double Mf4=G4PionZero::Definition()->GetPDGMass();

  G4double Mp=G4Proton::Definition()->GetPDGMass();
  // G4double Mpim=G4PionMinus::Definition()->GetPDGMass();
  // G4double Mn=G4Neutron::Definition()->GetPDGMass();
  // G4double Mpiz=G4PionZero::Definition()->GetPDGMass();

  G4ThreeVector LPos = GaussPosition_LqTarg(gConf.Get<G4double>("BeamX0"),
                                            gConf.Get<G4double>("BeamY0"),
                                            gGeom.GetGlobalPosition("Target").z(),
                                            gConf.Get<G4double>("BeamDX"),
                                            gConf.Get<G4double>("BeamDY"),
                                            gSize.Get("Target", ThreeVector::X),
                                            gSize.Get("Target", ThreeVector::Z));
  //Note!! env_target_width = Target_Size_z (height of target)

  G4ThreeVector LBeamDir =  GaussDirectionInUV(gConf.Get<G4double>("BeamU0"),
                                               gConf.Get<G4double>("BeamV0"),
                                               gConf.Get<G4double>("BeamDU"),
                                               gConf.Get<G4double>("BeamDV"));

  G4double pb = gConf.Get<G4double>("BeamMom")*CLHEP::GeV;
  G4double dpb = 0.;
  if(gConf.Get<G4double>("BeamMom")!=0.)
    dpb = G4RandGauss::shoot(0.,gConf.Get<G4double>("BeamWidth"))*CLHEP::GeV;
  pb += dpb;

  G4double Mm1 = BreitWigner(2.275, 0.162)*CLHEP::GeV; //E27

  G4double thetaK, theta_CM;
  G4ThreeVector LPm1, LPf1, LPf2, LPf3, LPf4;

  AGUniform gen1;
  //AGLambda1405 gen1;
  AGUniform gen2;
  G4ThreeVector LPini2=G4ThreeVector(0.,0.,0.);

  bool status = false;
  bool status2 = false;

  G4int n=0;
  while(1){
    if(++n>MaxTry){
      G4Exception(FUNC_NAME,
                  "Production under threshold",
                  RunMustBeAborted,
                  "Production under Threshold!!");
    }

    status=Scattering2Body_theta(Mi1, Mi2, Mf1, Mm1,
                                 pb*LBeamDir,LPini2,
                                 LPf1, LPm1,theta_CM, gen1);
    theta_CM = theta_CM*(180./(acos(-1.)));

    if(status ==true){
      thetaK = LPf1.theta()*(180./(acos(-1.)));
      G4cout<<"thetaK= " <<thetaK <<G4endl;
      if(thetaK<20.){

        status2=Decay3BodyPhaseSpace(Mm1, Mf2, Mf3, Mf4, LPm1, LPf2, LPf3, LPf4);
        if(status2 == true)
          break;
        else
          std::cout<<"Mm1="<<Mm1<<std::endl;
      }
    }
    pb = gConf.Get<G4double>("BeamMom")*CLHEP::GeV;
    dpb = 0.;
    if(gConf.Get<G4double>("BeamMom")!=0.)
      dpb = G4RandGauss::shoot(0.,gConf.Get<G4double>("BeamWidth"))*CLHEP::GeV;
    pb += dpb;
    Mm1 = BreitWigner(2.275, 0.162)*CLHEP::GeV; //E27
  }

  G4ThreeVector beam_mom = pb*LBeamDir;
  G4LorentzVector Lv_beam, Lv_targ_D, Lv_targ_P, Lv_K;
  Lv_beam.setVectM(pb*LBeamDir, Mi1);
  Lv_targ_D.setVectM(G4ThreeVector(0.,0.,0.), Mi2);
  Lv_targ_P.setVectM(G4ThreeVector(0.,0.,0.), Mp);
  Lv_K.setVectM(LPf1, Mf1);

  // gAnaMan.SetPrimaryBeam(beam_mom.x(),beam_mom.y(),beam_mom.z());
  auto KaonPlus = G4KaonPlus::Definition();
  m_particle_gun->SetParticleDefinition(KaonPlus);
  m_particle_gun->SetParticlePosition(LPos);
  m_particle_gun->SetParticleMomentum(LPf1);
  m_particle_gun->GeneratePrimaryVertex(anEvent);

  auto SigmaZero= G4SigmaZero::Definition();
  m_particle_gun->SetParticleDefinition(SigmaZero);
  m_particle_gun->SetParticlePosition(LPos);
  m_particle_gun->SetParticleMomentum(LPf2);
  m_particle_gun->GeneratePrimaryVertex(anEvent);

  auto Proton = G4Proton::Definition();
  m_particle_gun->SetParticleDefinition(Proton);
  m_particle_gun->SetParticlePosition(LPos);
  m_particle_gun->SetParticleMomentum(LPf3);
  m_particle_gun->GeneratePrimaryVertex(anEvent);

  auto PiZero = G4PionZero::Definition();
  m_particle_gun->SetParticleDefinition(PiZero);
  m_particle_gun->SetParticlePosition(LPos);
  m_particle_gun->SetParticleMomentum(LPf4);
  m_particle_gun->GeneratePrimaryVertex(anEvent);

  // double mm_d = (Lv_beam + Lv_targ_D + (-1.)*Lv_K).mag();
  // double mm_p = (Lv_beam + Lv_targ_P + (-1.)*Lv_K).mag();
  // double cos_theta_scat = (beam_mom.x()*LPf1.x()+beam_mom.y()*LPf1.y()+beam_mom.z()*LPf1.z())/(beam_mom.mag()*LPf1.mag());
  // double theta_scat = acos(cos_theta_scat)*(180./acos(-1.));
  // gAnaMan.SetPrimaryInfo(mm_d, mm_p, thetaK, theta_scat, theta_CM);
  // gAnaMan.SetPrimaryParticle(0,LPf1.x(),LPf1.y(),LPf1.z(),KaonPlus->GetPDGMass()/CLHEP::GeV);
  // gAnaMan.SetPrimaryParticle(1,LPf2.x(),LPf2.y(),LPf2.z(),SigmaZero->GetPDGMass()/CLHEP::GeV);
  // gAnaMan.SetPrimaryParticle(2,LPf3.x(),LPf3.y(),LPf3.z(),Proton->GetPDGMass()/CLHEP::GeV);
  // gAnaMan.SetPrimaryParticle(3,LPf4.x(),LPf4.y(),LPf4.z(),PiZero->GetPDGMass()/CLHEP::GeV);
  // gAnaMan.SetPrimaryVertex(0,LPos.x(),LPos.y(),LPos.z());
  // gAnaMan.SetPrimaryVertex(1,LPos.x(),LPos.y(),LPos.z());
  // gAnaMan.SetPrimaryVertex(2,LPos.x(),LPos.y(),LPos.z());
  // gAnaMan.SetPrimaryVertex(3,LPos.x(),LPos.y(),LPos.z());
}

//_____________________________________________________________________________
// reactio No #2707 pi+ d -> K+ K-pp, K-pp -> SigmaP pim p reaction
void
PrimaryGeneratorAction::GenerateE27KppFSigmaPPimP(G4Event* anEvent)
{
  G4double Mi1=G4PionPlus::Definition()->GetPDGMass();
  G4double Mi2=G4Deuteron::Definition()->GetPDGMass();
  G4double Mf1=G4KaonPlus::Definition()->GetPDGMass();
  G4double Mf2=G4SigmaPlus::Definition()->GetPDGMass();
  G4double Mf3=G4Proton::Definition()->GetPDGMass();
  G4double Mf4=G4PionMinus::Definition()->GetPDGMass();

  G4double Mp=G4Proton::Definition()->GetPDGMass();
  // G4double Mpim=G4PionMinus::Definition()->GetPDGMass();
  // G4double Mn=G4Neutron::Definition()->GetPDGMass();
  // G4double Mpiz=G4PionZero::Definition()->GetPDGMass();

  G4ThreeVector LPos = GaussPosition_LqTarg(gConf.Get<G4double>("BeamX0"),
                                            gConf.Get<G4double>("BeamY0"),
                                            gGeom.GetGlobalPosition("Target").z(),
                                            gConf.Get<G4double>("BeamDX"),
                                            gConf.Get<G4double>("BeamDY"),
                                            gSize.Get("Target", ThreeVector::X),
                                            gSize.Get("Target", ThreeVector::Z));
  //Note!! env_target_width = Target_Size_z (height of target)

  G4ThreeVector LBeamDir =  GaussDirectionInUV(gConf.Get<G4double>("BeamU0"),
                                               gConf.Get<G4double>("BeamV0"),
                                               gConf.Get<G4double>("BeamDU"),
                                               gConf.Get<G4double>("BeamDV"));

  G4double pb = gConf.Get<G4double>("BeamMom")*CLHEP::GeV;
  G4double dpb = 0.;
  if(gConf.Get<G4double>("BeamMom")!=0.)
    dpb = G4RandGauss::shoot(0.,gConf.Get<G4double>("BeamWidth"))*CLHEP::GeV;
  pb += dpb;

  G4double Mm1 = BreitWigner(2.275, 0.162)*CLHEP::GeV; //E27

  G4double thetaK, theta_CM;
  G4ThreeVector LPm1, LPf1, LPf2, LPf3, LPf4;

  AGUniform gen1;
  //AGLambda1405 gen1;
  AGUniform gen2;
  G4ThreeVector LPini2=G4ThreeVector(0.,0.,0.);

  bool status = false;
  bool status2 = false;

  G4int n=0;
  while(1){
    if(++n>MaxTry){
      G4Exception(FUNC_NAME,
                  "Production under threshold",
                  RunMustBeAborted,
                  "Production under Threshold!!");
    }

    status=Scattering2Body_theta(Mi1, Mi2, Mf1, Mm1,
                                 pb*LBeamDir,LPini2,
                                 LPf1, LPm1,theta_CM, gen1);
    theta_CM = theta_CM*(180./(acos(-1.)));

    if(status ==true){
      thetaK = LPf1.theta()*(180./(acos(-1.)));
      G4cout<<"thetaK= " <<thetaK <<G4endl;
      if(thetaK<20.){

        status2=Decay3BodyPhaseSpace(Mm1, Mf2, Mf3, Mf4, LPm1, LPf2, LPf3, LPf4);
        if(status2 == true)
          break;
        else
          std::cout<<"Mm1="<<Mm1<<std::endl;
      }
    }
    pb = gConf.Get<G4double>("BeamMom")*CLHEP::GeV;
    dpb = 0.;
    if(gConf.Get<G4double>("BeamMom")!=0.)
      dpb = G4RandGauss::shoot(0.,gConf.Get<G4double>("BeamWidth"))*CLHEP::GeV;
    pb += dpb;
    Mm1 = BreitWigner(2.275, 0.162)*CLHEP::GeV; //E27
  }

  G4ThreeVector beam_mom = pb*LBeamDir;
  G4LorentzVector Lv_beam, Lv_targ_D, Lv_targ_P, Lv_K;
  Lv_beam.setVectM(pb*LBeamDir, Mi1);
  Lv_targ_D.setVectM(G4ThreeVector(0.,0.,0.), Mi2);
  Lv_targ_P.setVectM(G4ThreeVector(0.,0.,0.), Mp);
  Lv_K.setVectM(LPf1, Mf1);

  // gAnaMan.SetPrimaryBeam(beam_mom.x(),beam_mom.y(),beam_mom.z());
  auto KaonPlus = G4KaonPlus::Definition();
  m_particle_gun->SetParticleDefinition(KaonPlus);
  m_particle_gun->SetParticlePosition(LPos);
  m_particle_gun->SetParticleMomentum(LPf1);
  m_particle_gun->GeneratePrimaryVertex(anEvent);

  auto SigmaPlus = G4SigmaPlus::Definition();
  m_particle_gun->SetParticleDefinition(SigmaPlus);
  m_particle_gun->SetParticlePosition(LPos);
  m_particle_gun->SetParticleMomentum(LPf2);
  m_particle_gun->GeneratePrimaryVertex(anEvent);

  auto Proton = G4Proton::Definition();
  m_particle_gun->SetParticleDefinition(Proton);
  m_particle_gun->SetParticlePosition(LPos);
  m_particle_gun->SetParticleMomentum(LPf3);
  m_particle_gun->GeneratePrimaryVertex(anEvent);

  auto PiMinus = G4PionMinus::Definition();
  m_particle_gun->SetParticleDefinition(PiMinus);
  m_particle_gun->SetParticlePosition(LPos);
  m_particle_gun->SetParticleMomentum(LPf4);
  m_particle_gun->GeneratePrimaryVertex(anEvent);

  // double mm_d = (Lv_beam + Lv_targ_D + (-1.)*Lv_K).mag();
  // double mm_p = (Lv_beam + Lv_targ_P + (-1.)*Lv_K).mag();
  // double cos_theta_scat = (beam_mom.x()*LPf1.x()+beam_mom.y()*LPf1.y()+beam_mom.z()*LPf1.z())/(beam_mom.mag()*LPf1.mag());
  // double theta_scat = acos(cos_theta_scat)*(180./acos(-1.));
  // gAnaMan.SetPrimaryInfo(mm_d, mm_p, thetaK, theta_scat, theta_CM);
  // gAnaMan.SetPrimaryParticle(0,LPf1.x(),LPf1.y(),LPf1.z(),KaonPlus->GetPDGMass()/CLHEP::GeV);
  // gAnaMan.SetPrimaryParticle(1,LPf2.x(),LPf2.y(),LPf2.z(),SigmaPlus->GetPDGMass()/CLHEP::GeV);
  // gAnaMan.SetPrimaryParticle(2,LPf3.x(),LPf3.y(),LPf3.z(),Proton->GetPDGMass()/CLHEP::GeV);
  // gAnaMan.SetPrimaryParticle(3,LPf4.x(),LPf4.y(),LPf4.z(),PiMinus->GetPDGMass()/CLHEP::GeV);
  // gAnaMan.SetPrimaryVertex(0,LPos.x(),LPos.y(),LPos.z());
  // gAnaMan.SetPrimaryVertex(1,LPos.x(),LPos.y(),LPos.z());
  // gAnaMan.SetPrimaryVertex(2,LPos.x(),LPos.y(),LPos.z());
  // gAnaMan.SetPrimaryVertex(3,LPos.x(),LPos.y(),LPos.z());
}

//_____________________________________________________________________________
// reactio No #2708 K- 12C -> p 11KB, 11KB -> Lambda 10Be reaction
void
PrimaryGeneratorAction::GenerateE27K11BLambda10Be(G4Event* anEvent)
{
  G4double Mi1=G4KaonPlus::Definition()->GetPDGMass();
  G4double Mi2=12.*AtomicMassUnit*CLHEP::GeV;//12C
  G4double Mf1=G4Proton::Definition()->GetPDGMass();
  G4double Mf2=G4Lambda::Definition()->GetPDGMass();
  G4double Mf3=10.0135338*AtomicMassUnit*CLHEP::GeV;//10Be

  G4double Mp=G4Proton::Definition()->GetPDGMass();
  // G4double Mpim=G4PionMinus::Definition()->GetPDGMass();
  // G4double Mn=G4Neutron::Definition()->GetPDGMass();
  // G4double Mpiz=G4PionZero::Definition()->GetPDGMass();

  G4ThreeVector LPos = GaussPosition_LqTarg(gConf.Get<G4double>("BeamX0"),
                                            gConf.Get<G4double>("BeamY0"),
                                            gGeom.GetGlobalPosition("Target").z(),
                                            gConf.Get<G4double>("BeamDX"),
                                            gConf.Get<G4double>("BeamDY"),
                                            gSize.Get("Target", ThreeVector::X),
                                            gSize.Get("Target", ThreeVector::Z));
  //Note!! env_target_width = Target_Size_z (height of target)

  G4ThreeVector LBeamDir =  GaussDirectionInUV(gConf.Get<G4double>("BeamU0"),
                                               gConf.Get<G4double>("BeamV0"),
                                               gConf.Get<G4double>("BeamDU"),
                                               gConf.Get<G4double>("BeamDV"));

  G4double pb = gConf.Get<G4double>("BeamMom")*CLHEP::GeV;
  G4double dpb = 0.;
  if(gConf.Get<G4double>("BeamMom")!=0.)
    dpb = G4RandGauss::shoot(0.,gConf.Get<G4double>("BeamWidth"))*CLHEP::GeV;
  pb += dpb;

  double Boron11Mass = 11.0093054 * AtomicMassUnit*CLHEP::GeV;
  double MK = G4KaonPlus::Definition()->GetPDGMass();
  G4double Mm1 = BreitWigner(Boron11Mass + MK - 132.5
			     , 183.); //J. Yamagata et al.,

  G4double thetap, theta_CM;
  G4ThreeVector LPm1, LPf1, LPf2, LPf3;

  AGUniform gen1;
  //AGLambda1405 gen1;
  AGUniform gen2;
  G4ThreeVector LPini2=G4ThreeVector(0.,0.,0.);

  bool status = false;
  bool status2 = false;

  G4int n=0;
  while(1){
    if(++n>MaxTry){
      G4Exception(FUNC_NAME,
                  "Production under threshold",
                  RunMustBeAborted,
                  "Production under Threshold!!");
    }

    status=Scattering2Body_theta(Mi1, Mi2, Mf1, Mm1,
                                 pb*LBeamDir,LPini2,
                                 LPf1, LPm1,theta_CM, gen1);
    theta_CM = theta_CM*(180./(acos(-1.)));

    if(status ==true){
      thetap = LPf1.theta()*(180./(acos(-1.)));
      G4cout<<"theta p= " <<thetap <<G4endl;
      if(thetap<20.){
        status2=Decay2Body(Mm1, Mf2, Mf3, LPm1, LPf2, LPf3, gen2);
        if(status2 == true)
          break;
        else
          std::cout<<"Mm1="<<Mm1<<std::endl;
      }
    }
    pb = gConf.Get<G4double>("BeamMom")*CLHEP::GeV;
    dpb = 0.;


    if(gConf.Get<G4double>("BeamWidth") != 0.){
      dpb = G4RandGauss::shoot(0., gConf.Get<G4double>("BeamWidth"))*CLHEP::GeV ;
    }
    pb += dpb;

    Mm1 = BreitWigner(Boron11Mass + MK - 132.5
                      , 183.); //J. Yamagata et al.,
  }

  G4ThreeVector beam_mom = pb*LBeamDir;
  G4LorentzVector Lv_beam, Lv_targ, Lv_targ_p, Lv_p;
  Lv_beam.setVectM(pb*LBeamDir, Mi1);
  Lv_targ.setVectM(G4ThreeVector(0.,0.,0.), Mi2);
  Lv_targ_p.setVectM(G4ThreeVector(0.,0.,0.), Mp);
  Lv_p.setVectM(LPf1, Mf1);

  // gAnaMan.SetPrimaryBeam(beam_mom.x(),beam_mom.y(),beam_mom.z());
  auto Proton = G4Proton::Definition();
  m_particle_gun->SetParticleDefinition(Proton);
  m_particle_gun->SetParticlePosition(LPos);
  m_particle_gun->SetParticleMomentum(LPf1);
  m_particle_gun->GeneratePrimaryVertex(anEvent);

  auto Lambda = G4Lambda::Definition();
  m_particle_gun->SetParticleDefinition(Lambda);
  m_particle_gun->SetParticlePosition(LPos);
  m_particle_gun->SetParticleMomentum(LPf2);
  m_particle_gun->GeneratePrimaryVertex(anEvent);


  G4IonTable* ionTable = G4IonTable::GetIonTable();
  G4int Z = 4;
  G4int A = 10;
  G4ParticleDefinition *Be10 = ionTable->GetIon(Z, A, 0.);
  if (0 == Be10) {
    G4cerr << " ERROR: Can't create nucleus with Z=" << Z << " A=" << A
	   << G4endl;
    ::exit(1);
  }
  // G4ParticleDefinition* ;
  // Proton= particleTable->FindParticle("proton");
  m_particle_gun->SetParticleDefinition(Be10);
  m_particle_gun->SetParticlePosition(LPos);
  m_particle_gun->SetParticleMomentum(LPf3);
  m_particle_gun->GeneratePrimaryVertex(anEvent);

  // double mm1 = (Lv_beam + Lv_targ + (-1.)*Lv_p).mag();
  // double mm_p = (Lv_beam + Lv_targ_p + (-1.)*Lv_p).mag();
  // double cos_theta_scat = (beam_mom.x()*LPf1.x()+beam_mom.y()*LPf1.y()+beam_mom.z()*LPf1.z())/(beam_mom.mag()*LPf1.mag());
  // double theta_scat = acos(cos_theta_scat)*(180./acos(-1.));
  // std::cout<<"theta_scat="<<theta_scat<<std::endl;
  // gAnaMan.SetPrimaryInfo(mm1, mm_p, thetap, theta_scat, theta_CM);
  // gAnaMan.SetPrimaryParticle(0,LPf1.x(),LPf1.y(),LPf1.z(),Proton->GetPDGMass()/CLHEP::GeV);
  // gAnaMan.SetPrimaryParticle(1,LPf2.x(),LPf2.y(),LPf2.z(),Lambda->GetPDGMass()/CLHEP::GeV);
  // gAnaMan.SetPrimaryParticle(2,LPf3.x(),LPf3.y(),LPf3.z(),Be10->GetPDGMass()/CLHEP::GeV);
  // gAnaMan.SetPrimaryVertex(0,LPos.x(),LPos.y(),LPos.z());
  // gAnaMan.SetPrimaryVertex(1,LPos.x(),LPos.y(),LPos.z());
  // gAnaMan.SetPrimaryVertex(2,LPos.x(),LPos.y(),LPos.z());
}

//_____________________________________________________________________________
// reactio No #2709 K+ gun for test
void
PrimaryGeneratorAction::GenerateE27Kptest2(G4Event* anEvent)
{
  int nev = anEvent->GetEventID();
  G4double pbeam;

  for(int i=0; i<14; ++i){
    int fac = nev%14;
    pbeam = (0.1+0.1*(double)fac);
  }

  // G4double Maxang = 30.;
  // G4double cost = cos(Maxang*G4UniformRand()*acos(-1.)/180.);
  G4double cost = cos(0.*acos(-1.)/180.);
  G4double sint =sqrt(1.0 - cost *cost);
  //G4double phi = 360.*(0.5-G4UniformRand())*degree;
  G4double phi = 0.;


  G4double Energy_Kp,  mom_Kp_x, mom_Kp_y, mom_Kp_z;

  auto KaonPlus = G4KaonPlus::Definition();

  mom_Kp_x=pbeam * sint*cos(phi);
  mom_Kp_y=pbeam * sint*sin(phi);
  mom_Kp_z=pbeam * cost;
  Energy_Kp=sqrt(KaonPlus->GetPDGMass()/CLHEP::GeV*KaonPlus->GetPDGMass()/CLHEP::GeV+pbeam*pbeam);

  double theta = acos(cost)*180./acos(-1.);
  std::cout<<"theta ="<<theta<<", phi(rad)="<<phi<<std::endl;
  std::cout<<"Mom =("<<mom_Kp_x<<", "<<mom_Kp_y<<", "<<mom_Kp_z<<std::endl;


  // gAnaMan.SetPrimaryBeam(mom_Kp_x,mom_Kp_y,mom_Kp_z);

  // G4double Ebeam = sqrt(pbeam*pbeam+KaonPlus->GetPDGMass()/CLHEP::GeV*KaonPlus->GetPDGMass()/CLHEP::GeV);


  //  G4double vtx = CLHEP::RandFlat::shoot(-15.,15.)*CLHEP::mm;
  //  G4double vty = CLHEP::RandFlat::shoot(-5.,5.)*CLHEP::mm;
  //  G4double vtz= CLHEP::RandFlat::shoot(env_target_pos_z-env_target_width/2,env_target_pos_z+env_target_width/2)*CLHEP::mm;

  // G4double vtx = CLHEP::RandGauss::shoot(0,10.)*CLHEP::mm;
  // G4double vty = CLHEP::RandFlat::shoot(0.,3.2)*CLHEP::mm;

  G4double vtx = 0.*CLHEP::mm;
  G4double vty = 0.*CLHEP::mm;
  //  G4double vtz= CLHEP::RandFlat::shoot(gGeom.GetGlobalPosition("Target").z()-gSize.Get("Target", ThreeVector::Z)/2,gGeom.GetGlobalPosition("Target").z()+gSize.Get("Target", ThreeVector::Z)/2)*CLHEP::mm-250.*CLHEP::mm;
  G4double vtz= gGeom.GetGlobalPosition("Target").z();
  std::cout<<"pbeam = "<<pbeam<<std::endl;
  // getchar();
  //scat K+
  m_particle_gun->SetParticleDefinition(KaonPlus);
  m_particle_gun->SetParticleMomentumDirection(G4ThreeVector(mom_Kp_x,mom_Kp_y,mom_Kp_z));
  m_particle_gun->SetParticleEnergy((Energy_Kp - KaonPlus->GetPDGMass()/CLHEP::GeV)*CLHEP::GeV);
  m_particle_gun->SetParticlePosition(G4ThreeVector(vtx,vty,vtz));
  m_particle_gun->GeneratePrimaryVertex(anEvent);
  // gAnaMan.SetPrimaryParticle(0,mom_Kp_x,mom_Kp_y,mom_Kp_z,KaonPlus->GetPDGMass()/CLHEP::GeV);
  // gAnaMan.SetPrimaryVertex(0,vtx,vty,vtz);
}

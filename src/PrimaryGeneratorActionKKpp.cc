// -*- C++ -*-

#include "PrimaryGeneratorAction.hh"

#include <G4Event.hh>
#include <G4IonConstructor.hh>
#include <G4LorentzVector.hh>
#include <G4ParticleDefinition.hh>
#include <G4ParticleGun.hh>
#include <G4ParticleTable.hh>
#include <G4ParticleTypes.hh>
#include <G4UImanager.hh>
#include <Randomize.hh>

#include <TTree.h>

#include "AnaManager.hh"
#include "AngDisGenerator.hh"
#include "ConfMan.hh"
#include "DCGeomMan.hh"
#include "DetSizeMan.hh"
#include "GeneratorHelper.hh"
#include "Kinema3Body.hh"
#include "Kinema3Resonance.hh"
#include "Kinema4Body.hh"
#include "KinemaFermi.hh"
#include "KinemaHResonance.hh"
#include "KinemaHweak.hh"
#include "KinemaHybrid.hh"
#include "KinemaKstar.hh"

namespace
{
auto& gAnaMan = AnaManager::GetInstance();
const int MaxTry=1000;
const auto& gConf = ConfMan::GetInstance();
const auto& gGeom = DCGeomMan::GetInstance();
const auto& gSize = DetSizeMan::GetInstance();
}

//_____________________________________________________________________________
// reaction #3001 K- d -> K0 K-K-pp, K-K-pp -> LL -> p pi p pi reaction
void
PrimaryGeneratorAction::GenerateKKppLL1(G4Event* anEvent)
{
  G4double Mp   = G4Proton::Definition()->GetPDGMass();
  G4double Mpim = G4PionMinus::Definition()->GetPDGMass();
  G4double Mpip = G4PionPlus::Definition()->GetPDGMass();
  // G4double Mn   = G4Neutron::Definition()->GetPDGMass();
  // G4double Mpiz = G4PionZero::Definition()->GetPDGMass();
  G4double ML   = G4Lambda::Definition()->GetPDGMass();
  G4double MKz  = G4KaonZero::Definition()->GetPDGMass();
  G4double Mi1  = G4KaonMinus::Definition()->GetPDGMass();
  G4double Mi2  = G4Deuteron::Definition()->GetPDGMass();

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
    dpb = G4RandGauss::shoot(0., gConf.Get<G4double>("BeamWidth"))*CLHEP::GeV;
  pb += dpb;

  double KKppM0 = 1.405 + 1.405; //assuming L(1405)+L(1405)
  double KKppG0 = 0.1; //assumption
  G4double Mm1 = BreitWigner(KKppM0, KKppG0)*CLHEP::GeV; //KKpp

  G4double thetaK, theta_CM;
  G4ThreeVector LPKz, LPm1, LPL1, LPL2, LPf1, LPf2, LPf3, LPf4, LPf5, LPf6;

  AGUniform gen1;
  AGUniform gen2;
  AGUniform gen3;
  AGUniform gen4;

  G4ThreeVector LPini2=G4ThreeVector(0.,0.,0.);

  bool status = false;
  bool status2 = false;
  bool status3 = false;
  bool status4 = false;
  bool status5 = false;


  G4int n=0;
  while(1){
    if(++n>MaxTry){
      G4Exception("PrimaryGeneratorAction::GenerateKKppLL1",
		  "Production under threshold",
		  RunMustBeAborted,
		  "PrimaryGeneratorAction::GenerateProduction under Threshold!!");
    }

    status=Scattering2Body_theta(Mi1, Mi2, MKz, Mm1,
                                 pb*LBeamDir,LPini2,
                                 LPKz, LPm1,theta_CM, gen1);
    theta_CM = theta_CM*(180./(acos(-1.)));

    if(status ==true&&Mm1>0.){
      thetaK = LPKz.theta()*(180./(acos(-1.)));
      G4cout<<"thetaK= " <<thetaK <<G4endl;
      // K0 -> pi+ pi-
      status2=Decay2Body(MKz, Mpip, Mpim, LPKz, LPf1, LPf2, gen2);
      // KKpp -> L L
      status3=Decay2Body(Mm1, ML, ML, LPm1, LPL1, LPL2, gen2);
      if(status2 == true && status3 == true){
        //L -> p pi-
        status4=Decay2Body(ML, Mp, Mpim, LPL1, LPf3, LPf4, gen2);
        //L -> p pi-
        status5=Decay2Body(ML, Mp, Mpim, LPL2, LPf5, LPf6, gen2);

        if(status4 == true && status5 == true)
          break;
      }
      else
        std::cout<<"Mm1="<<Mm1<<std::endl;
    }

    pb = gConf.Get<G4double>("BeamMom")*CLHEP::GeV;
    dpb = 0.;
    if(gConf.Get<G4double>("BeamMom")!=0.)
      dpb = G4RandGauss::shoot(0.,gConf.Get<G4double>("BeamWidth"))*CLHEP::GeV;
    pb += dpb;
    Mm1 = BreitWigner(KKppM0, KKppG0)*CLHEP::GeV; //KKpp
  }

  double BE_KKpp = Mm1 - 2.*Mi1 - 2.*Mp;
  std::cout<<"BE_KKpp="<<BE_KKpp<<", p_K0="<<LPKz.mag()
	   <<", p_L1="<<LPL1.mag()<<", p_L2="<<LPL2.mag()<<std::endl;
  std::cout<<"K0 decay: p_pi+="<<LPf1.mag()<<", p_pi-="<<LPf2.mag()<<std::endl;
  std::cout<<"L decay1: p_p="<<LPf3.mag()<<", p_pi-="<<LPf4.mag()<<std::endl;
  std::cout<<"L decay2: p_p="<<LPf5.mag()<<", p_pi-="<<LPf6.mag()<<std::endl;

  G4ThreeVector beam_mom = pb*LBeamDir;

  G4LorentzVector Lv_beam, Lv_targ_D, Lv_targ_P, Lv_K;
  Lv_beam.setVectM(pb*LBeamDir, Mi1);
  Lv_targ_D.setVectM(G4ThreeVector(0.,0.,0.), Mi2);
  Lv_targ_P.setVectM(G4ThreeVector(0.,0.,0.), Mp);
  Lv_K.setVectM(LPKz, MKz);

  // gAnaMan.SetPrimaryBeam(beam_mom.x(),beam_mom.y(),beam_mom.z());
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* PionPlus;
  PionPlus = particleTable->FindParticle("pi+");

  m_particle_gun->SetParticleDefinition(PionPlus);
  m_particle_gun->SetParticlePosition(LPos);
  m_particle_gun->SetParticleMomentum(LPf1);
  m_particle_gun->GeneratePrimaryVertex(anEvent);


  G4ParticleDefinition* PionMinus;
  PionMinus = particleTable->FindParticle("pi-");
  m_particle_gun->SetParticleDefinition(PionMinus);
  m_particle_gun->SetParticlePosition(LPos);
  m_particle_gun->SetParticleMomentum(LPf2);
  m_particle_gun->GeneratePrimaryVertex(anEvent);

  G4ParticleDefinition* Proton;
  Proton = particleTable->FindParticle("proton");
  m_particle_gun->SetParticleDefinition(Proton);
  m_particle_gun->SetParticlePosition(LPos);
  m_particle_gun->SetParticleMomentum(LPf3);
  m_particle_gun->GeneratePrimaryVertex(anEvent);

  m_particle_gun->SetParticleDefinition(PionMinus);
  m_particle_gun->SetParticlePosition(LPos);
  m_particle_gun->SetParticleMomentum(LPf4);
  m_particle_gun->GeneratePrimaryVertex(anEvent);

  m_particle_gun->SetParticleDefinition(Proton);
  m_particle_gun->SetParticlePosition(LPos);
  m_particle_gun->SetParticleMomentum(LPf5);
  m_particle_gun->GeneratePrimaryVertex(anEvent);

  m_particle_gun->SetParticleDefinition(PionMinus);
  m_particle_gun->SetParticlePosition(LPos);
  m_particle_gun->SetParticleMomentum(LPf6);
  m_particle_gun->GeneratePrimaryVertex(anEvent);

  // double mm_d = (Lv_beam + Lv_targ_D + (-1.)*Lv_K).mag();
  // double mm_p = (Lv_beam + Lv_targ_P + (-1.)*Lv_K).mag();
  // double cos_theta_scat = (beam_mom.x()*LPKz.x()+beam_mom.y()*LPKz.y()+beam_mom.z()*LPKz.z())/(beam_mom.mag()*LPKz.mag());
  // double theta_scat = acos(cos_theta_scat)*(180./acos(-1.));
  // std::cout<<"theta_scat="<<theta_scat<<std::endl;
  // gAnaMan.SetPrimaryInfo(mm_d, mm_p, thetaK, theta_scat, theta_CM);
  // gAnaMan.SetPrimaryParticle(0,LPf1.x(),LPf1.y(),LPf1.z(),
  //                            PionPlus->GetPDGMass()/CLHEP::GeV, PionPlus->GetPDGEncoding());
  // gAnaMan.SetPrimaryParticle(1,LPf2.x(),LPf2.y(),LPf2.z(),
  //                            PionMinus->GetPDGMass()/CLHEP::GeV, PionPlus->GetPDGEncoding());
  // gAnaMan.SetPrimaryParticle(2,LPf3.x(),LPf3.y(),LPf3.z(),
  //                            Proton->GetPDGMass()/CLHEP::GeV, Proton->GetPDGEncoding());
  // gAnaMan.SetPrimaryParticle(3,LPf4.x(),LPf4.y(),LPf4.z(),
  //                            PionMinus->GetPDGMass()/CLHEP::GeV, PionMinus->GetPDGEncoding());
  // gAnaMan.SetPrimaryParticle(4,LPf5.x(),LPf5.y(),LPf5.z(),
  //                            Proton->GetPDGMass()/CLHEP::GeV, Proton->GetPDGEncoding());
  // gAnaMan.SetPrimaryParticle(5,LPf6.x(),LPf6.y(),LPf6.z(),
  //                            PionMinus->GetPDGMass()/CLHEP::GeV, PionMinus->GetPDGEncoding());
  // gAnaMan.SetPrimaryVertex(0,LPos.x(),LPos.y(),LPos.z());
  // gAnaMan.SetPrimaryVertex(1,LPos.x(),LPos.y(),LPos.z());
  // gAnaMan.SetPrimaryVertex(2,LPos.x(),LPos.y(),LPos.z());
  // gAnaMan.SetPrimaryVertex(3,LPos.x(),LPos.y(),LPos.z());
  // gAnaMan.SetPrimaryVertex(4,LPos.x(),LPos.y(),LPos.z());
  // gAnaMan.SetPrimaryVertex(5,LPos.x(),LPos.y(),LPos.z());

}

//_____________________________________________________________________________
// reaction #3002 K- d -> K0 K-K-pp, K-K-pp -> LL reaction
//  K0S and Lambda is directory gunned
void
PrimaryGeneratorAction::GenerateKKppLL2(G4Event* anEvent)
{
  G4double Mp   = G4Proton::Definition()->GetPDGMass();
  // G4double Mpim = G4PionMinus::Definition()->GetPDGMass();
  // G4double Mpip = G4PionPlus::Definition()->GetPDGMass();
  // G4double Mn   = G4Neutron::Definition()->GetPDGMass();
  // G4double Mpiz = G4PionZero::Definition()->GetPDGMass();
  G4double ML   = G4Lambda::Definition()->GetPDGMass();
  G4double MKz  = G4KaonZero::Definition()->GetPDGMass();
  G4double Mi1  = G4KaonMinus::Definition()->GetPDGMass();
  G4double Mi2  = G4Deuteron::Definition()->GetPDGMass();

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

  double KKppM0 = 1.405 + 1.405; //assuming L(1405)+L(1405)
  double KKppG0 = 0.1; //assumption
  G4double Mm1= BreitWigner(KKppM0, KKppG0)*CLHEP::GeV; //KKpp


  G4double thetaK, theta_CM;
  G4ThreeVector LPKz, LPm1, LPL1, LPL2, LPf1, LPf2, LPf3, LPf4, LPf5, LPf6;

  AGUniform gen1;
  AGUniform gen2;
  AGUniform gen3;
  AGUniform gen4;

  G4ThreeVector LPini2=G4ThreeVector(0.,0.,0.);

  bool status = false;
  // bool status2 = false;
  bool status3 = false;
  // bool status4 = false;
  // bool status5 = false;

  G4int n=0;
  while(1){
    if(++n>MaxTry){
      G4Exception("PrimaryGeneratorAction::GenerateKKppLL1",
		  "Production under threshold",
		  RunMustBeAborted,
		  "PrimaryGeneratorAction::GenerateProduction under Threshold!!");
    }

    status=Scattering2Body_theta(Mi1, Mi2, MKz, Mm1,
                                 pb*LBeamDir,LPini2,
                                 LPKz, LPm1,theta_CM, gen1);
    theta_CM = theta_CM*(180./(acos(-1.)));

    if(status ==true&&Mm1>0.){
      thetaK = LPKz.theta()*(180./(acos(-1.)));
      G4cout<<"thetaK= " <<thetaK <<G4endl;

      // K0 -> pi+ pi-
      //status2=Decay2Body(MKz, Mpip, Mpim, LPKz, LPf1, LPf2, gen2);
      // KKpp -> L L
      status3=Decay2Body(Mm1, ML, ML, LPm1, LPL1, LPL2, gen2);
      //if(status2 == true && status3 == true){
      if(status3 == true)
        break;

      //L -> p pi-
      //status4=Decay2Body(ML, Mp, Mpim, LPL1, LPf3, LPf4, gen2);
      //L -> p pi-
      //status5=Decay2Body(ML, Mp, Mpim, LPL2, LPf5, LPf6, gen2);

      //if(status4 == true && status5 == true)
      else
        std::cout<<"Mm1="<<Mm1<<std::endl;
    }

    pb = gConf.Get<G4double>("BeamMom")*CLHEP::GeV;
    dpb = 0.;
    if(gConf.Get<G4double>("BeamMom")!=0.)
      dpb = G4RandGauss::shoot(0.,gConf.Get<G4double>("BeamWidth"))*CLHEP::GeV;
    pb += dpb;
    Mm1 = BreitWigner(KKppM0, KKppG0)*CLHEP::GeV; //KKpp
  }

  double BE_KKpp = Mm1 - 2.*Mi1 - 2.*Mp;
  std::cout<<"BE_KKpp="<<BE_KKpp<<", p_K0="<<LPKz.mag()
	   <<", p_L1="<<LPL1.mag()<<", p_L2="<<LPL2.mag()<<std::endl;
  // std::cout<<"K0 decay: p_pi+="<<LPf1.mag()<<", p_pi-="<<LPf2.mag()<<std::endl;
  // std::cout<<"L decay1: p_p="<<LPf3.mag()<<", p_pi-="<<LPf4.mag()<<std::endl;
  // std::cout<<"L decay2: p_p="<<LPf5.mag()<<", p_pi-="<<LPf6.mag()<<std::endl;

  G4ThreeVector beam_mom = pb*LBeamDir;

  G4LorentzVector Lv_beam, Lv_targ_D, Lv_targ_P, Lv_K;
  Lv_beam.setVectM(pb*LBeamDir, Mi1);
  Lv_targ_D.setVectM(G4ThreeVector(0.,0.,0.), Mi2);
  Lv_targ_P.setVectM(G4ThreeVector(0.,0.,0.), Mp);
  Lv_K.setVectM(LPKz, MKz);

  // gAnaMan.SetPrimaryBeam(beam_mom.x(),beam_mom.y(),beam_mom.z());
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* Kaon0S;
  Kaon0S = particleTable->FindParticle("kaon0S");
  m_particle_gun->SetParticleDefinition(Kaon0S);
  m_particle_gun->SetParticlePosition(LPos);
  m_particle_gun->SetParticleMomentum(LPKz);
  m_particle_gun->GeneratePrimaryVertex(anEvent);


  G4ParticleDefinition* Lambda;
  Lambda = particleTable->FindParticle("lambda");
  m_particle_gun->SetParticleDefinition(Lambda);
  m_particle_gun->SetParticlePosition(LPos);
  m_particle_gun->SetParticleMomentum(LPL1);
  m_particle_gun->GeneratePrimaryVertex(anEvent);

  m_particle_gun->SetParticleDefinition(Lambda);
  m_particle_gun->SetParticlePosition(LPos);
  m_particle_gun->SetParticleMomentum(LPL2);
  m_particle_gun->GeneratePrimaryVertex(anEvent);

  // double mm_d = (Lv_beam + Lv_targ_D + (-1.)*Lv_K).mag();
  // double mm_p = (Lv_beam + Lv_targ_P + (-1.)*Lv_K).mag();
  // double cos_theta_scat = (beam_mom.x()*LPKz.x()+beam_mom.y()*LPKz.y()+beam_mom.z()*LPKz.z())/(beam_mom.mag()*LPKz.mag());
  // double theta_scat = acos(cos_theta_scat)*(180./acos(-1.));
  // std::cout<<"theta_scat="<<theta_scat<<std::endl;
  // gAnaMan.SetPrimaryInfo(mm_d, mm_p, thetaK, theta_scat, theta_CM);
  // gAnaMan.SetPrimaryParticle(0,LPKz.x(),LPKz.y(),LPKz.z(),
  //                            Kaon0S->GetPDGMass()/CLHEP::GeV, Kaon0S->GetPDGEncoding());
  // gAnaMan.SetPrimaryParticle(1,LPL1.x(),LPL1.y(),LPL1.z(),
  //                            Lambda->GetPDGMass()/CLHEP::GeV, Lambda->GetPDGEncoding());
  // gAnaMan.SetPrimaryParticle(2,LPL2.x(),LPL2.y(),LPL2.z(),
  //                            Lambda->GetPDGMass()/CLHEP::GeV, Lambda->GetPDGEncoding());

  // gAnaMan.SetPrimaryVertex(0,LPos.x(),LPos.y(),LPos.z());
  // gAnaMan.SetPrimaryVertex(1,LPos.x(),LPos.y(),LPos.z());
  // gAnaMan.SetPrimaryVertex(2,LPos.x(),LPos.y(),LPos.z());

}

//_____________________________________________________________________________
// reaction #3003 K- d -> K0 K-K-pp, K-K-pp -> LS-pi+ reaction
//  K0S, Lambda, Sigma is directory gunned
void
PrimaryGeneratorAction::GenerateKKppLSmPip(G4Event* anEvent)
{
  G4double Mp   = G4Proton::Definition()->GetPDGMass();
  // G4double Mpim = G4PionMinus::Definition()->GetPDGMass();
  G4double Mpip = G4PionPlus::Definition()->GetPDGMass();
  // G4double Mn   = G4Neutron::Definition()->GetPDGMass();
  // G4double Mpiz = G4PionZero::Definition()->GetPDGMass();
  G4double ML   = G4Lambda::Definition()->GetPDGMass();
  G4double MSm  = G4SigmaMinus::Definition()->GetPDGMass();
  // G4double MSp  = G4SigmaPlus::Definition()->GetPDGMass();
  G4double MKz  = G4KaonZero::Definition()->GetPDGMass();
  G4double Mi1  = G4KaonMinus::Definition()->GetPDGMass();
  G4double Mi2  = G4Deuteron::Definition()->GetPDGMass();

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

  double KKppM0 = 1.405 + 1.405; //assuming L(1405)+L(1405)
  double KKppG0 = 0.1; //assumption
  G4double Mm1 = BreitWigner(KKppM0, KKppG0)*CLHEP::GeV; //KKpp

  G4double thetaK, theta_CM;
  G4ThreeVector LPKz, LPm1, LPL, LPS, LPpi;

  AGUniform gen1;
  AGUniform gen2;
  AGUniform gen3;
  AGUniform gen4;

  G4ThreeVector LPini2=G4ThreeVector(0.,0.,0.);

  bool status = false;
  // bool status2 = false;
  bool status3 = false;
  // bool status4 = false;
  // bool status5 = false;


  G4int n=0;
  while(1){
    if(++n>MaxTry){
      G4Exception("PrimaryGeneratorAction::GenerateKKppLSmPip",
		  "Production under threshold",
		  RunMustBeAborted,
		  "PrimaryGeneratorAction::GenerateProduction under Threshold!!");
    }

    status=Scattering2Body_theta(Mi1, Mi2, MKz, Mm1,
                                 pb*LBeamDir,LPini2,
                                 LPKz, LPm1,theta_CM, gen1);
    theta_CM = theta_CM*(180./(acos(-1.)));

    if(status ==true&&Mm1>0.){
      std::cout<<"Mm1="<<Mm1<<", ML+MS+Mpi="<<ML+MSm+Mpip<<std::endl;
      thetaK = LPKz.theta()*(180./(acos(-1.)));
      G4cout<<"thetaK= " <<thetaK <<G4endl;

      // K0 -> pi+ pi-
      //status2=Decay2Body(MKz, Mpip, Mpim, LPKz, LPf1, LPf2, gen2);
      // KKpp -> L Sm pi+
      status3=Decay3BodyPhaseSpace(Mm1, ML, MSm, Mpip, LPm1, LPL, LPS, LPpi);
      //if(status2 == true && status3 == true){
      if(status3 == true)
        break;

      //L -> p pi-
      //status4=Decay2Body(ML, Mp, Mpim, LPL1, LPf3, LPf4, gen2);
      //L -> p pi-
      //status5=Decay2Body(ML, Mp, Mpim, LPL2, LPf5, LPf6, gen2);

      //if(status4 == true && status5 == true)
      else
        std::cout<<"Mm1="<<Mm1<<std::endl;
    }

    pb = gConf.Get<G4double>("BeamMom")*CLHEP::GeV;
    dpb = 0.;
    if(gConf.Get<G4double>("BeamMom")!=0.)
      dpb = G4RandGauss::shoot(0.,gConf.Get<G4double>("BeamWidth"))*CLHEP::GeV;
    pb += dpb;
    Mm1 = BreitWigner(KKppM0, KKppG0)*CLHEP::GeV; //KKpp
  }

  double BE_KKpp = Mm1 - 2.*Mi1 - 2.*Mp;
  std::cout<<"BE_KKpp="<<BE_KKpp<<", p_K0="<<LPKz.mag()
	   <<", p_L="<<LPL.mag()<<", p_S="<<LPS.mag()<<", p_pi="<<LPpi.mag()<<std::endl;
  // std::cout<<"K0 decay: p_pi+="<<LPf1.mag()<<", p_pi-="<<LPf2.mag()<<std::endl;
  // std::cout<<"L decay1: p_p="<<LPf3.mag()<<", p_pi-="<<LPf4.mag()<<std::endl;
  // std::cout<<"L decay2: p_p="<<LPf5.mag()<<", p_pi-="<<LPf6.mag()<<std::endl;

  G4ThreeVector beam_mom = pb*LBeamDir;
  G4LorentzVector Lv_beam, Lv_targ_D, Lv_targ_P, Lv_K;
  Lv_beam.setVectM(pb*LBeamDir, Mi1);
  Lv_targ_D.setVectM(G4ThreeVector(0.,0.,0.), Mi2);
  Lv_targ_P.setVectM(G4ThreeVector(0.,0.,0.), Mp);
  Lv_K.setVectM(LPKz, MKz);

  // gAnaMan.SetPrimaryBeam(beam_mom.x(),beam_mom.y(),beam_mom.z());
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* Kaon0S;
  Kaon0S = particleTable->FindParticle("kaon0S");
  m_particle_gun->SetParticleDefinition(Kaon0S);
  m_particle_gun->SetParticlePosition(LPos);
  m_particle_gun->SetParticleMomentum(LPKz);
  m_particle_gun->GeneratePrimaryVertex(anEvent);


  G4ParticleDefinition* Lambda;
  Lambda = particleTable->FindParticle("lambda");
  m_particle_gun->SetParticleDefinition(Lambda);
  m_particle_gun->SetParticlePosition(LPos);
  m_particle_gun->SetParticleMomentum(LPL);
  m_particle_gun->GeneratePrimaryVertex(anEvent);

  G4ParticleDefinition* Sigma;
  Sigma = particleTable->FindParticle("sigma-");
  m_particle_gun->SetParticleDefinition(Sigma);
  m_particle_gun->SetParticlePosition(LPos);
  m_particle_gun->SetParticleMomentum(LPS);
  m_particle_gun->GeneratePrimaryVertex(anEvent);

  G4ParticleDefinition* Pi;
  Pi = particleTable->FindParticle("pi+");
  m_particle_gun->SetParticleDefinition(Pi);
  m_particle_gun->SetParticlePosition(LPos);
  m_particle_gun->SetParticleMomentum(LPpi);
  m_particle_gun->GeneratePrimaryVertex(anEvent);

  // double mm_d = (Lv_beam + Lv_targ_D + (-1.)*Lv_K).mag();
  // double mm_p = (Lv_beam + Lv_targ_P + (-1.)*Lv_K).mag();
  // double cos_theta_scat = (beam_mom.x()*LPKz.x()+beam_mom.y()*LPKz.y()+beam_mom.z()*LPKz.z())/(beam_mom.mag()*LPKz.mag());
  // double theta_scat = acos(cos_theta_scat)*(180./acos(-1.));
  // std::cout<<"theta_scat="<<theta_scat<<std::endl;
  // gAnaMan.SetPrimaryInfo(mm_d, mm_p, thetaK, theta_scat, theta_CM);
  // gAnaMan.SetPrimaryParticle(0,LPKz.x(),LPKz.y(),LPKz.z(),
  //                            Kaon0S->GetPDGMass()/CLHEP::GeV, Kaon0S->GetPDGEncoding());
  // gAnaMan.SetPrimaryParticle(1,LPL.x(),LPL.y(),LPL.z(),
  //                            Lambda->GetPDGMass()/CLHEP::GeV, Lambda->GetPDGEncoding());
  // gAnaMan.SetPrimaryParticle(2,LPS.x(),LPS.y(),LPS.z(),
  //                            Sigma->GetPDGMass()/CLHEP::GeV, Sigma->GetPDGEncoding());
  // gAnaMan.SetPrimaryParticle(3,LPpi.x(),LPpi.y(),LPpi.z(),
  //                            Pi->GetPDGMass()/CLHEP::GeV, Pi->GetPDGEncoding());

  // gAnaMan.SetPrimaryVertex(0,LPos.x(),LPos.y(),LPos.z());
  // gAnaMan.SetPrimaryVertex(1,LPos.x(),LPos.y(),LPos.z());
  // gAnaMan.SetPrimaryVertex(2,LPos.x(),LPos.y(),LPos.z());
  // gAnaMan.SetPrimaryVertex(3,LPos.x(),LPos.y(),LPos.z());
}

//_____________________________________________________________________________
// reaction #3004 K- d -> K0 K-K-pp, K-K-pp -> LS+pi- reaction
// K0S, Lambda, Sigma is directory gunned
void
PrimaryGeneratorAction::GenerateKKppLSpPim(G4Event* anEvent)
{
  G4double Mp   = G4Proton::Definition()->GetPDGMass();
  G4double Mpim = G4PionMinus::Definition()->GetPDGMass();
  // G4double Mpip = G4PionPlus::Definition()->GetPDGMass();
  // G4double Mn   = G4Neutron::Definition()->GetPDGMass();
  // G4double Mpiz = G4PionZero::Definition()->GetPDGMass();
  G4double ML   = G4Lambda::Definition()->GetPDGMass();
  // G4double MSm  = G4SigmaMinus::Definition()->GetPDGMass();
  G4double MSp  = G4SigmaPlus::Definition()->GetPDGMass();
  G4double MKz  = G4KaonZero::Definition()->GetPDGMass();
  G4double Mi1  = G4KaonMinus::Definition()->GetPDGMass();
  G4double Mi2  = G4Deuteron::Definition()->GetPDGMass();

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

  double KKppM0 = 1.405 + 1.405; //assuming L(1405)+L(1405)
  double KKppG0 = 0.1; //assumption
  G4double Mm1 = BreitWigner(KKppM0, KKppG0)*CLHEP::GeV; //KKpp

  G4double thetaK, theta_CM;
  G4ThreeVector LPKz, LPm1, LPL, LPS, LPpi;

  AGUniform gen1;
  AGUniform gen2;
  AGUniform gen3;
  AGUniform gen4;

  G4ThreeVector LPini2=G4ThreeVector(0.,0.,0.);

  bool status = false;
  // bool status2 = false;
  bool status3 = false;
  // bool status4 = false;
  // bool status5 = false;

  G4int n=0;
  while(1){
    if(++n>MaxTry){
      G4Exception("PrimaryGeneratorAction::GenerateKKppLSpPim",
		  "Production under threshold",
		  RunMustBeAborted,
		  "PrimaryGeneratorAction::GenerateProduction under Threshold!!");
    }

    status=Scattering2Body_theta(Mi1, Mi2, MKz, Mm1,
                                 pb*LBeamDir,LPini2,
                                 LPKz, LPm1,theta_CM, gen1);
    theta_CM = theta_CM*(180./(acos(-1.)));

    if(status ==true&&Mm1>0.){
      thetaK = LPKz.theta()*(180./(acos(-1.)));
      G4cout<<"thetaK= " <<thetaK <<G4endl;

      // K0 -> pi+ pi-
      //status2=Decay2Body(MKz, Mpip, Mpim, LPKz, LPf1, LPf2, gen2);
      // KKpp -> L Sm pi+
      status3=Decay3BodyPhaseSpace(Mm1, ML, MSp, Mpim, LPm1, LPL, LPS, LPpi);
      //if(status2 == true && status3 == true){
      if(status3 == true)
        break;

      //L -> p pi-
      //status4=Decay2Body(ML, Mp, Mpim, LPL1, LPf3, LPf4, gen2);
      //L -> p pi-
      //status5=Decay2Body(ML, Mp, Mpim, LPL2, LPf5, LPf6, gen2);

      //if(status4 == true && status5 == true)
      else
        std::cout<<"Mm1="<<Mm1<<std::endl;
    }

    pb = gConf.Get<G4double>("BeamMom")*CLHEP::GeV;
    dpb = 0.;
    if(gConf.Get<G4double>("BeamMom")!=0.)
      dpb = G4RandGauss::shoot(0.,gConf.Get<G4double>("BeamWidth"))*CLHEP::GeV;
    pb += dpb;
    Mm1 = BreitWigner(KKppM0, KKppG0)*CLHEP::GeV; //KKpp
  }

  double BE_KKpp = Mm1 - 2.*Mi1 - 2.*Mp;
  std::cout<<"BE_KKpp="<<BE_KKpp<<", p_K0="<<LPKz.mag()
	   <<", p_L="<<LPL.mag()<<", p_S="<<LPS.mag()<<", p_pi="<<LPpi.mag()<<std::endl;
  // std::cout<<"K0 decay: p_pi+="<<LPf1.mag()<<", p_pi-="<<LPf2.mag()<<std::endl;
  // std::cout<<"L decay1: p_p="<<LPf3.mag()<<", p_pi-="<<LPf4.mag()<<std::endl;
  // std::cout<<"L decay2: p_p="<<LPf5.mag()<<", p_pi-="<<LPf6.mag()<<std::endl;

  G4ThreeVector beam_mom = pb*LBeamDir;
  G4LorentzVector Lv_beam, Lv_targ_D, Lv_targ_P, Lv_K;
  Lv_beam.setVectM(pb*LBeamDir, Mi1);
  Lv_targ_D.setVectM(G4ThreeVector(0.,0.,0.), Mi2);
  Lv_targ_P.setVectM(G4ThreeVector(0.,0.,0.), Mp);
  Lv_K.setVectM(LPKz, MKz);

  // gAnaMan.SetPrimaryBeam(beam_mom.x(),beam_mom.y(),beam_mom.z());
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* Kaon0S;
  Kaon0S = particleTable->FindParticle("kaon0S");
  m_particle_gun->SetParticleDefinition(Kaon0S);
  m_particle_gun->SetParticlePosition(LPos);
  m_particle_gun->SetParticleMomentum(LPKz);
  m_particle_gun->GeneratePrimaryVertex(anEvent);


  G4ParticleDefinition* Lambda;
  Lambda = particleTable->FindParticle("lambda");
  m_particle_gun->SetParticleDefinition(Lambda);
  m_particle_gun->SetParticlePosition(LPos);
  m_particle_gun->SetParticleMomentum(LPL);
  m_particle_gun->GeneratePrimaryVertex(anEvent);

  G4ParticleDefinition* Sigma;
  Sigma = particleTable->FindParticle("sigma+");
  m_particle_gun->SetParticleDefinition(Sigma);
  m_particle_gun->SetParticlePosition(LPos);
  m_particle_gun->SetParticleMomentum(LPS);
  m_particle_gun->GeneratePrimaryVertex(anEvent);

  G4ParticleDefinition* Pi;
  Pi = particleTable->FindParticle("pi-");
  m_particle_gun->SetParticleDefinition(Pi);
  m_particle_gun->SetParticlePosition(LPos);
  m_particle_gun->SetParticleMomentum(LPpi);
  m_particle_gun->GeneratePrimaryVertex(anEvent);

  // double mm_d = (Lv_beam + Lv_targ_D + (-1.)*Lv_K).mag();
  // double mm_p = (Lv_beam + Lv_targ_P + (-1.)*Lv_K).mag();
  // double cos_theta_scat = (beam_mom.x()*LPKz.x()+beam_mom.y()*LPKz.y()+beam_mom.z()*LPKz.z())/(beam_mom.mag()*LPKz.mag());
  // double theta_scat = acos(cos_theta_scat)*(180./acos(-1.));
  // std::cout<<"theta_scat="<<theta_scat<<std::endl;
  // gAnaMan.SetPrimaryInfo(mm_d, mm_p, thetaK, theta_scat, theta_CM);
  // gAnaMan.SetPrimaryParticle(0,LPKz.x(),LPKz.y(),LPKz.z(),
  //                            Kaon0S->GetPDGMass()/CLHEP::GeV, Kaon0S->GetPDGEncoding());
  // gAnaMan.SetPrimaryParticle(1,LPL.x(),LPL.y(),LPL.z(),
  //                            Lambda->GetPDGMass()/CLHEP::GeV, Lambda->GetPDGEncoding());
  // gAnaMan.SetPrimaryParticle(2,LPS.x(),LPS.y(),LPS.z(),
  //                            Sigma->GetPDGMass()/CLHEP::GeV, Sigma->GetPDGEncoding());
  // gAnaMan.SetPrimaryParticle(3,LPpi.x(),LPpi.y(),LPpi.z(),
  //                            Pi->GetPDGMass()/CLHEP::GeV, Pi->GetPDGEncoding());

  // gAnaMan.SetPrimaryVertex(0,LPos.x(),LPos.y(),LPos.z());
  // gAnaMan.SetPrimaryVertex(1,LPos.x(),LPos.y(),LPos.z());
  // gAnaMan.SetPrimaryVertex(2,LPos.x(),LPos.y(),LPos.z());
  // gAnaMan.SetPrimaryVertex(3,LPos.x(),LPos.y(),LPos.z());
}

int Nbeam_JAM=0;

//_____________________________________________________________________________
// reaction #3101 JAM input
void
PrimaryGeneratorAction::GenerateJAMInput(G4Event* anEvent, TTree*t1)
{
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  // G4ParticleDefinition* kaonMinus = particleTable->FindParticle("kaon-");
  G4double pbeam = gConf.Get<G4double>("BeamMom");
  //  pbeam=1.8;
  // G4double mom_kp_x = 0;
  // G4double mom_kp_y = 0;
  // G4double mom_kp_z = pbeam;
  // gAnaMan.SetPrimaryBeam(0,0,pbeam);

  int N_first = gConf.Get<G4int>("NbeamFirst");

  int Nbeam_JAMInput = Nbeam_JAM+N_first;

  int ev_max= t1->GetEntries();
  int np;
  int pid[20];
  G4double  px[20], py[20], pz[20];
  t1->SetBranchAddress("np", &np);
  t1->SetBranchAddress("pid", pid);
  t1->SetBranchAddress("px", px);
  t1->SetBranchAddress("py", py);
  t1->SetBranchAddress("pz", pz);

  //  t1->GetEntry(Nbeam_JAMInput%ev_max);
  if(Nbeam_JAMInput<=ev_max){

    t1->GetEntry(Nbeam_JAMInput);
    G4int np_JAM = np;

    //  G4int np_JAM = Getnp_JAM(Nbeam_JAMInput);
    // gAnaMan.SetPrimaryInfo(0., 0., 0., 0., 0.);



    std::cout<<"Nbeam_JAMInpu="<<Nbeam_JAMInput
             <<", Num of Primary Particle="<<np_JAM<<std::endl;;
    for(int inp=0; inp<np_JAM; ++inp){
      G4ParticleDefinition *ptmp;
      G4int pid_JAM = pid[inp];
      //G4int pid_JAM = GetPID_JAM(Nbeam_JAMInput, inp);
      if(pid_JAM==2212) ptmp=particleTable->FindParticle("proton");
      else if(pid_JAM==2112) ptmp=particleTable->FindParticle("neutron");
      else if(pid_JAM==3122) ptmp=particleTable->FindParticle("lambda");
      else if(pid_JAM==3112) ptmp=particleTable->FindParticle("sigma-");
      else if(pid_JAM==3212) ptmp=particleTable->FindParticle("sigma0");
      else if(pid_JAM==3222) ptmp=particleTable->FindParticle("sigma+");
      else if(pid_JAM==3312) ptmp=particleTable->FindParticle("xi-");
      else if(pid_JAM==3322) ptmp=particleTable->FindParticle("xi0");
      else if(pid_JAM==-211) ptmp=particleTable->FindParticle("pi-");
      else if(pid_JAM==111) ptmp=particleTable->FindParticle("pi0");
      else if(pid_JAM==211) ptmp=particleTable->FindParticle("pi+");
      else if(pid_JAM==221) ptmp=particleTable->FindParticle("eta");
      else if(pid_JAM==22) ptmp=particleTable->FindParticle("gamma");
      else if(pid_JAM==-321) ptmp=particleTable->FindParticle("kaon-");
      else if(pid_JAM==321) ptmp=particleTable->FindParticle("kaon+");
      else if(pid_JAM==311){
        if(G4UniformRand()>0.5)
          ptmp=particleTable->FindParticle("kaon0S");
        else
          ptmp=particleTable->FindParticle("kaon0L");
      }
      else if(pid_JAM==-311){
        if(G4UniformRand()>0.5)
          ptmp=particleTable->FindParticle("kaon0S");
        else
          ptmp=particleTable->FindParticle("kaon0L");
      }
      else{
        std::cout<<"pid_JAM="<<pid_JAM<<std::endl;
        //getchar();
        continue;
      }
      G4ThreeVector LPos = G4ThreeVector(0.,0.,gGeom.GetGlobalPosition("Target").z());
      //G4ThreeVector LP = GetP_JAM(Nbeam_JAMInput, inp);
      G4ThreeVector LP = G4ThreeVector(px[inp]*CLHEP::GeV, py[inp]*CLHEP::GeV, pz[inp]*CLHEP::GeV);

      m_particle_gun->SetParticleDefinition(ptmp);
      m_particle_gun->SetParticlePosition(LPos);
      m_particle_gun->SetParticleMomentum(LP);
      m_particle_gun->GeneratePrimaryVertex(anEvent);

      //    std::cout<<"encording="<<ptmp->GetPDGEncoding()<<std::endl;
      // if(pid_JAM==-311)
      //   gAnaMan.SetPrimaryParticle(inp,LP.x(),LP.y(),LP.z(),
      //                              -1.*ptmp->GetPDGMass()/CLHEP::GeV, ptmp->GetPDGEncoding());
      // else
      //   gAnaMan.SetPrimaryParticle(inp,LP.x(),LP.y(),LP.z(),
      //                              ptmp->GetPDGMass()/CLHEP::GeV, ptmp->GetPDGEncoding());
      // //gAnaMan.SetPrimaryVertex(0,LPos.x(),LPos.y(),LPos.z());
      // gAnaMan.SetPrimaryVertex(inp,LPos.x(),LPos.y(),LPos.z());

      //delete ptmp;
    }
  }

  ++Nbeam_JAM;
}

//_____________________________________________________________________________
//reactio No #3102 K- beam through gaus without correlation
void
PrimaryGeneratorAction::GenerateKKppBeamThrough1(G4Event* anEvent)
{
  // G4double Mp=G4Proton::Definition()->GetPDGMass();
  // G4double Mpim=G4PionMinus::Definition()->GetPDGMass();
  // G4double Mpip=G4PionPlus::Definition()->GetPDGMass();
  // G4double Mn=G4Neutron::Definition()->GetPDGMass();
  // G4double Mpiz=G4PionZero::Definition()->GetPDGMass();
  // G4double ML=G4Lambda::Definition()->GetPDGMass();
  // G4double MKz=G4KaonZero::Definition()->GetPDGMass();

  // G4double Mi1=G4KaonMinus::Definition()->GetPDGMass();
  // G4double Mi2=G4Deuteron::Definition()->GetPDGMass();


  double beam_x = G4RandGauss::shoot(gConf.Get<G4double>("BeamX0"),gConf.Get<G4double>("BeamDX"));
  double beam_y = G4RandGauss::shoot(gConf.Get<G4double>("BeamY0"),gConf.Get<G4double>("BeamDY"));
  //double beam_z = -129.;
  double beam_z = -300.;
  G4ThreeVector LPos = G4ThreeVector(beam_x, beam_y, beam_z);

  G4ThreeVector LBeamDir =  GaussDirectionInUV(gConf.Get<G4double>("BeamU0"),
                                               gConf.Get<G4double>("BeamV0"),
                                               gConf.Get<G4double>("BeamDU"),
                                               gConf.Get<G4double>("BeamDV"));

  G4double pb = gConf.Get<G4double>("BeamMom")*CLHEP::GeV;
  //G4double pb = 0.3*CLHEP::GeV;
  G4double dpb = 0.;
  if(gConf.Get<G4double>("BeamMom")!=0.)
    dpb = G4RandGauss::shoot(0.,gConf.Get<G4double>("BeamWidth"))*CLHEP::GeV;
  pb += dpb;


  G4ThreeVector beam_mom = pb*LBeamDir;

  // gAnaMan.SetPrimaryBeam(beam_mom.x(),beam_mom.y(),beam_mom.z());
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* KaonM;
  KaonM = particleTable->FindParticle("kaon-");
  m_particle_gun->SetParticleDefinition(KaonM);
  m_particle_gun->SetParticlePosition(LPos);
  m_particle_gun->SetParticleMomentum(beam_mom);
  m_particle_gun->GeneratePrimaryVertex(anEvent);

  // gAnaMan.SetPrimaryInfo(0., 0., 0., 0., 0.);
  // gAnaMan.SetPrimaryParticle(0,beam_mom.x(),beam_mom.y(),beam_mom.z(),
  //                            KaonM->GetPDGMass()/CLHEP::GeV, KaonM->GetPDGEncoding());
  // gAnaMan.SetPrimaryVertex(0,LPos.x(),LPos.y(),LPos.z());
}

//_____________________________________________________________________________
//reaction No #3103 JAM input K0
void
PrimaryGeneratorAction::GenerateJAMInputK0(G4Event* anEvent, TTree*t1)
{
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  // G4ParticleDefinition* kaonMinus = particleTable->FindParticle("kaon-");
  G4double pbeam=gConf.Get<G4double>("BeamMom");
  //  pbeam=1.8;
  // G4double mom_kp_x = 0;
  // G4double mom_kp_y = 0;
  // G4double mom_kp_z = pbeam;
  // gAnaMan.SetPrimaryBeam(0,0,pbeam);

  int N_first = gConf.Get<G4int>("NbeamFirst");

  int Nbeam_JAMInput = Nbeam_JAM+N_first;

  int ev_max= t1->GetEntries();
  int np;
  int pid[20];
  G4double  px[20], py[20], pz[20];
  t1->SetBranchAddress("np", &np);
  t1->SetBranchAddress("pid", pid);
  t1->SetBranchAddress("px", px);
  t1->SetBranchAddress("py", py);
  t1->SetBranchAddress("pz", pz);

  t1->GetEntry(Nbeam_JAMInput%ev_max);
  G4int np_JAM = np;

  //  G4int np_JAM = Getnp_JAM(Nbeam_JAMInput);
  // gAnaMan.SetPrimaryInfo(0., 0., 0., 0., 0.);
  std::cout<<"Nbeam_JAMInpu="<<Nbeam_JAMInput
	   <<", Num of Primary Particle="<<np_JAM<<std::endl;;
  G4int K0flag=0;
  for(int inp=0; inp<np_JAM; ++inp){
    G4int pid_JAM = pid[inp];
    if(pid_JAM==311)
      ++K0flag;
  }

  if(K0flag>0){
    for(int inp=0; inp<np_JAM; ++inp){
      G4ParticleDefinition *ptmp;
      G4int pid_JAM = pid[inp];
      //G4int pid_JAM = GetPID_JAM(Nbeam_JAMInput, inp);
      if(pid_JAM==2212) ptmp=particleTable->FindParticle("proton");
      else if(pid_JAM==2112) ptmp=particleTable->FindParticle("neutron");
      else if(pid_JAM==3122) ptmp=particleTable->FindParticle("lambda");
      else if(pid_JAM==3112) ptmp=particleTable->FindParticle("sigma-");
      else if(pid_JAM==3212) ptmp=particleTable->FindParticle("sigma0");
      else if(pid_JAM==3222) ptmp=particleTable->FindParticle("sigma+");
      else if(pid_JAM==3312) ptmp=particleTable->FindParticle("xi-");
      else if(pid_JAM==3322) ptmp=particleTable->FindParticle("xi0");
      else if(pid_JAM==-211) ptmp=particleTable->FindParticle("pi-");
      else if(pid_JAM==111) ptmp=particleTable->FindParticle("pi0");
      else if(pid_JAM==211) ptmp=particleTable->FindParticle("pi+");
      else if(pid_JAM==221) ptmp=particleTable->FindParticle("eta");
      else if(pid_JAM==22) ptmp=particleTable->FindParticle("gamma");
      else if(pid_JAM==-321) ptmp=particleTable->FindParticle("kaon-");
      else if(pid_JAM==321) ptmp=particleTable->FindParticle("kaon+");
      else if(pid_JAM==311){
	if(G4UniformRand()>0.5)
	  ptmp=particleTable->FindParticle("kaon0S");
	else
	  ptmp=particleTable->FindParticle("kaon0L");
      }
      else if(pid_JAM==-311){
	if(G4UniformRand()>0.5)
	  ptmp=particleTable->FindParticle("kaon0S");
	else
	  ptmp=particleTable->FindParticle("kaon0L");
      }
      else{
	std::cout<<"pid_JAM="<<pid_JAM<<std::endl;
	//getchar();
	continue;
      }
      G4ThreeVector LPos = G4ThreeVector(0.,0.,gGeom.GetGlobalPosition("Target").z());
      //G4ThreeVector LP = GetP_JAM(Nbeam_JAMInput, inp);
      G4ThreeVector LP = G4ThreeVector(px[inp]*CLHEP::GeV, py[inp]*CLHEP::GeV, pz[inp]*CLHEP::GeV);

      m_particle_gun->SetParticleDefinition(ptmp);
      m_particle_gun->SetParticlePosition(LPos);
      m_particle_gun->SetParticleMomentum(LP);
      m_particle_gun->GeneratePrimaryVertex(anEvent);

      //    std::cout<<"encording="<<ptmp->GetPDGEncoding()<<std::endl;
      // if(pid_JAM==-311)
      //   gAnaMan.SetPrimaryParticle(inp,LP.x(),LP.y(),LP.z(),
      //                              -1.*ptmp->GetPDGMass()/CLHEP::GeV, ptmp->GetPDGEncoding());
      // else
      //   gAnaMan.SetPrimaryParticle(inp,LP.x(),LP.y(),LP.z(),
      //                              ptmp->GetPDGMass()/CLHEP::GeV, ptmp->GetPDGEncoding());
      // //gAnaMan.SetPrimaryVertex(0,LPos.x(),LPos.y(),LPos.z());
      // gAnaMan.SetPrimaryVertex(inp,LPos.x(),LPos.y(),LPos.z());

      //delete ptmp;
    }
  }

  ++Nbeam_JAM;
}

//_____________________________________________________________________________
//reactio No #3104 JAM input K0bar
void
PrimaryGeneratorAction::GenerateJAMInputK0bar(G4Event* anEvent, TTree*t1)
{
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  // G4ParticleDefinition* kaonMinus = particleTable->FindParticle("kaon-");
  G4double pbeam=gConf.Get<G4double>("BeamMom");
  //  pbeam=1.8;
  // G4double mom_kp_x = 0;
  // G4double mom_kp_y = 0;
  // G4double mom_kp_z = pbeam;
  // gAnaMan.SetPrimaryBeam(0,0,pbeam);

  int N_first = gConf.Get<G4int>("NbeamFirst");

  int Nbeam_JAMInput = Nbeam_JAM+N_first;

  int ev_max= t1->GetEntries();
  int np;
  int pid[20];
  G4double  px[20], py[20], pz[20];
  t1->SetBranchAddress("np", &np);
  t1->SetBranchAddress("pid", pid);
  t1->SetBranchAddress("px", px);
  t1->SetBranchAddress("py", py);
  t1->SetBranchAddress("pz", pz);

  t1->GetEntry(Nbeam_JAMInput%ev_max);
  G4int np_JAM = np;

  //  G4int np_JAM = Getnp_JAM(Nbeam_JAMInput);
  // gAnaMan.SetPrimaryInfo(0., 0., 0., 0., 0.);

  std::cout<<"Nbeam_JAMInpu="<<Nbeam_JAMInput
	   <<", Num of Primary Particle="<<np_JAM<<std::endl;;
  G4int K0bflag=0;
  for(int inp=0; inp<np_JAM; ++inp){
    G4int pid_JAM = pid[inp];
    if(pid_JAM==311)
      ++K0bflag;
  }

  if(K0bflag>0){
    for(int inp=0; inp<np_JAM; ++inp){
      G4ParticleDefinition *ptmp;
      G4int pid_JAM = pid[inp];
      //G4int pid_JAM = GetPID_JAM(Nbeam_JAMInput, inp);
      if(pid_JAM==2212) ptmp=particleTable->FindParticle("proton");
      else if(pid_JAM==2112) ptmp=particleTable->FindParticle("neutron");
      else if(pid_JAM==3122) ptmp=particleTable->FindParticle("lambda");
      else if(pid_JAM==3112) ptmp=particleTable->FindParticle("sigma-");
      else if(pid_JAM==3212) ptmp=particleTable->FindParticle("sigma0");
      else if(pid_JAM==3222) ptmp=particleTable->FindParticle("sigma+");
      else if(pid_JAM==3312) ptmp=particleTable->FindParticle("xi-");
      else if(pid_JAM==3322) ptmp=particleTable->FindParticle("xi0");
      else if(pid_JAM==-211) ptmp=particleTable->FindParticle("pi-");
      else if(pid_JAM==111) ptmp=particleTable->FindParticle("pi0");
      else if(pid_JAM==211) ptmp=particleTable->FindParticle("pi+");
      else if(pid_JAM==221) ptmp=particleTable->FindParticle("eta");
      else if(pid_JAM==22) ptmp=particleTable->FindParticle("gamma");
      else if(pid_JAM==-321) ptmp=particleTable->FindParticle("kaon-");
      else if(pid_JAM==321) ptmp=particleTable->FindParticle("kaon+");
      else if(pid_JAM==311){
	if(G4UniformRand()>0.5)
	  ptmp=particleTable->FindParticle("kaon0S");
	else
	  ptmp=particleTable->FindParticle("kaon0L");
      }
      else if(pid_JAM==-311){
	if(G4UniformRand()>0.5)
	  ptmp=particleTable->FindParticle("kaon0S");
	else
	  ptmp=particleTable->FindParticle("kaon0L");
      }
      else{
	std::cout<<"pid_JAM="<<pid_JAM<<std::endl;
	//getchar();
	continue;
      }
      G4ThreeVector LPos = G4ThreeVector(0.,0.,gGeom.GetGlobalPosition("Target").z());
      //G4ThreeVector LP = GetP_JAM(Nbeam_JAMInput, inp);
      G4ThreeVector LP = G4ThreeVector(px[inp]*CLHEP::GeV, py[inp]*CLHEP::GeV, pz[inp]*CLHEP::GeV);

      m_particle_gun->SetParticleDefinition(ptmp);
      m_particle_gun->SetParticlePosition(LPos);
      m_particle_gun->SetParticleMomentum(LP);
      m_particle_gun->GeneratePrimaryVertex(anEvent);

      //    std::cout<<"encording="<<ptmp->GetPDGEncoding()<<std::endl;
      // if(pid_JAM==-311)
      //   gAnaMan.SetPrimaryParticle(inp,LP.x(),LP.y(),LP.z(),
      //                              -1.*ptmp->GetPDGMass()/CLHEP::GeV, ptmp->GetPDGEncoding());
      // else
      //   gAnaMan.SetPrimaryParticle(inp,LP.x(),LP.y(),LP.z(),
      //                              ptmp->GetPDGMass()/CLHEP::GeV, ptmp->GetPDGEncoding());
      // //gAnaMan.SetPrimaryVertex(0,LPos.x(),LPos.y(),LPos.z());
      // gAnaMan.SetPrimaryVertex(inp,LPos.x(),LPos.y(),LPos.z());

      //delete ptmp;
    }
  }

  ++Nbeam_JAM;
}

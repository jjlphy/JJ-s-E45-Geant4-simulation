// -*- C++ -*-

#include "SteppingAction.hh"

#include <unordered_map>

#include <G4Material.hh>
#include <G4ParticleDefinition.hh>
#include <G4ParticleTypes.hh>
#include <G4SteppingManager.hh>
#include <G4Step.hh>
#include <G4StepPoint.hh>
#include <G4Track.hh>
#include <G4TrackStatus.hh>
#include <G4VPhysicalVolume.hh>
#include <G4RunManager.hh> 

#include "ConfMan.hh"
#include "PrintHelper.hh"
#include "AnaManager.hh"
#include "PrimaryGeneratorAction.hh"

namespace
{
auto& gAnaMan = AnaManager::GetInstance();
const auto& gConf = ConfMan::GetInstance();
}

//_____________________________________________________________________________
SteppingAction::SteppingAction()
  : G4UserSteppingAction()
{
}

//_____________________________________________________________________________
SteppingAction::~SteppingAction()
{
}

//_____________________________________________________________________________
void
SteppingAction::UserSteppingAction(const G4Step* theStep)
{
  static const G4bool KillStepInIron = gConf.Get<G4bool>("KillStepInIron");

  auto theTrack = theStep->GetTrack();
  auto theParticle = theTrack->GetParticleDefinition();
  auto parentID = theTrack->GetParentID();
  auto particleName = theParticle->GetParticleName();
  auto particlePdgCode = theParticle->GetPDGEncoding();
  auto particleMass = theParticle->GetPDGMass();
  auto prePoint = theStep->GetPreStepPoint();
  auto prePV = prePoint->GetPhysicalVolume();
  auto prePVName = prePV->GetName();
  auto postPoint = theStep->GetPostStepPoint();
  auto theProcess = postPoint->GetProcessDefinedStep()->GetProcessName();
  auto stepLength = theTrack->GetStepLength();
  G4ThreeVector stepMiddlePosition = (prePoint->GetPosition() + postPoint->GetPosition())/2.0;
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

  // -- cal effective thickness -----
  if (prePVName == "TargetPV") {
    G4int generator = gAnaMan.GetNextGenerator();
    if (generator == 7201 && particleName == "kaon-") {
      G4double effective_thickness = gAnaMan.GetEffectiveThickness();
      if (effective_thickness == -1.0) gAnaMan.SetEffectiveThickness(stepLength);
      else gAnaMan.SetEffectiveThickness( (G4double) effective_thickness+stepLength);
    }
  }
  
  // -- check decay particle -----
  std::pair<G4String, G4String> previous_particle = gAnaMan.GetPreviousParticle();
  G4ThreeVector previous_step_pos = gAnaMan.GetDecayPosition();
  G4int generator = gAnaMan.GetNextGenerator();
  if (previous_particle.second == "Decay" && previous_particle.first == gAnaMan.GetFocusParticle(generator) ){
    if ( gAnaMan.IsInsideHtof(previous_step_pos) ) gAnaMan.SetDecayParticleCode( particlePdgCode );
    gAnaMan.SetFocusParentID( parentID );
  }
  
  gAnaMan.SetPreviousParticle(particleName, theProcess);
  gAnaMan.SetDecayPosition(stepMiddlePosition);
  
  // -- Get Seconday Vertex info --
  PrimaryGeneratorAction* generatorAction = (PrimaryGeneratorAction*) G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction();
  
  for(int i=0;i<10;i++){
    if(particlePdgCode != generatorAction->m_primary_pdg[i]){
      continue;
    }else if(particlePdgCode == generatorAction->m_primary_pdg[i]){
      if(theTrack->GetTrackStatus() == fStopAndKill){
	const std::vector<const G4Track*>* secTracks = theStep->GetSecondaryInCurrentStep();
	if (!secTracks->empty()) {
	  for (const auto& secTrack : *secTracks) {
            if (secTrack->GetCreatorProcess()) {
	      G4String secProcessName = secTrack->GetCreatorProcess()->GetProcessName();
	      G4int motherPdgCode = particlePdgCode;
	      G4int daughterPdgCode = secTrack->GetDefinition()->GetPDGEncoding();
	      G4ThreeVector mom_se = secTrack->GetMomentum();
	      G4LorentzVector v_se(secTrack->GetPosition(), 0);
	      G4LorentzVector p_se(mom_se, std::sqrt(std::pow(particleMass,2)+std::pow(mom_se.mag(),2)));
	      
	      gAnaMan.SetSecondaryVertex(daughterPdgCode,motherPdgCode,p_se,v_se);
	      
            } else {
	      G4cout << "Secondary particle has no creator process!" << G4endl;
            }
	  }
	}
      }
    }
  }
  
#ifdef DEBUG
  PrintHelper helper(3, std::ios::fixed, G4cout);
  auto time = prePoint->GetGlobalTime();
  auto track_id = theTrack->GetTrackID();
  std::stringstream particle_ss;
  particle_ss << particleName << "(" << track_id << ") ";
  if(theProcess != "eIoni" &&
     theProcess != "hIoni" &&
     theProcess != "msc" &&
     theProcess != "eBeam" &&
     theProcess != "Transportation"){
    G4cout << "   " << time/CLHEP::ns << " ns : "
	   << particle_ss.str() << theProcess << G4endl;
  }

  auto secondary = theStep->GetSecondaryInCurrentStep();
  for(const auto& s : *secondary){
    auto particle = s->GetDefinition();
    auto name = particle->GetParticleName();
    // auto type = particle->GetParticleType();
    if(true
       || (particleName == "lambda" && name == "proton")){
      G4cout << "   " << particle_ss.str() << "\tP" << prePoint->GetMomentum()
	     << " X" << prePoint->GetPosition()
	     << "\t-> " << name << " P" << s->GetMomentum() << G4endl;
    }
  }

  if(false
     && particleName == "proton"){
    auto preMaterial = prePoint->GetMaterial();
    G4double edep = theStep->GetTotalEnergyDeposit();
    G4cout << "   " << particleName << " " << theProcess
	   << " " << theTrack->GetTrackStatus()
	   << " " << preMaterial->GetName()
	   << " x=" << prePoint->GetPosition()
	   << " p=" << prePoint->GetMomentum()
	   << " edep=" << edep
	   << G4endl;
  }
#endif

  if(KillStepInIron){
    auto preMaterial = prePoint->GetMaterial();
    if(preMaterial->GetName() == "Iron"){
      theTrack->SetTrackStatus(fStopAndKill);
      return;
    }
  }

  if(particleName == "e-" || particleName == "e+"){
    theTrack->SetTrackStatus(fStopAndKill);
    return;
  }

  // if(prePVName.contains("Coil") || prePVName.contains("Guard")){
  //   theTrack->SetTrackStatus(fStopAndKill);
  //   return;
  // }
}

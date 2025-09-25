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

// === [added] for spawning secondaries at target entry =========
#include <G4EventManager.hh>
#include <G4StackManager.hh>
#include <Randomize.hh>           // G4UniformRand
// ===============================================================

#include "ConfMan.hh"
#include "PrintHelper.hh"
#include "AnaManager.hh"
#include "PrimaryGeneratorAction.hh"

#include "TrackTag.hh"
#include "DCGeomMan.hh"
#include "DetSizeMan.hh"

#include "TLorentzVector.h"//jaejin
#include "TGenPhaseSpace.h"//jaejin

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

// ------------------------- [added helper only] ------------------------------
// 3-body phase space spawner at a given vertex (lab frame).
// channel=0 : pi+ pi- n   , channel=1 : pi- pi0 p
namespace {
  void Spawn3BodyAt(const G4Track* motherTrack, int channel,
                    const G4ThreeVector& vtx, G4StackManager* stackMan)
  {
    using namespace CLHEP;

    const G4ThreeVector p3 = motherTrack->GetMomentum(); // (Geant4 unit)

    auto pt = G4ParticleTable::GetParticleTable();
    const double mPim = pt->FindParticle("pi-")->GetPDGMass()/GeV;
    const double mPip = pt->FindParticle("pi+")->GetPDGMass()/GeV;
    const double mPi0 = pt->FindParticle("pi0")->GetPDGMass()/GeV;
    const double mN   = pt->FindParticle("neutron")->GetPDGMass()/GeV;
    const double mP   = pt->FindParticle("proton")->GetPDGMass()/GeV;

    const double pGeV = p3.mag()/GeV;
    TLorentzVector LVpi ( p3.x()/GeV, p3.y()/GeV, p3.z()/GeV, std::hypot(pGeV, mPim) );
    TLorentzVector LVpro( 0.,0.,0., mP );
    TLorentzVector W = LVpi + LVpro;

    double masses[3];
    if(channel==0){ masses[0]=mPip; masses[1]=mPim; masses[2]=mN; }
    else          { masses[0]=mPim; masses[1]=mPi0; masses[2]=mP; }

    if(W.M() <= masses[0] + masses[1] + masses[2]) return;

    TGenPhaseSpace gen;
    gen.SetDecay(W, 3, masses);

    // small downstream bias (optional)
    int tries=0;
    for(;;){
      gen.Generate();
      bool ok=false;
      for(int i=0;i<3;i++){
        const TLorentzVector* d = gen.GetDecay(i);
        if(d && d->Pz()>0){ ok=true; break; }
      }
      if(ok || ++tries>100) break;
    }

    struct Out { const char* name; } outs0[3]={{"pi+"},{"pi-"},{"neutron"}};
    struct Out  outs1[3]={{"pi-"},{"pi0"},{"proton"}};
    const auto* outs = (channel==0)? outs0 : outs1;

    const double t0 = motherTrack->GetGlobalTime();
    for(int i=0;i<3;i++){
      const TLorentzVector* d = gen.GetDecay(i);
      auto def = pt->FindParticle(outs[i].name);
      auto dyn = new G4DynamicParticle(def,
                    G4ThreeVector(d->Px()*GeV, d->Py()*GeV, d->Pz()*GeV));
      auto trk = new G4Track(dyn, t0, vtx);
      trk->SetParentID(motherTrack->GetTrackID());
      trk->SetTrackStatus(fAlive);
      stackMan->PushOneTrack(trk);
    }
  }
}
// ----------------------------------------------------------------------------

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

//---------수정 부분 25.9.25--------------
// ------------------------------------------------------------------
// [UPSTREAM_HELPER kill at VP4 plane without geometry]
//  - SimpleTag("UPSTREAM_HELPER") 가 달린 트랙만 대상
//  - VP4(z=z_vp4) 평면을 통과할 때, 타겟 원기둥 반경 밖이면 StopAndKill
// ------------------------------------------------------------------
{
  // 1) 태그 확인 (아직 StackingAction에서 안 달았으면 nullptr라서 그냥 패스)
  const auto* tag = dynamic_cast<const SimpleTag*>(theTrack->GetUserInformation());
  if (tag && tag->why_ && std::string(tag->why_) == "UPSTREAM_HELPER") {

    // 2) VP4 z 위치와 타겟 치수 얻기
    static const auto& gGeom = DCGeomMan::GetInstance();
    static const auto& gSize = DetSizeMan::GetInstance();

    // VP4 중심 좌표 (전역좌표, mm)
    static const G4ThreeVector vp4_pos = gGeom.GetGlobalPosition("VP4"); // (0,0,-143) mm 설정 기준
    const double z_vp4 = vp4_pos.z();

    // 타겟 중심(= SHSTarget)과 치수
    static const G4ThreeVector tgt_c = gGeom.GetGlobalPosition("SHSTarget"); // 보통 (0,0,-143) mm
    // DetSizeMan: Target 항목이 "Rin Rout DZ" (full length) 포맷이라고 가정
    static const G4ThreeVector tgt_sz = gSize.GetSize("Target");
    const double Rin = tgt_sz.x();        // 보통 0.0
    const double Rout = tgt_sz.y();       // (단위: mm) 외반경(‘full’ 표기지만 여기선 반경 값으로 쓰임)
    const double DZ_full = tgt_sz.z();    // (단위: mm) full 길이
    const double half_dz = 0.5 * DZ_full;

    // 3) 이번 스텝이 VP4 평면을 '가로지르는지' 판정
    const double z1 = prePoint->GetPosition().z();
    const double z2 = postPoint->GetPosition().z();

    // 선형보간으로 교차점 계산 (z=z_vp4를 지날 때만)
    if ( (z1 - z_vp4) * (z2 - z_vp4) <= 0.0 && (z1 != z2) ) {

      const double t = (z_vp4 - z1) / (z2 - z1); // [0,1] 사이
      const G4ThreeVector x1 = prePoint->GetPosition();
      const G4ThreeVector x2 = postPoint->GetPosition();
      const G4ThreeVector x_at = x1 + t*(x2 - x1); // VP4에서의 교차점

      // 타겟 원기둥 단순 in/out 체크 (지오메트리 없이)
      const G4ThreeVector rel = x_at - tgt_c;
      const double r_xy = std::hypot(rel.x(), rel.y());
      const bool inside_rad = (r_xy <= Rout + 1e-6);        // 약간의 수치 여유
      const bool inside_z   = (std::abs(rel.z()) <= half_dz + 1e-6);

      const bool hits_target_cylinder = (inside_rad && inside_z);

      // 4) 타겟 원기둥 영역을 지나지 않으면, 이 트랙은 분석용 helper이므로 Kill
      if (!hits_target_cylinder) {
        theTrack->SetTrackStatus(fStopAndKill);
        // (선택) 분석기록: 원하면 여기에 gAnaMan.SetDoHitTGT(false) 같은 플래그를 달 수 있음
        return;
      } else {
        // 통과 = 진짜 타겟 통과로 간주 → 계속 추적 (아무 것도 하지 않음)
      }
    }
  }
}
//---------수정 부분 25.9.25--------------

//---------------------- [added] Target first-entry → 3-body spawn ----------------------
// - UPSTREAM_HELPER(=primary pi-)만 대상
// - "첫 진입": 경계(fGeomBoundary) && pre!=TargetPV && post==TargetPV
{
  const auto* tag = dynamic_cast<const SimpleTag*>(theTrack->GetUserInformation());
  if (tag && tag->why_ && std::string(tag->why_)=="UPSTREAM_HELPER") {
    auto postPV = postPoint->GetPhysicalVolume();
    auto prePV  = prePoint->GetPhysicalVolume();
    if (postPV && prePV) {
      const auto& postName = postPV->GetName();
      const auto& preName  = prePV->GetName();

      if (postPoint->GetStepStatus()==fGeomBoundary
          && postName=="TargetPV" && preName!="TargetPV")
      {
        // 진입점에서 3체 생성
        const G4ThreeVector vtx = postPoint->GetPosition();
        int channel = (G4UniformRand()<0.5)? 0:1; // 0: pi+ pi- n  , 1: pi- pi0 p

        auto stackMan = G4EventManager::GetEventManager()->GetStackManager();
        if (stackMan) {
          Spawn3BodyAt(theTrack, channel, vtx, stackMan);
        }

        // 원 빔 Kill → downstream으로는 오직 딸입자만
        theTrack->SetTrackStatus(fStopAndKill);
        return;
      }
    }
  }
}
//---------------------------------------------------------------------------------------

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

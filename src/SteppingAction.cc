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

// === [added] for spawning secondaries & fail-safe planes =========
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

// ==== SCENARIO SWITCH =======================================================
// 1) 강제 2π 생성 모드 (현재 우리가 쓰던 모드)
//    - 타겟 첫 진입에서 (π+ π− n) 또는 (π− π0 p) 생성하고 원빔 kill
//    - 타겟 원기둥을 빗나간 원빔은 VP4에서 kill, VP5 fail-safe도 동작
//
// 2) 빔-through 모드
//    - 강제 생성 없음, VP4/VP5 가드 없음 → 원빔은 물리 그대로 진행
//    - (분석용으로 “진짜 빔-through 이벤트”를 수집할 때 사용)
// ============================================================================

// #define E45_SCENARIO_FORCE_2PI  // ← 이 줄이 있으면 “강제 2π 생성 모드”
// (주석 처리하면 빔-through 모드)
#define E45_SCENARIO_FORCE_2PI   // <-- 기본값: 강제 2π 생성

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

// ------------------------- helper: 3-body spawner ---------------------------
namespace {
  // channel=0 : pi+ pi- n   , channel=1 : pi- pi0 p
  void Spawn3BodyAt(const G4Track* motherTrack, int channel,
                    const G4ThreeVector& vtx, G4StackManager* stackMan)
  {
    using namespace CLHEP;

    const G4ThreeVector p3 = motherTrack->GetMomentum(); // Geant4 unit

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

  // [helper] “원빔” 판정자: ParentID==0 && PDG==-211 (primary π−)
  const bool isPrimaryBeam = (theTrack->GetParentID()==0 && particlePdgCode==-211);

  // --------- (원본 코드) 지표 수집/기존 기능 ---------
  if (prePVName == "TargetPV") {
    G4int generator = gAnaMan.GetNextGenerator();
    if (generator == 7201 && particleName == "kaon-") {
      G4double effective_thickness = gAnaMan.GetEffectiveThickness();
      if (effective_thickness == -1.0) gAnaMan.SetEffectiveThickness(stepLength);
      else gAnaMan.SetEffectiveThickness((G4double)effective_thickness + stepLength);
    }
  }

  std::pair<G4String, G4String> previous_particle = gAnaMan.GetPreviousParticle();
  G4ThreeVector previous_step_pos = gAnaMan.GetDecayPosition();
  G4int generator = gAnaMan.GetNextGenerator();
  if (previous_particle.second == "Decay" && previous_particle.first == gAnaMan.GetFocusParticle(generator) ){
    if ( gAnaMan.IsInsideHtof(previous_step_pos) ) gAnaMan.SetDecayParticleCode( particlePdgCode );
    gAnaMan.SetFocusParentID( parentID );
  }
  gAnaMan.SetPreviousParticle(particleName, theProcess);
  gAnaMan.SetDecayPosition(stepMiddlePosition);

  // -- Secondary Vertex 기록(원본) --
  PrimaryGeneratorAction* generatorAction =
    (PrimaryGeneratorAction*) G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction();

  for(int i=0;i<10;i++){
    if(particlePdgCode != generatorAction->m_primary_pdg[i]) continue;
    if(theTrack->GetTrackStatus() == fStopAndKill){
      const std::vector<const G4Track*>* secTracks = theStep->GetSecondaryInCurrentStep();
      if (!secTracks->empty()) {
        for (const auto& secTrack : *secTracks) {
          if (secTrack->GetCreatorProcess()) {
            G4int motherPdgCode = particlePdgCode;
            G4int daughterPdgCode = secTrack->GetDefinition()->GetPDGEncoding();
            G4ThreeVector mom_se = secTrack->GetMomentum();
            G4LorentzVector v_se(secTrack->GetPosition(), 0);
            G4LorentzVector p_se(mom_se, std::sqrt(std::pow(particleMass,2)+std::pow(mom_se.mag(),2)));
            gAnaMan.SetSecondaryVertex(daughterPdgCode,motherPdgCode,p_se,v_se);
          }
        }
      }
    }
  }

#ifdef DEBUG
  // ... (DEBUG 블록은 원본 그대로 유지)
#endif

//=========================== 커스텀 로직 (토글) ===========================

#ifdef E45_SCENARIO_FORCE_2PI
  // ==================== [모드 1: 강제 2π 생성] ====================

  // (A) VP4 평면 miss → 원빔 kill
  {
    const auto* tag = dynamic_cast<const SimpleTag*>(theTrack->GetUserInformation());
    if ( (tag && tag->why_ && std::string(tag->why_) == "UPSTREAM_HELPER") || isPrimaryBeam ) {

      static const auto& gGeom = DCGeomMan::GetInstance();
      static const auto& gSize = DetSizeMan::GetInstance();

      static const G4ThreeVector vp4_pos = gGeom.GetGlobalPosition("VP4");
      const double z_vp4 = vp4_pos.z();

      static const G4ThreeVector tgt_c  = gGeom.GetGlobalPosition("SHSTarget");
      static const G4ThreeVector tgt_sz = gSize.GetSize("Target");
      const double Rout = tgt_sz.y();
      const double half_dz = 0.5 * tgt_sz.z();

      const double z1 = prePoint->GetPosition().z();
      const double z2 = postPoint->GetPosition().z();

      if ( (z1 - z_vp4) * (z2 - z_vp4) <= 0.0 && (z1 != z2) ) {
        const double t = (z_vp4 - z1) / (z2 - z1);
        const G4ThreeVector x1 = prePoint->GetPosition();
        const G4ThreeVector x2 = postPoint->GetPosition();
        const G4ThreeVector x_at = x1 + t*(x2 - x1);

        const G4ThreeVector rel = x_at - tgt_c;
        const double r_xy = std::hypot(rel.x(), rel.y());
        const bool inside_rad = (r_xy <= Rout + 1e-6);
        const bool inside_z   = (std::abs(rel.z()) <= half_dz + 1e-6);

        if (!(inside_rad && inside_z)) {
          theTrack->SetTrackStatus(fStopAndKill);
          return;
        }
      }
    }
  }

  // (B) 타겟 첫 진입에서 3체 생성 후 원빔 kill
  {
    const auto* tag = dynamic_cast<const SimpleTag*>(theTrack->GetUserInformation());
    if ( (tag && tag->why_ && std::string(tag->why_)=="UPSTREAM_HELPER") || isPrimaryBeam ) {
      auto postPV = postPoint->GetPhysicalVolume();
      auto prePV  = prePoint->GetPhysicalVolume();
      if (postPV && prePV) {
        const G4Material* preMat  = prePoint->GetMaterial();
        const G4Material* postMat = postPoint->GetMaterial();
        const auto preMatName  = preMat  ? preMat ->GetName() : "";
        const auto postMatName = postMat ? postMat->GetName() : "";
        const bool isGeomBoundary = (postPoint->GetStepStatus()==fGeomBoundary);

        const bool enterLH2 =
          isGeomBoundary && postMat && (postMatName=="LH2") && (!preMat || preMatName!="LH2");

        const auto& postName = postPV->GetName();
        const auto& preName  = prePV ->GetName();
        const bool enterTargetLike =
          isGeomBoundary &&
          (postName.find("Target")!=std::string::npos) &&
          (preName.find("Target")==std::string::npos);

        if (enterLH2 || enterTargetLike)
        {
          const G4ThreeVector vtx = postPoint->GetPosition();
          int channel = (G4UniformRand()<0.5)? 0:1; // 0: π+π−n, 1: π−π0p

          auto stackMan = G4EventManager::GetEventManager()->GetStackManager();
          if (stackMan) Spawn3BodyAt(theTrack, channel, vtx, stackMan);

          theTrack->SetTrackStatus(fStopAndKill); // 원빔 제거
          return;
        }
      }
    }
  }

  // (C) VP4 downstream 하드 가드 (LH2 바깥에서만)
  {
    if (isPrimaryBeam) {
      static const auto& gGeom = DCGeomMan::GetInstance();
      static const G4ThreeVector vp4_pos = gGeom.GetGlobalPosition("VP4");
      const double z_vp4  = vp4_pos.z();
      const double z_post = postPoint->GetPosition().z();

      const G4Material* postMat = postPoint->GetMaterial();
      const auto postMatName = postMat ? postMat->GetName() : "";

      if (postMatName!="LH2" && z_post > z_vp4 + 0.5*CLHEP::mm) {
        theTrack->SetTrackStatus(fStopAndKill);
        return;
      }
    }
  }

  // (D) VP5 fail-safe
  {
    const auto* tag = dynamic_cast<const SimpleTag*>(theTrack->GetUserInformation());
    if ( (tag && tag->why_ && std::string(tag->why_) == "UPSTREAM_HELPER") || isPrimaryBeam ) {
      static const auto& gGeom = DCGeomMan::GetInstance();
      static const G4ThreeVector vp5_pos = gGeom.GetGlobalPosition("VP5");
      const double z_vp5 = vp5_pos.z();

      const double z1 = prePoint->GetPosition().z();
      const double z2 = postPoint->GetPosition().z();
      if ( (z1 - z_vp5) * (z2 - z_vp5) <= 0.0 && (z1 != z2) ) {
        theTrack->SetTrackStatus(fStopAndKill);
        return;
      }
    }
  }

#else
  // ==================== [모드 2: 빔-through] ====================
  // 주의:
  //  - 강제 3체 생성 없음
  //  - VP4 miss-kill / VP4 하드가드 / VP5 fail-safe 없음
  //  - 따라서 원빔(π−, ParentID=0)은 타겟을 통과/산란/흡수 모두 “물리 모델 그대로” 진행
  //  - 전자/양전자/철 내부 kill 같은 기존 일반 컷은 아래 원본 코드에서 계속 적용
#endif

//========================= 커스텀 로직 끝 ==========================

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

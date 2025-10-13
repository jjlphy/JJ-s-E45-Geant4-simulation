// -*- C++ -*-
// SteppingAction.cc  (E45)
// - 강제 2π 생성 모드(옵션) + 채널 선택(0/1/2)을 conf/환경변수로 런타임 반영
// - 원빔/자식 출처 태깅: OriginTag 적용 (kPrimaryBeam / kForced2PiChild)
// - (하위호환) SimpleTag("UPSTREAM_HELPER")도 fallback으로 인식
// jaejin 2025-10-13: Spawn3BodyAt에서 자식들에 OriginTag 부여,
//                    IsPrimaryBeamLike() 공통화, 로직 안전성/가독성 개선.

#include "SteppingAction.hh"

#include <unordered_map>
#include <string>
#include <cmath>    // std::hypot
#include <cstdlib>  // std::getenv

#include <G4Material.hh>
#include <G4ParticleDefinition.hh>
#include <G4ParticleTypes.hh>
#include <G4ParticleTable.hh>
#include <G4SteppingManager.hh>
#include <G4Step.hh>
#include <G4StepPoint.hh>
#include <G4Track.hh>
#include <G4TrackStatus.hh>
#include <G4VPhysicalVolume.hh>
#include <G4RunManager.hh>
#include <G4EventManager.hh>
#include <G4StackManager.hh>

#include <Randomize.hh>           // G4UniformRand
#include <CLHEP/Units/SystemOfUnits.h>

#include "ConfMan.hh"
#include "PrintHelper.hh"
#include "AnaManager.hh"
#include "PrimaryGeneratorAction.hh"

#include "TrackTag.hh"   // OriginTag / SimpleTag
#include "DCGeomMan.hh"
#include "DetSizeMan.hh"

#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"

// ==== SCENARIO SWITCH =======================================================
// 1) 강제 2π 생성 모드 (우리가 쓰던 모드)
// 2) 빔-through 모드
#define E45_SCENARIO_FORCE_2PI
// ============================================================================

namespace
{
  auto&       gAnaMan = AnaManager::GetInstance();
  const auto& gConf   = ConfMan::GetInstance();

  // 런타임 동적 읽기: 환경변수(FORCE_2PI_MODE) > conf > 기본값(2: mixture)
  inline G4int GetForce2PiMode()
  {
    if(const char* env = std::getenv("FORCE_2PI_MODE")){
      int v = std::atoi(env);
      if(v==0 || v==1) return v;
      return 2;
    }
    int raw = 2;
    try { raw = gConf.Get<G4int>("Force2PiMode"); }
    catch(...) { raw = 2; }
    return (raw==0 || raw==1) ? raw : 2;
  }

  // 실행 시작 시 1회만 실제 모드 출력
  inline void PrintForce2PiModeOnce()
  {
    static bool s_printOnce = true;
    if(!s_printOnce) return;
    const int mode = GetForce2PiMode();
    G4cout << "[E45] Force2PiMode(read)=" << mode
#ifdef E45_SCENARIO_FORCE_2PI
           << "  (SCENARIO_FORCE_2PI=ON)"
#else
           << "  (SCENARIO_FORCE_2PI=OFF)"
#endif
           << G4endl;
    s_printOnce = false;
  }

  // ──────────────────────────────────────────────────────────────────────────
  // 공통 판정: 이 track이 "원빔처럼" 취급되어야 하는가?
  //  1) OriginTag(kPrimaryBeam) 우선
  //  2) (하위호환) SimpleTag("UPSTREAM_HELPER")
  //  3) (백업) parentID==0 && pdg==-211
  // ──────────────────────────────────────────────────────────────────────────
  inline bool IsPrimaryBeamLike(const G4Track* trk){
    if(!trk) return false;
    if(HasOrigin(trk, EOrigin::kPrimaryBeam)) return true;

    // fallback: SimpleTag
    if(auto* s = dynamic_cast<const SimpleTag*>(trk->GetUserInformation())){
      if(s->why_ && std::string(s->why_)=="UPSTREAM_HELPER") return true;
    }

    // backup: 물리적 정의
    const auto* pd = trk->GetParticleDefinition();
    const int pdg = pd ? pd->GetPDGEncoding() : 0;
    return (trk->GetParentID()==0 && pdg==-211);
  }

} // namespace

//_____________________________________________________________________________
SteppingAction::SteppingAction()
  : G4UserSteppingAction()
{
  // conf 반영 여부를 눈으로 확인
  PrintForce2PiModeOnce();
}

SteppingAction::~SteppingAction()
{
}

// ------------------------- helper: 3-body spawner ---------------------------
// channel=0 : pi+ pi- n   , channel=1 : pi- pi0 p
namespace {

// jaejin 2025-10-13: Spawn3BodyAt에서 생성 자식에 OriginTag(kForced2PiChild)를 부여.
//                    name/채널 저장으로 오프라인 필터가 쉬워진다.
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
  gen.Generate(); // unbiased single sample

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

    // NEW: 출처 태그 부여 (강제 2π 자식)
    trk->SetUserInformation(new OriginTag(EOrigin::kForced2PiChild, channel, outs[i].name));

    stackMan->PushOneTrack(trk);
  }
}

} // namespace

// ----------------------------------------------------------------------------

void
SteppingAction::UserSteppingAction(const G4Step* theStep)
{
  using namespace CLHEP;

  static const G4bool KillStepInIron = gConf.Get<G4bool>("KillStepInIron");

  auto theTrack    = theStep->GetTrack();
  auto theParticle = theTrack->GetParticleDefinition();
  auto parentID    = theTrack->GetParentID();
  auto particleName= theParticle->GetParticleName();
  auto particlePdg = theParticle->GetPDGEncoding();
  auto particleMass= theParticle->GetPDGMass();

  auto prePoint    = theStep->GetPreStepPoint();
  auto postPoint   = theStep->GetPostStepPoint();

  auto prePV       = prePoint->GetPhysicalVolume();
  auto prePVName   = prePV ? prePV->GetName() : "";

  auto theProcess  = postPoint->GetProcessDefinedStep()
                   ? postPoint->GetProcessDefinedStep()->GetProcessName()
                   : "";

  auto stepLength  = theTrack->GetStepLength();
  G4ThreeVector stepMiddlePosition =
      (prePoint->GetPosition() + postPoint->GetPosition())/2.0;

  const bool isPrimaryBeam = (parentID==0 && particlePdg==-211); // backup 정의

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
  if (previous_particle.second == "Decay" &&
      previous_particle.first  == gAnaMan.GetFocusParticle(generator) ){
    if ( gAnaMan.IsInsideHtof(previous_step_pos) )
      gAnaMan.SetDecayParticleCode( theParticle->GetPDGEncoding() );
    gAnaMan.SetFocusParentID( parentID );
  }
  gAnaMan.SetPreviousParticle(particleName, theProcess);
  gAnaMan.SetDecayPosition(stepMiddlePosition);

  // -- Secondary Vertex 기록(원본) --
  PrimaryGeneratorAction* generatorAction =
    (PrimaryGeneratorAction*) G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction();

  for(int i=0;i<10;i++){
    if(particlePdg != generatorAction->m_primary_pdg[i]) continue;
    if(theTrack->GetTrackStatus() == fStopAndKill){
      const std::vector<const G4Track*>* secTracks = theStep->GetSecondaryInCurrentStep();
      if (secTracks && !secTracks->empty()) {
        for (const auto& secTrack : *secTracks) {
          if (secTrack->GetCreatorProcess()) {
            G4int motherPdgCode   = particlePdg;
            G4int daughterPdgCode = secTrack->GetDefinition()->GetPDGEncoding();
            G4ThreeVector mom_se  = secTrack->GetMomentum();
            G4LorentzVector v_se(secTrack->GetPosition(), 0);
            G4LorentzVector p_se(mom_se, std::sqrt(std::pow(particleMass,2)+std::pow(mom_se.mag(),2)));
            gAnaMan.SetSecondaryVertex(daughterPdgCode,motherPdgCode,p_se,v_se);
          }
        }
      }
    }
  }

//=========================== 커스텀 로직 (토글) ===========================
#ifdef E45_SCENARIO_FORCE_2PI
  // ==================== [모드 1: 강제 2π 생성] ====================

  // (A) VP4 평면 miss → 원빔 kill
  {
    if (IsPrimaryBeamLike(theTrack)) {
      static const auto& gGeom = DCGeomMan::GetInstance();
      static const auto& gSize = DetSizeMan::GetInstance();

      static const G4ThreeVector vp4_pos = gGeom.GetGlobalPosition("VP4");
      const double z_vp4 = vp4_pos.z();

      static const G4ThreeVector tgt_c  = gGeom.GetGlobalPosition("SHSTarget");
      static const G4ThreeVector tgt_sz = gSize.GetSize("Target");
      const double Rout    = tgt_sz.y();
      const double half_dz = 0.5 * tgt_sz.z();

      const double z1 = prePoint->GetPosition().z();
      const double z2 = postPoint->GetPosition().z();

      if ( (z1 - z_vp4) * (z2 - z_vp4) <= 0.0 && (z1 != z2) ) {
        const double t = (z_vp4 - z1) / (z2 - z1);
        const G4ThreeVector x1 = prePoint->GetPosition();
        const G4ThreeVector x2 = postPoint->GetPosition();
        const G4ThreeVector x_at = x1 + t*(x2 - x1);

        const G4ThreeVector rel = x_at - tgt_c;
        const double r_xy   = std::hypot(rel.x(), rel.y());
        const bool inside_r = (r_xy <= Rout + 1e-6);
        const bool inside_z = (std::abs(rel.z()) <= half_dz + 1e-6);

        if (!(inside_r && inside_z)) {
          theTrack->SetTrackStatus(fStopAndKill);
          return;
        }
      }
    }
  }

  // (B) 타겟 첫 진입에서 3체 생성 후 원빔 kill
  {
    if (IsPrimaryBeamLike(theTrack)) {
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

          // 채널 선택
          const int mode = GetForce2PiMode();
          int channel = 0;
          if      (mode == 0) channel = 0;
          else if (mode == 1) channel = 1;
          else                channel = (G4UniformRand()<0.5)? 0 : 1;  // 1:1

          auto stackMan = G4EventManager::GetEventManager()->GetStackManager();
          if (stackMan) Spawn3BodyAt(theTrack, channel, vtx, stackMan);

          // (선택) 이벤트 브랜치에 채널 저장
          // gAnaMan.SetChannelID(channel); // 브랜치 추가 후 활성화

          theTrack->SetTrackStatus(fStopAndKill); // 원빔 제거
          return;
        }
      }
    }
  }

  // (C) VP4 downstream 하드 가드 (LH2 바깥에서만)
  {
    if (IsPrimaryBeamLike(theTrack)) {
      static const auto& gGeom = DCGeomMan::GetInstance();
      static const G4ThreeVector vp4_pos = gGeom.GetGlobalPosition("VP4");
      const double z_vp4  = vp4_pos.z();
      const double z_post = postPoint->GetPosition().z();

      const G4Material* postMat = postPoint->GetMaterial();
      const auto postMatName = postMat ? postMat->GetName() : "";

      if (postMatName!="LH2" && z_post > z_vp4 + 0.5*mm) {
        theTrack->SetTrackStatus(fStopAndKill);
        return;
      }
    }
  }

  // (D) VP5 fail-safe
  {
    if (IsPrimaryBeamLike(theTrack)) {
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
  //  - 강제 3체 생성 없음
  //  - VP4 miss-kill / VP4 하드가드 / VP5 fail-safe 없음
  //  - 따라서 원빔(π−, ParentID=0)은 물리 그대로 진행
#endif
//========================= 커스텀 로직 끝 ==========================

  // ---------------- 일반 컷 ----------------
  if(KillStepInIron){
    auto preMaterial = prePoint->GetMaterial();
    if(preMaterial && preMaterial->GetName() == "Iron"){
      theTrack->SetTrackStatus(fStopAndKill);
      return;
    }
  }

  if(particleName == "e-" || particleName == "e+"){
    theTrack->SetTrackStatus(fStopAndKill);
    return;
  }
}

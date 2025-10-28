// -*- C++ -*-
// SteppingAction.cc  (E45)
// 2025-10-26 (user modification):
//   - π+ 원빔( PDG=+211 )도 '원빔처럼' 판정하도록 확대
//   - E72 π+ beam-through용 제너레이터(7227) 지원
//   - 2π 강제 스폰 로직에 π+ 전용 채널 추가
//       * Force2PiMode=0 : π+ π+ n
//       * Force2PiMode=1 : π+ π0 p
//       * Force2PiMode=2 : 위의 1:1 mixture
//
// 기존(π−)과의 대응:
//   * (π−, mode=0): π+ π− n
//   * (π−, mode=1): π− π0 p
//   * (π−, mode=2): 1:1 mixture
//
// - 원빔/자식 출처 태깅: OriginTag 적용 (kPrimaryBeam / kForced2PiChild)
// - (하위호환) SimpleTag("UPSTREAM_HELPER")도 fallback으로 인식
// - 2025-10-24: 길이-가중 리저버 샘플링 도입 (타겟 내부 1회 스폰)

#include "SteppingAction.hh"

#include <string>
#include <cmath>
#include <cstdlib>

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

#include <Randomize.hh>
#include <CLHEP/Units/SystemOfUnits.h>

#include "ConfMan.hh"
#include "PrintHelper.hh"
#include "AnaManager.hh"
#include "PrimaryGeneratorAction.hh"

#include "TrackTag.hh"
#include "DCGeomMan.hh"
#include "DetSizeMan.hh"

#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"

// ==== SCENARIO SWITCH =======================================================
// 1) 강제 2π 생성 모드
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
  //  3) (백업) parentID==0 && pdg==±211 (π±)
  // ──────────────────────────────────────────────────────────────────────────
  inline bool IsPrimaryBeamLike(const G4Track* trk){
    if(!trk) return false;
    if(HasOrigin(trk, EOrigin::kPrimaryBeam)) return true;

    if(auto* s = dynamic_cast<const SimpleTag*>(trk->GetUserInformation())){
      if(s->why_ && std::string(s->why_)=="UPSTREAM_HELPER") return true;
    }

    const auto* pd = trk->GetParticleDefinition();
    const int pdg = pd ? pd->GetPDGEncoding() : 0;
    return (trk->GetParentID()==0 && (pdg==+211 || pdg==-211));
  }

  // ------------------- helper: 3-body spawner (p override) ------------------
  // beamPdg = +211 (π+) or -211 (π−)
  // For π−: mode=0 → π+ π− n,  mode=1 → π− π0 p
  // For π+: mode=0 → π+ π+ n,  mode=1 → π+ π0 p
  void Spawn3BodyAtWithP_ChargedPi(const G4ThreeVector& p_override,
                                   const G4Track* motherTrack, int beamPdg,
                                   int mode01,  // 0/1, mixture는 호출부에서 결정
                                   const G4ThreeVector& vtx, G4StackManager* stackMan,
                                   G4double t_override = -1.0)
  {
    using namespace CLHEP;

    auto pt = G4ParticleTable::GetParticleTable();
    const double mPim = pt->FindParticle("pi-")->GetPDGMass()/GeV;
    const double mPip = pt->FindParticle("pi+")->GetPDGMass()/GeV;
    const double mPi0 = pt->FindParticle("pi0")->GetPDGMass()/GeV;
    const double mN   = pt->FindParticle("neutron")->GetPDGMass()/GeV;
    const double mP   = pt->FindParticle("proton")->GetPDGMass()/GeV;

    const double pGeV = p_override.mag()/GeV;
    // beam hadron(π±) + stationary proton(p) 로 근사 (기존 로직 유지)
    TLorentzVector LVpi ( p_override.x()/GeV, p_override.y()/GeV, p_override.z()/GeV,
                          std::hypot(pGeV, (beamPdg==-211? mPim : mPip)) );
    TLorentzVector LVpro( 0.,0.,0., mP );
    TLorentzVector W = LVpi + LVpro;

    // 출력 입자 셋 채널 정의
    double masses[3];
    struct Out { const char* name; } outs[3];

    if (beamPdg==-211) {
      if (mode01==0) { // π−: π+ π− n
        masses[0]=mPip; masses[1]=mPim; masses[2]=mN;
        outs[0]={"pi+"}; outs[1]={"pi-"}; outs[2]={"neutron"};
      } else {          // π−: π− π0 p
        masses[0]=mPim; masses[1]=mPi0; masses[2]=mP;
        outs[0]={"pi-"}; outs[1]={"pi0"}; outs[2]={"proton"};
      }
    } else { // beamPdg==+211
      if (mode01==0) { // π+: π+ π+ n
        masses[0]=mPip; masses[1]=mPip; masses[2]=mN;
        outs[0]={"pi+"}; outs[1]={"pi+"}; outs[2]={"neutron"};
      } else {          // π+: π+ π0 p
        masses[0]=mPip; masses[1]=mPi0; masses[2]=mP;
        outs[0]={"pi+"}; outs[1]={"pi0"}; outs[2]={"proton"};
      }
    }

    if(W.M() <= masses[0] + masses[1] + masses[2]) return;

    TGenPhaseSpace gen;
    gen.SetDecay(W, 3, masses);
    gen.Generate(); // single unbiased sample

    const double t0 = (t_override>=0.0) ? t_override : motherTrack->GetGlobalTime();

    for(int i=0;i<3;i++){
      const TLorentzVector* d = gen.GetDecay(i);
      auto def = pt->FindParticle(outs[i].name);

      // (A) 자식 트랙 push
      auto dyn = new G4DynamicParticle(def,
                    G4ThreeVector(d->Px()*GeV, d->Py()*GeV, d->Pz()*GeV));
      auto trk = new G4Track(dyn, t0, vtx);
      trk->SetParentID(motherTrack->GetTrackID());
      trk->SetTrackStatus(fAlive);
      trk->SetUserInformation(new OriginTag(EOrigin::kForced2PiChild,
                                            mode01, outs[i].name));
      stackMan->PushOneTrack(trk);

      // (B) SEC 브랜치 기록
      const int  motherPdg   = beamPdg;
      const int  daughterPdg = def->GetPDGEncoding();
      G4LorentzVector p4(d->Px()*GeV, d->Py()*GeV, d->Pz()*GeV, d->E()*GeV);
      G4LorentzVector v4(vtx, 0.0);
      gAnaMan.SetSecondaryVertex(daughterPdg, motherPdg, p4, v4);
    }
  }

} // namespace

// ============================================================================

SteppingAction::SteppingAction()
  : G4UserSteppingAction()
{
  PrintForce2PiModeOnce();
}

SteppingAction::~SteppingAction() {}

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

#ifdef E45_SCENARIO_FORCE_2PI
  // ==================== [모드 1: 강제 2π 생성] ====================

  // (A) VP4 miss → 원빔 kill (원빔처럼 보이는 트랙만)
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

  // (B') 타겟 첫 진입: 리저버 초기화(즉시 spawn 없음)
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
          gAnaMan.MarkTargetTouch(); // 통계

          const auto tid = theTrack->GetTrackID();
          auto &st = m_uat[tid];
          st.active  = true;
          st.hasCand = false;
          st.S       = 0.0;
        }
      }
    }
  }

  // (B'') LH2 내부 스텝마다 후보 버텍스/모멘텀 갱신 (길이-가중 리저버)
  {
    if (IsPrimaryBeamLike(theTrack)) {
      const G4Material* postMat = postPoint->GetMaterial();
      const auto postMatName = postMat ? postMat->GetName() : "";
      const bool inLH2_now = (postMatName=="LH2");

      const auto tid = theTrack->GetTrackID();
      auto it  = m_uat.find(tid);
      if (inLH2_now && it!=m_uat.end() && it->second.active) {
        auto &st = it->second;
        const double ell = theTrack->GetStepLength();
        if (ell > 0.0) {
          st.S += ell;
          if (G4UniformRand() < (ell / st.S)) {
            const double u = G4UniformRand();

            const auto xpre  = prePoint ->GetPosition();
            const auto xpost = postPoint->GetPosition();
            st.vtx = xpre + u*(xpost - xpre);

            const auto ppre  = prePoint ->GetMomentum();
            const auto ppost = postPoint->GetMomentum();
            st.pcand = ppre + u*(ppost - ppre);

            const auto tpre  = prePoint ->GetGlobalTime();
            const auto tpost = postPoint->GetGlobalTime();
            st.tcand = tpre + u*(tpost - tpre);

            st.hasCand = true;
          }
        }
      }
    }
  }

  // (B''') LH2 이탈 시 단 한 번 spawn(후보 지점 p,vtx,t) 하고 원빔 kill
  {
    if (IsPrimaryBeamLike(theTrack)) {
      const G4Material* preMat  = prePoint ->GetMaterial();
      const G4Material* postMat = postPoint->GetMaterial();
      const auto preMatName  = preMat  ? preMat ->GetName() : "";
      const auto postMatName = postMat ? postMat->GetName() : "";
      const bool leaveLH2 =
        (preMatName=="LH2") && (postMatName!="LH2") &&
        (postPoint->GetStepStatus()==fGeomBoundary);

      if (leaveLH2) {
        const auto tid = theTrack->GetTrackID();
        auto it  = m_uat.find(tid);
        if (it!=m_uat.end() && it->second.active) {
          auto &st = it->second;

          if (st.hasCand) {
            // 채널 선택 (0/1 or 50:50)
            int mode = GetForce2PiMode();
            if (mode!=0 && mode!=1) mode = (G4UniformRand()<0.5 ? 0 : 1);

            auto* stackMan = G4EventManager::GetEventManager()->GetStackManager();
            if (stackMan) {
              const int beamPdg = theTrack->GetParticleDefinition()->GetPDGEncoding(); // ±211
              Spawn3BodyAtWithP_ChargedPi(st.pcand, theTrack, beamPdg, mode,
                                          st.vtx, stackMan, st.tcand);
            }

            gAnaMan.MarkForced2Pi(mode);
            theTrack->SetTrackStatus(fStopAndKill); // 원빔 제거
          }
          m_uat.erase(it);
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
  //  - VP4/VP5 가드 없음
  //  - 원빔(π±, ParentID=0)은 물리 그대로 진행
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

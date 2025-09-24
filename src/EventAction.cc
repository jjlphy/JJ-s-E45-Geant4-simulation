// -*- C++ -*-

#include "EventAction.hh"

#include <fstream>

#include <G4RunManager.hh>
#include <G4Event.hh>
#include <G4SDManager.hh>
#include <G4TrajectoryContainer.hh>
#include <G4VVisManager.hh>
#include <G4TrajectoryContainer.hh>
#include <G4Trajectory.hh>
#include <G4VTrajectory.hh>
#include <G4UIterminal.hh>
#include <G4UItcsh.hh>

#include "AnaManager.hh"
#include "TPCPadHelper.hh"
#include "BACSD.hh"
#include "ConfMan.hh"
#include "DetectorConstruction.hh"
#include "FuncName.hh"
#include "BH2SD.hh"
#include "FTOFSD.hh"
#include "HTOFSD.hh"
#include "BACSD.hh"
#include "TPCSD.hh"
#include "TPCPadSD.hh"
#include "TPCEdepSD.hh"
#include "SCHSD.hh"
#include "SDCSD.hh"
#include "TargetSD.hh"
#include "VPSD.hh"
#include "KVCSD.hh"
#include "DetSizeMan.hh"

#include "BVH_UHit.hh"
#include "BVH_DHit.hh"
#include "T0Hit.hh"



#define SDCLASS(x) x ## Hit

namespace
{
  auto& gAnaMan = AnaManager::GetInstance();
  const auto& gConf = ConfMan::GetInstance();
  const auto& gSize = DetSizeMan::GetInstance();
}

//_____________________________________________________________________________
EventAction::EventAction()
  : G4UserEventAction()
{
}

//_____________________________________________________________________________
EventAction::~EventAction()
{
}

//_____________________________________________________________________________
void
EventAction::BeginOfEventAction(const G4Event*)
{
  #ifdef DEBUG
  G4cout << FUNC_NAME << G4endl;
  #endif
  gAnaMan.BeginOfEventAction();
}

//_____________________________________________________________________________
void
EventAction::EndOfEventAction(const G4Event* anEvent)
{
  G4int eventID = anEvent->GetEventID();
  if (eventID % 100 == 0) {
    G4cout << FUNC_NAME << G4endl
	   << "   Event number = " << eventID << G4endl;
  }
  auto SDManager = G4SDManager::GetSDMpointer();
  auto HCTE = anEvent->GetHCofThisEvent();
  if (!HCTE) return;

  {
    static const auto id = SDManager->GetCollectionID("BH2/hit");
    if (id >= 0) {
      auto HC = dynamic_cast<G4THitsCollection<BH2Hit>*>(HCTE->GetHC(id));
      for (G4int i=0, n=HC->entries(); i<n; ++i) {
	gAnaMan.SetHitData((*HC)[i]);
      }
      gAnaMan.SetNhits("BH2", HC->entries());
    }
  }

  {
    static const auto id = SDManager->GetCollectionID("TGT/hit");
    if (id >= 0) {
      auto HC = dynamic_cast<G4THitsCollection<TargetHit>*>(HCTE->GetHC(id));
      for (G4int i=0, n=HC->entries(); i<n; ++i) {
	gAnaMan.SetHitData((*HC)[i]);
      }
      gAnaMan.SetNhits("TGT", HC->entries());
    }
  }

  {
    static const auto id = SDManager->GetCollectionID("HTOF/hit");
    if (id >= 0) {
      auto HC = dynamic_cast<G4THitsCollection<HTOFHit>*>(HCTE->GetHC(id));
      for (G4int i=0, n=HC->entries(); i<n; ++i) {
	gAnaMan.SetHitData((*HC)[i]);
      }
      gAnaMan.SetNhits("HTOF", HC->entries());
    }
  }

  {
    static const auto id = SDManager->GetCollectionID("FTOF/hit");
    if (id >= 0) {
      auto HC = dynamic_cast<G4THitsCollection<FTOFHit>*>(HCTE->GetHC(id));
      for (G4int i=0, n=HC->entries(); i<n; ++i) {
	gAnaMan.SetHitData((*HC)[i]);
      }
      gAnaMan.SetNhits("FTOF", HC->entries());
    }
  }

  {
    static const auto id = SDManager->GetCollectionID("BAC/hit");
    if (id >= 0) {
      auto HC = dynamic_cast<G4THitsCollection<BACHit>*>(HCTE->GetHC(id));
      for(G4int i=0, n=HC->entries(); i<n; ++i){
	gAnaMan.SetHitData((*HC)[i]);
      }
      gAnaMan.SetNhits("BAC", HC->entries());
    }
  }

  {
    static const auto id = SDManager->GetCollectionID("KVC/hit");
    if (id >= 0) {
      auto HC = dynamic_cast<G4THitsCollection<KVCHit>*>(HCTE->GetHC(id));
      for (G4int i=0, n=HC->entries(); i<n; ++i) {
	gAnaMan.SetHitData((*HC)[i]);
      }
      gAnaMan.SetNhits("KVC", HC->entries());
    }
  }
// --- BVH_U Hits ---
{
  static const auto id_bvhu = SDManager->GetCollectionID("BVH_U/hit");
  if (id_bvhu >= 0) {
    auto HC = dynamic_cast<G4THitsCollection<BVH_UHit>*>(HCTE->GetHC(id_bvhu));
    if (HC) {
      for (G4int i=0, n=HC->entries(); i<n; ++i) {
        gAnaMan.SetHitData((*HC)[i]);
      }
      gAnaMan.SetNhits("BVH_U", HC->entries());
    }
  }
}

// --- BVH_D Hits ---
{
  static const auto id_bvhd = SDManager->GetCollectionID("BVH_D/hit");
  if (id_bvhd >= 0) {
    auto HC = dynamic_cast<G4THitsCollection<BVH_DHit>*>(HCTE->GetHC(id_bvhd));
    if (HC) {
      for (G4int i=0, n=HC->entries(); i<n; ++i) {
        gAnaMan.SetHitData((*HC)[i]);
      }
      gAnaMan.SetNhits("BVH_D", HC->entries());
    }
  }
}

// --- T0 Hits ---
{
  static const auto id = SDManager->GetCollectionID("T0/hit");
  if (id >= 0) {
    auto HC = dynamic_cast<G4THitsCollection<T0Hit>*>(HCTE->GetHC(id));
    if (HC) {
      for (G4int i = 0, n = HC->entries(); i < n; ++i) {
        gAnaMan.SetHitData((*HC)[i]);
      }
      gAnaMan.SetNhits("T0", HC->entries());
    }
  }
}

//---SCH Hits---
{
  auto SDManager = G4SDManager::GetSDMpointer();
  static int id_sch = -1;
  if (id_sch < 0) {
    // 과거/현재 표기 모두 시도 (방어코드)
    for (auto key : { "SCH/hit", "/SCH/hit", "SCH/SCHHits", "/SCH/SCHHits" }) {
      int tmp = SDManager->GetCollectionID(key);
      if (tmp >= 0) { id_sch = tmp; break; }
    }
  }
  if (id_sch >= 0) {
    auto HC = dynamic_cast<G4THitsCollection<SCHHit>*>(HCTE->GetHC(id_sch));
    if (HC) {
      for (G4int i=0, n=HC->entries(); i<n; ++i) {
        gAnaMan.SetHitData((*HC)[i]);   // VHitInfo 기반이면 바로 동작
      }
      gAnaMan.SetNhits("SCH", HC->entries());
    } else {
      gAnaMan.SetNhits("SCH", 0);
    }
  }
}

  {
    static const auto id = SDManager->GetCollectionID("VP/hit");
    if (id >= 0) {
      auto HC = dynamic_cast<G4THitsCollection<VPHit>*>(HCTE->GetHC(id));
      for (G4int i=0, n=HC->entries(); i<n; ++i) {
	gAnaMan.SetHitData((*HC)[i]);
      }
      gAnaMan.SetNhits("VP", HC->entries());
    }
  }

  {
    static const auto id = SDManager->GetCollectionID("TPC/hit");
    if(id > 0){
      auto HC = dynamic_cast<G4THitsCollection<TPCHit>*>(HCTE->GetHC(id));
      for (G4int i=0, n=HC->entries(); i<n; ++i) {
	gAnaMan.SetHitData((*HC)[i]);
      }
      gAnaMan.SetNhits("TPC", HC->entries());
    }
  }

  {
    static const G4int id = SDManager->GetCollectionID("TPCPad/hit");
    static const G4int id_edep = SDManager->GetCollectionID("TPCEdep/hit");
    if( id > 0 ){
      //test
      G4int tid_check;
      auto HC = (G4THitsCollection<TPCPadHit>*)(HCTE->GetHC(id));
      auto HC_edep = (G4THitsCollection<TPCEdepHit>*)(HCTE->GetHC(id_edep));

      G4int nhits= HC -> entries();
      G4int nhits_edep = HC_edep -> entries();

      G4int pidtr[MaxHitsTPC]={0};
      G4int ptidtpc[MaxHitsTPC]={0};
      G4int ptidtpc_pid[MaxHitsTPC]={0};
      G4double pmtpc[MaxHitsTPC]={0};
      G4int qqtpc[MaxHitsTPC]={0};
      G4double pxtpc[MaxHitsTPC]={0};
      G4double pytpc[MaxHitsTPC]={0};
      G4double pztpc[MaxHitsTPC]={0};
      G4double pptpc[MaxHitsTPC]={0};
      G4double vtxxtpc[MaxHitsTPC]={0};
      G4double vtxytpc[MaxHitsTPC]={0};
      G4double vtxztpc[MaxHitsTPC]={0};
      G4double vtxpxtpc[MaxHitsTPC]={0};
      G4double vtxpytpc[MaxHitsTPC]={0};
      G4double vtxpztpc[MaxHitsTPC]={0};
      // G4double vtxpptpc[MaxHitsTPC]={0};
      G4double vtxenetpc[MaxHitsTPC]={0};
      G4double detpc[MaxHitsTPC]={0};
      G4int laytpc[MaxHitsTPC]={0};
      G4double lentpc[MaxHitsTPC]={0};
      G4int nparticle=0;
      // G4cout << "TPC  " << nhits << G4endl;


      const G4int pad_configure =  gSize.Get("TpcPadConfigure");
      std::vector<G4ThreeVector> temp_pos;
      std::vector<double> temp_edep;
      std::vector<int> temp_ilay;

      G4int temp_pid;
      G4int temp_tid;
      G4int temp_parentid;
      G4int temp_parentpid;

      std::vector<std::vector<G4ThreeVector>> remain_pos;
      std::vector<std::vector<double>> remain_edep;
      std::vector<std::vector<int>> remain_ilay;
      std::vector<int> remain_pid;
      std::vector<int> remain_tid;
      std::vector<int> remain_parentid;
      std::vector<int> remain_parentpid;


      std::vector<double> cal_posx;
      std::vector<double> cal_posy;
      std::vector<double> cal_posz;
      std::vector<double> cal_edep;

      if(pad_configure==4){
	//Edep Hit sort
	for(int i=0; i<nhits_edep;i++){
	  G4int pid_edep = (*HC_edep)[i]-> GetParticleID();
	  G4int parentid_edep = (*HC_edep)[i]-> GetParentID();
	  G4int parentpid_edep = (*HC_edep)[i]-> GetParentID_pid();
	  G4double edep_edep = (*HC_edep)[i]-> GetEdep();
	  G4int tid_edep = (*HC_edep)[i]-> GetTrackID();
	  G4ThreeVector xyz_edep = (*HC_edep)[i]-> GetPosition();
	  G4int ilay_edep = (*HC_edep)[i]-> GetPadLay();

	  if(i==0){
	    cal_posx.push_back(xyz_edep.x());
	    cal_posy.push_back(xyz_edep.y());
	    cal_posz.push_back(xyz_edep.z());
	    cal_edep.push_back(edep_edep);

	    temp_pid = pid_edep;
	    temp_tid = tid_edep;
	    temp_parentid = parentid_edep;
	    temp_parentpid = parentpid_edep;

	    temp_ilay.push_back(ilay_edep);

	    if(nhits_edep==1){
	      G4ThreeVector ave_pos(CalculateAverage(cal_posx), CalculateAverage(cal_posy), CalculateAverage(cal_posz));
	      temp_pos.push_back(ave_pos);
	      temp_edep.push_back(CalculateSum(cal_edep));

	      //put track info into remain
	      remain_pos.push_back(temp_pos);
	      remain_edep.push_back(temp_edep);
	      remain_ilay.push_back(temp_ilay);
	      remain_pid.push_back(temp_pid);
	      remain_tid.push_back(temp_tid);
	      remain_parentid.push_back(temp_parentid);
	      remain_parentpid.push_back(temp_parentpid);
	    }

	  }

	  else if(i>0){
	    if(tid_edep == temp_tid){
	      if(ilay_edep == temp_ilay.back()){
		cal_posx.push_back(xyz_edep.x());
		cal_posy.push_back(xyz_edep.y());
		cal_posz.push_back(xyz_edep.z());
		cal_edep.push_back(edep_edep);
	      }

	      else if(ilay_edep != temp_ilay.back()){
		//put previous layer info into temp
		G4ThreeVector ave_pos(CalculateAverage(cal_posx), CalculateAverage(cal_posy), CalculateAverage(cal_posz));
		temp_pos.push_back(ave_pos);
		temp_edep.push_back(CalculateSum(cal_edep));

		//clear cal
		cal_posx.clear();
		cal_posy.clear();
		cal_posz.clear();
		cal_edep.clear();

		//put new hit info to cal
		cal_posx.push_back(xyz_edep.x());
		cal_posy.push_back(xyz_edep.y());
		cal_posz.push_back(xyz_edep.z());
		cal_edep.push_back(edep_edep);

		temp_ilay.push_back(ilay_edep);
	      }

	      if(i == nhits_edep-1){
		//End of the edep hits
		//Finalize previous track!!!
		//put current layer info into temp (last layer of previous track)
		G4ThreeVector ave_pos(CalculateAverage(cal_posx), CalculateAverage(cal_posy), CalculateAverage(cal_posz));
		temp_pos.push_back(ave_pos);
		temp_edep.push_back(CalculateSum(cal_edep));

		//put track info into remain
		remain_pos.push_back(temp_pos);
		remain_edep.push_back(temp_edep);
		remain_ilay.push_back(temp_ilay);
		remain_pid.push_back(temp_pid);
		remain_tid.push_back(temp_tid);
		remain_parentid.push_back(temp_parentid);
		remain_parentpid.push_back(temp_parentpid);
	      }

	    }
	    else if(tid_edep != temp_tid){
	      //Finalize previous track!!!
	      //put previous layer info into temp (last layer of previous track)
	      G4ThreeVector ave_pos(CalculateAverage(cal_posx), CalculateAverage(cal_posy), CalculateAverage(cal_posz));
	      temp_pos.push_back(ave_pos);
	      temp_edep.push_back(CalculateSum(cal_edep));

	      //put track info into remain
	      remain_pos.push_back(temp_pos);
	      remain_edep.push_back(temp_edep);
	      remain_ilay.push_back(temp_ilay);
	      remain_pid.push_back(temp_pid);
	      remain_tid.push_back(temp_tid);
	      remain_parentid.push_back(temp_parentid);
	      remain_parentpid.push_back(temp_parentpid);

	      //clear cal

	      cal_posx.clear();
	      cal_posy.clear();
	      cal_posz.clear();
	      cal_edep.clear();

	      //clear temp

	      temp_pos.clear();
	      temp_edep.clear();
	      temp_ilay.clear();

	      temp_pid = -9999;
	      temp_tid = -9999;
	      temp_parentid = -9999;
	      temp_parentpid = -9999;

	      //put new hit of new track to temp and cal
	      cal_posx.push_back(xyz_edep.x());
	      cal_posy.push_back(xyz_edep.y());
	      cal_posz.push_back(xyz_edep.z());
	      cal_edep.push_back(edep_edep);

	      temp_pid = pid_edep;
	      temp_tid = tid_edep;
	      temp_parentid = parentid_edep;
	      temp_parentpid = parentpid_edep;

	      temp_ilay.push_back(ilay_edep);



	      if(i == nhits_edep-1){
		//End of the edep hits
		//Finalize previous track!!!
		//put current layer info into temp (last layer of previous track)
		G4ThreeVector ave_pos(CalculateAverage(cal_posx), CalculateAverage(cal_posy), CalculateAverage(cal_posz));
		temp_pos.push_back(ave_pos);
		temp_edep.push_back(CalculateSum(cal_edep));

		//put track info into remain
		remain_pos.push_back(temp_pos);
		remain_edep.push_back(temp_edep);
		remain_ilay.push_back(temp_ilay);
		remain_pid.push_back(temp_pid);
		remain_tid.push_back(temp_tid);
		remain_parentid.push_back(temp_parentid);
		remain_parentpid.push_back(temp_parentpid);
	      }
	    }
	  }
	}


	for(int n=0;n<remain_tid.size();n++){
	  for(int u=0;u<remain_tid.size();u++){
	    if(remain_tid[n]==remain_parentid[u]){
	      for(int j=0;j<remain_edep[u].size();j++){
		for(int k=0;k<remain_edep[n].size();k++){
		  if(remain_ilay[n][k]==remain_ilay[u][j] && (remain_pos[n][k]-remain_pos[u][j]).mag()<30){
		    remain_edep[n][k]+=remain_edep[u][j];
		  }

		}
	      }
	    }
	  }
	}




      }


      for( G4int i=0; i<nhits; ++i ){
	if(nparticle >20) continue;


	G4ThreeVector vtxpos = (*HC)[i]-> GetVtxPosition();
	G4ThreeVector vtxmom = (*HC)[i]-> GetVtxMomentum();
	G4double vtxene =(*HC)[i]-> GetVtxEnergy();
	G4ThreeVector xyz = (*HC)[i]-> GetPosition();
	G4ThreeVector mom = (*HC)[i]-> GetMomentum();
	G4double tof= (*HC)[i]-> GetTOF();
	G4int tid = (*HC)[i]-> GetTrackID();
	G4int ptid = (*HC)[i]-> GetParentID();
	G4int ptid_pid = (*HC)[i]-> GetParentID_pid();
	G4int pid = (*HC)[i]-> GetParticleID();
	G4double mass = (*HC)[i]-> GetMass();
	G4int charge = (*HC)[i]-> GetCharge();
	// std::cout<<"pid="<<pid<<", mass="<<mass<<", charge="<<charge<<std::endl;
	// getchar();
	G4int ilay = (*HC)[i]-> GetPadLay();
	//      G4double mass = (*HC)[i]-> GetPDGMass(); //mass(GeV)
	G4int parentid = (*HC)[i]-> GetParentID();
	//G4int parentpid = (*HC)[i]-> GetParentID_pid();
	//Get Parent pid
	G4int parentpid = -9999;
	if(parentid>0){
	  const G4Track* parentTrack = nullptr;

	  for (G4int k = 0; k < anEvent->GetNumberOfPrimaryVertex(); k++) {
	    G4PrimaryVertex* primaryVertex = anEvent->GetPrimaryVertex(k);
	    for (G4int j = 0; j < primaryVertex->GetNumberOfParticle(); j++) {
	      G4PrimaryParticle* primaryParticle = primaryVertex->GetPrimary(j);
	      if (primaryParticle->GetTrackID() == parentid) {
		parentpid = primaryParticle->GetPDGcode();
		break;
	      }
	    }
	  }
	}

	G4double tlength = (*HC)[i]-> GettLength();
	G4int irow=(*HC)[i]-> GetPadRow();
	G4double beta = (*HC)[i]-> GetBeta();

	//test edep info save
	G4double edep;


	//test end

	if(pad_configure==3){
	  //const G4int edep_configure =  gSize.Get("TpcEdep");
	  int edep_configure=0;
	  if(edep_configure==0)
	    edep = (*HC)[i]-> GetEdep();

	  else if(edep_configure==1){

	    if(i==0){
	      G4ThreeVector mom_pre = (*HC)[i+1]-> GetMomentum();
	      G4double beta_pre = (*HC)[i+1]-> GetBeta();
	      G4double E_pre = mom_pre.mag() / beta_pre;
	      G4double E_post = mom.mag() / beta;
	      Double_t MeVToeV = TMath::Power(10.,6);
	      edep = (E_post - E_pre) * MeVToeV;
	    }
	    else{
	      G4ThreeVector xyz_pre = (*HC)[i-1]-> GetPosition();
	      if(abs((xyz - xyz_pre).mag())>30){
		G4ThreeVector mom_pre = (*HC)[i+1]-> GetMomentum();
		G4double beta_pre = (*HC)[i+1]-> GetBeta();
		G4double E_pre = mom_pre.mag() / beta_pre;
		G4double E_post = mom.mag() / beta;
		Double_t MeVToeV = TMath::Power(10.,6);
		edep = (E_post - E_pre) * MeVToeV;
	      }
	      else if(tid!=(*HC)[i-1]-> GetTrackID()){
		G4ThreeVector mom_pre = (*HC)[i+1]-> GetMomentum();
		G4double beta_pre = (*HC)[i+1]-> GetBeta();
		G4double E_pre = mom_pre.mag() / beta_pre;
		G4double E_post = mom.mag() / beta;
		Double_t MeVToeV = TMath::Power(10.,6);
		edep = (E_post - E_pre) * MeVToeV;
	      }
	      else{
		G4ThreeVector mom_pre = (*HC)[i-1]-> GetMomentum();
		G4double beta_pre = (*HC)[i-1]-> GetBeta();
		G4double E_pre = mom_pre.mag() / beta_pre;
		G4double E_post = mom.mag() / beta;
		Double_t MeVToeV = TMath::Power(10.,6);
		edep = (E_pre - E_post) * MeVToeV;
		std::cout<<mom_pre.mag()<<" , "<<mom.mag()<<", diff : "<<mom_pre.mag()-mom.mag()<<std::endl;
	      }
	    }
	  }

	}

	G4bool find_track = false;
	G4bool find_hit = false;
	if(pad_configure == 4){
	  for(int j=0;j<remain_tid.size();j++){
	    find_hit = false;
	    if(tid == remain_tid[j]){
	      find_track = true;
	      for(int k=0;k<remain_ilay[j].size();k++){
		if(ilay == remain_ilay[j][k]){
		  //check position difference
		  G4double pos_diff = (remain_pos[j][k] - xyz).mag();

		  if(pos_diff < 50){
		    edep = remain_edep[j][k];
		    find_hit = true;
		    break;
		  }

		}
	      }
	      if(!find_hit){

		std::cout<<"no same hit"<<std::endl;
		std::cout<<"Pos SD ilay : "<<ilay<<std::endl;
	      }
	      break;
	    }
	  }
	  if(!find_track)std::cout<<"no same track"<<std::endl;
	
	  if(!find_hit)edep = -9999;
	  //std::cout<<"Layer : "<<ilay<<", Cal : "<<(*HC)[i]->GetEdep()<<", Exp : "<<edep<<std::endl;


	  int hit_check = 0;
	  if(remain_ilay.size()>0){
	    for(int f=0;f<remain_ilay.size();++f){
	      if(remain_ilay[f].size()>0){
		hit_check+=remain_ilay[f].size();
		for(int k=0;k<remain_ilay[f].size();++k){
		  //std::cout<<"layer : "<<remain_ilay[f][k]<<std::endl;
		}
	      }
	    }
	  }
	  //std::cout<<"nhit Cal : "<<nhits<<"nhit Edep : "<<hit_check<<std::endl;
	}






	G4double slength = (*HC)[i]-> GetsLength();
	//      G4VTrajectoryPoint *tp -> HC->GetPoint(i);
	//    G4int nhits= HC -> entries();
	if( nparticle==0 ){
	  qqtpc[nparticle]=charge;
	  pmtpc[nparticle]=mass;
	  detpc[nparticle]=detpc[nparticle]+edep;
	  // vtxpptpc[nparticle]=sqrt(pow(vtxmom[0],2)+pow(vtxmom[1],2)+pow(vtxmom[2],2));
	  if(ilay>-1){
	    laytpc[nparticle]=laytpc[nparticle]+1;
	  }
	  pidtr[nparticle]=pid;
	  pxtpc[nparticle]=mom[0];
	  pytpc[nparticle]=mom[1];
	  pztpc[nparticle]=mom[2];


	  //////////////////////vertex information /////////////////////////
	  vtxpxtpc[nparticle]=vtxmom[0];
	  vtxpytpc[nparticle]=vtxmom[1];
	  vtxpztpc[nparticle]=vtxmom[2];
	  // vtxpptpc[nparticle]=sqrt(pow(vtxmom[0],2)+pow(vtxmom[1],2)+pow(vtxmom[2],2));

	  vtxxtpc[nparticle]=vtxpos[0];
	  vtxytpc[nparticle]=vtxpos[1];
	  vtxztpc[nparticle]=vtxpos[2];

	  vtxenetpc[nparticle]=vtxene;

	  pptpc[nparticle]=sqrt(pow(mom[0],2.)+pow(mom[1],2.)+pow(mom[2],2.));
	  lentpc[nparticle]=tlength;
	  ptidtpc[nparticle]=ptid;
	  ptidtpc_pid[nparticle]=ptid_pid;
	  nparticle=nparticle+1;

	}else if( nparticle>0 ){
	  //G4cout<<nparticle<<G4endl;
	}

	if( (pidtr[nparticle-1] != pid) || (pidtr[nparticle-1] == pid && vtxpxtpc[nparticle-1] != vtxmom[0] && vtxpytpc[nparticle-1] != vtxmom[1] && vtxpztpc[nparticle-1] != vtxmom[2])){
	  qqtpc[nparticle]=charge;
	  pmtpc[nparticle]=mass;
	  detpc[nparticle]=detpc[nparticle]+edep;
	  if(ilay>-1){
	    laytpc[nparticle]=laytpc[nparticle]+1;
	  }
	  pidtr[nparticle]=pid;
	  pxtpc[nparticle]=mom[0];
	  pytpc[nparticle]=mom[1];
	  pztpc[nparticle]=mom[2];

	  vtxpxtpc[nparticle]=vtxmom[0];
	  vtxpytpc[nparticle]=vtxmom[1];
	  vtxpztpc[nparticle]=vtxmom[2];
	  // vtxpptpc[nparticle]=sqrt(pow(vtxmom[0],2)+pow(vtxmom[1],2)+pow(vtxmom[2],2));
	  vtxenetpc[nparticle]=vtxene;

	  vtxxtpc[nparticle]=vtxpos[0];
	  vtxytpc[nparticle]=vtxpos[1];
	  vtxztpc[nparticle]=vtxpos[2];

	  pptpc[nparticle]=sqrt(pow(mom[0],2.)+pow(mom[1],2.)+pow(mom[2],2.));
	  lentpc[nparticle]=tlength;
	  ptidtpc[nparticle]=ptid;
	  ptidtpc_pid[nparticle]=ptid_pid;
	  nparticle=nparticle+1;
	}else if (pidtr[nparticle-1] == pid && vtxpxtpc[nparticle-1] == vtxmom[0] && vtxpytpc[nparticle-1] == vtxmom[1] && vtxpztpc[nparticle-1] == vtxmom[2]){
	  if( ptidtpc[nparticle-1] != ptid){
	    qqtpc[nparticle]=charge;
	    pmtpc[nparticle]=mass;
	    detpc[nparticle]=detpc[nparticle]+edep;
	    if(ilay>-1){
	      laytpc[nparticle]=laytpc[nparticle]+1;
	    }
	    pidtr[nparticle]=pid;
	    pxtpc[nparticle]=mom[0];
	    pytpc[nparticle]=mom[1];
	    pztpc[nparticle]=mom[2];

	    vtxpxtpc[nparticle]=vtxmom[0];
	    vtxpytpc[nparticle]=vtxmom[1];
	    vtxpztpc[nparticle]=vtxmom[2];
	    // vtxpptpc[nparticle]=sqrt(pow(vtxmom[0],2)+pow(vtxmom[1],2)+pow(vtxmom[2],2));

	    vtxenetpc[nparticle]=vtxene;

	    vtxxtpc[nparticle]=vtxpos[0];
	    vtxytpc[nparticle]=vtxpos[1];
	    vtxztpc[nparticle]=vtxpos[2];

	    pptpc[nparticle]=sqrt(pow(mom[0],2.)+pow(mom[1],2.)+pow(mom[2],2.));
	    lentpc[nparticle]=tlength;
	    ptidtpc[nparticle]=ptid;
	    ptidtpc_pid[nparticle]=ptid_pid;
	    nparticle=nparticle+1;
	  }
	  else{
	    detpc[nparticle-1]=detpc[nparticle-1]+edep;
	    if(ilay>-1){
	      laytpc[nparticle-1]=laytpc[nparticle-1]+1;
	    }
	    lentpc[nparticle-1]=tlength;
	  }

	}

	if(ilay>-1){ //-->  -1 : TPC, layer is from 0 to 38. 2012.10.30
	  //Dead Channel Delete
	  if(TPCPadHelper::GetDeadCon((*HC)[i]->GetPadLay(), (*HC)[i]->GetPadRow())){
	    continue;
	  }
	  int restype = gConf.Get<G4int>("ResType");
	  switch(restype){
	  case 0:
	    gAnaMan.SetCounterDataSimple( nparticle-1,tof, xyz, mom, tid, pid, ilay,
					  irow, beta, edep/CLHEP::MeV, parentid, parentpid, tlength,slength );
	    break;
	  case 1:
	    gAnaMan.SetCounterDataExp( nparticle-1,tof, xyz, mom, tid, pid, ilay,
				       irow, beta, edep/CLHEP::MeV, parentid, parentpid, tlength,slength );
	    break;
	  default:
	    gAnaMan.SetCounterDataSimple( nparticle-1,tof, xyz, mom, tid, pid, ilay,
					  irow, beta, edep/CLHEP::MeV, parentid, parentpid, tlength,slength );
	    G4cout<<"TPC Resolution type is not determined. Now using constant resolution"<<G4endl;
	    break;
	  }

	}
      }

      for(G4int i=0;i<nparticle;i++){
	gAnaMan.SetTPCData(i, pidtr[i], ptidtpc[i], ptidtpc_pid[i], pxtpc[i],pytpc[i],pztpc[i],pptpc[i], qqtpc[i], pmtpc[i], detpc[i], lentpc[i],laytpc[i],
			   vtxpxtpc[i],vtxpytpc[i],vtxpztpc[i],
			   vtxxtpc[i],vtxytpc[i],vtxztpc[i], vtxenetpc[i]);
      }


    }



  }


  #if 0
  auto trajectoryContainer = anEvent->GetTrajectoryContainer();
  if(trajectoryContainer && G4VVisManager::GetConcreteInstance()){
    G4int n_trajectories = trajectoryContainer->entries();
    for(G4int i=0; i<n_trajectories; ++i){
      auto trj = (G4Trajectory*)((*(anEvent->GetTrajectoryContainer()))[i]);
      trj->DrawTrajectory();
    }
  }
  #endif

  gAnaMan.EndOfEventAction();
}

G4double
EventAction::CalculateAverage(const std::vector<double>& pos)
{
  if(pos.size()==0)return -1000;

  G4double pos_sum = 0;
  for(size_t i = 0; i < pos.size(); ++i){
    pos_sum+=pos[i];
  }

  return pos_sum / pos.size();
}

G4double
EventAction::CalculateSum(const std::vector<double>& edep)
{
  G4double edep_sum = 0;

  for(size_t i = 0; i < edep.size(); ++i){
    edep_sum+=edep[i];
  }

  return edep_sum;
}

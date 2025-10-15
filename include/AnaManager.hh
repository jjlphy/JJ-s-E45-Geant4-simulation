// -*- C++ -*-

#ifndef ANA_MANAGER_HH
#define ANA_MANAGER_HH

#include <vector>
#include <unordered_map>
#include <set>
#include <algorithm>

#include <G4LorentzVector.hh>
#include <G4ThreeVector.hh>

#include <TParticle.h>
#include <TVector3.h>

class G4ParticleDefinition;

class TFile;
class TTree;

struct Track;

class VHitInfo;

//_____________________________________________________________________________
static const G4int MaxHits    = 500;
static const G4int MaxHitsTPC = 500;
static const G4int MaxPrimaryParticle = 10;

//const int MaxTrack = 1560*4;
//const int MaxTrack = 78*4;
const G4int MaxTrack = 54*20;

const G4int MaxNthLay = 40;
const G4int MaxNthPad = 250;


void initTrack(Track* aTrack);
void initTrack_ku(Track* aTrack);
int setInitialPara(Track* aTrack, double* initPara);
int setVirtualPlane(Track* aTrack);
void minuitInit(double printLevel);

static const int MAXtpctrNum=30;
static const int MAXtpctrhitNum=500;

//_____________________________________________________________________________
struct CounterData
{
  G4int ntrk;
  G4double resoX;
  G4int trackID;
  G4int particleID;
  G4double time;
  G4double beta;
  G4double edep;
  G4double dedx;
  G4double slength;
  G4double tlength;
  G4double mass;
  G4double pos0[3];
  G4double pos[3];
  G4double mom[4];
  G4int iLay;
  G4int iPad;
  G4int iRow;
  G4int parentID; //parent particle track id
  G4int parentPID;
  G4double res[3];
};

//_____________________________________________________________________________
struct TPCData
{
  G4int tpctr;
  G4int tpcpid;
  G4int tpcparentid;
  G4int tpcparentid_pid;
  G4double tpcpx;
  G4double tpcpy;
  G4double tpcpz;
  G4double tpcpp;

  G4double tpcpxfit;
  G4double tpcpyfit;
  G4double tpcpzfit;
  G4double tpcppfit;
  G4double tpcptfit;

  G4int tpcqq;
  G4double tpcpm;
  G4double tpcde;
  G4double tpclen;
  G4double tpcdedx;
  G4int tpclay;
  G4double tpcvtxpx;
  G4double tpcvtxpy;
  G4double tpcvtxpz;
  G4double tpcvtxx;
  G4double tpcvtxy;
  G4double tpcvtxz;
  //  G4double tpcene2
};

//_____________________________________________________________________________
struct Event
{
  Int_t evnum; // Event number
  Int_t generator;
  Int_t mode;        // mode number
  Int_t inc;        // INC id number

  Int_t HitNum_K;
  //  int tpctrNum_K;

  Int_t HitNum_p;
  //  int tpctrNum_p;

  std::map<TString, std::vector<TParticle>> hits;

  /* number of ntrks in TPC by shhwang*/
  Int_t ntrtpc;
  Double_t trpptpc[MaxHitsTPC];
  Double_t trpxtpc[MaxHitsTPC];
  Double_t trpytpc[MaxHitsTPC];
  Double_t trpztpc[MaxHitsTPC];
  Double_t trpttpc[MaxHitsTPC];

  Double_t trpptpcfit[MaxHitsTPC];
  Double_t trpxtpcfit[MaxHitsTPC];
  Double_t trpytpcfit[MaxHitsTPC];
  Double_t trpztpcfit[MaxHitsTPC];
  Double_t trpttpcfit[MaxHitsTPC];

  Int_t trqqtpc[MaxHitsTPC];
  Int_t trpidtpc[MaxHitsTPC];
  Int_t trparentidtpc[MaxHitsTPC];
  Int_t trparentid_pid_tpc[MaxHitsTPC];
  Double_t trpmtpc[MaxHitsTPC];
  Double_t trdetpc[MaxHitsTPC];
  Double_t trlentpc[MaxHitsTPC];
  Double_t trdedxtpc[MaxHitsTPC];
  Double_t trdedxtrtpc[MaxHitsTPC]; //trancated mean, but now just mean

  Int_t trlaytpc[MaxHitsTPC];

  Double_t vtpxtpc[MaxHitsTPC];
  Double_t vtpytpc[MaxHitsTPC];
  Double_t vtpztpc[MaxHitsTPC];
  Double_t vtpptpc[MaxHitsTPC];

  Double_t vtxtpc[MaxHitsTPC];
  Double_t vtytpc[MaxHitsTPC];
  Double_t vtztpc[MaxHitsTPC];

  Double_t vtxtpcfit[MaxHitsTPC];
  Double_t vtytpcfit[MaxHitsTPC];
  Double_t vtztpcfit[MaxHitsTPC];

  /////PAD multiplicity & ASAD multiplicy
  Int_t nthlay[MaxTrack];
  Int_t nthpad[MaxTrack];
  Int_t laypad[MaxTrack][MaxNthLay][MaxNthPad]; //[layer][pad number]


  ///////////////
  Int_t nhittpc;                 // Number of Hit in Pads
  Int_t ntrk[MaxTrack];        // Number of Track

  Int_t ititpc[MaxTrack];      // Track ID
  Int_t idtpc[MaxTrack];       // Particle ID
  Double_t xtpc[MaxTrack];     // coordinates
  Double_t ytpc[MaxTrack];     // coordinates
  Double_t ztpc[MaxTrack];     // coordinates

  Double_t xtpc_pad[MaxTrack];     // coordinates
  Double_t ytpc_pad[MaxTrack];     // coordinates
  Double_t ztpc_pad[MaxTrack];     // coordinates

  Double_t dxtpc_pad[MaxTrack];     // coordinates
  Double_t dytpc_pad[MaxTrack];     // coordinates
  Double_t dztpc_pad[MaxTrack];     // coordinates

  Double_t x0tpc[MaxTrack];    // coordinates
  Double_t y0tpc[MaxTrack];    // coordinates
  Double_t z0tpc[MaxTrack];    // coordinates

  Double_t resoX[MaxTrack];    // coordinates
  Double_t resxtpc[MaxTrack];  // coordinates 
  Double_t resytpc[MaxTrack];  // coordinates 
  Double_t resztpc[MaxTrack];  // coordinates 


  Double_t pxtpc[MaxTrack];    // momentum
  Double_t pytpc[MaxTrack];    // momentum
  Double_t pztpc[MaxTrack];    // momentum
  Double_t pptpc[MaxTrack];    // momentum
  Double_t masstpc[MaxTrack];    // mass


  Double_t timetpc[MaxTrack];    // global time
  Double_t tlengthtpc[MaxTrack];    // global time

  Double_t betatpc[MaxTrack];    // beta

  Double_t edeptpc[MaxTrack];    // Energy deposit
  Double_t dedxtpc[MaxTrack];    // Energy deposit/dx
  Double_t slengthtpc[MaxTrack];    // Energy deposit/dx

  Int_t iPadtpc[MaxTrack];      // number of pad
  Int_t laytpc[MaxTrack];      // number of pad layer
  Int_t rowtpc[MaxTrack];      // number of pad raw
  Double_t toftpc[MaxTrack];   // tof
  Int_t parentID[MaxTrack];      // parent id
  Int_t parentPID[MaxTrack];      // parent id
  Double_t cir_r[MaxTrack];   // fit radius
  Double_t cir_x[MaxTrack];   // fit center x
  Double_t cir_z[MaxTrack];   // fit center z
  Double_t cir_fit[MaxTrack];   // fit center fit
  Int_t vtx_flag[MaxTrack]; // flag, how to estimate vtx
  Double_t a_fory[MaxTrack]; // co-efficient a for linear track (y, theta)
  Double_t b_fory[MaxTrack]; // co-efficient b for linear track (y, theta)
};

//_____________________________________________________________________________
class AnaManager
{
public:
  static G4String ClassName();
  static AnaManager& GetInstance();
  ~AnaManager();

private:
  AnaManager();
  AnaManager(const AnaManager&);
  AnaManager& operator=(const AnaManager&);

private:
  TFile* m_file;
  TTree* m_tree;
  TTree* m_tree_light;
  G4int m_on_off_helm;
  G4int m_pad_config;
  G4int m_experiment;

  G4double m_effective_thickness;
  G4double m_mom_kaon_lab;
  G4double m_cos_theta;
  G4double m_cos_theta_lambda;
  
  CounterData counterData[MaxTrack];
  TPCData tpcData[MAXtpctrNum];

  int HitNum;
  int tpctrNum;
  int HitNum_K;
  int HitNum_p;

  G4double mean[MAXtpctrNum];//read fit parameters
  G4double trmean[MAXtpctrNum];//read fit parameters
  G4double cir_r[MAXtpctrNum];//read fit parameters
  G4double error[MAXtpctrNum];//read fit parameters
  G4double chi2[MAXtpctrNum];//read fit parameters
  G4double ndf[MAXtpctrNum];//read fit parameters
  G4double Pz[MAXtpctrNum];//read fit parameters

  ////////////////////getenv parameters
  G4double pad_length_in;
  G4double pad_length_out;
  G4double pad_gap;
  G4double pad_in_width;
  G4double pad_out_width;
  G4double pad_in_num;
  G4double pad_out_num;
  G4double truncated_mean_cut;

  G4double angle[40];
  G4double seg_angle[40];
  G4double seg_width[40];
  G4int numpads[40];
  G4double pad_in[40];
  G4double pad_out[40];
  G4double tpc_rad;

  // --------------------------------
  // combine beam and reaction generator  
  G4bool m_do_hit_tgt;
  G4bool m_do_generate_beam;
  G4bool m_do_combine;
  G4bool m_threshold_con;
  G4int m_effective_evnum;
  G4int m_next_generator;
  G4int m_first_generator;
  G4int m_second_generator;
  G4ThreeVector m_next_pos;
  G4ThreeVector m_next_mom;
  G4ThreeVector m_vertex_pos;
  G4ThreeVector m_debug_pos;
  // --------------------------------
  // === NEW (jaejin 2025-10-15): minimal event-level truth flags for 2pi study ===
  int  m_forced2pi_flag    = 0;   // 0: no forced-2pi, 1: forced-2pi occurred
  int  m_forced2pi_channel = -1;  // -1: n/a, 0: (pi+ pi- n), 1: (pi- pi0 p)
  int  m_tgt_touch_flag    = 0;   // 0: beam did NOT enter LH2, 1: entered LH2
  // ==============================================================================


  // --------------------------------
  // for checking decay particle
  std::pair<G4String, G4String> m_previous_particle; // particle name, process name
  const std::unordered_map<G4int, G4String> m_focus_particle = {
    // generator id, particle name
    {7201, "kaon-"},
    {7202, "lambda"},
    {7203, "lambda"},
    {7204, "sigma-"},
    {7205, "lambda"},
    {7206, "sigma+"},
    {7207, "kaon-"},
    {7208, "kaon0S"}
  };
  G4int m_decay_particle_code;
  G4ThreeVector m_decay_position;
  // --------------------------------

  // --------------------------------
  // for acceptance study
  const G4double m_edep_threshold = 0.2; // MeV
  const std::unordered_map<G4int, std::vector<G4int>> m_tpc_check_list = {
  //  gen   { check parentID, PDG codes of check list }
    { 7202, {1, 2212, -211} }, // eta Lambda
    { 7203, {1, 2212, -211} }, // pi0 Lambda
    { 7204, {0, +211, -211} }, // pi+ Sigma-
    { 7205, {1, 2212, -211} }, // pi0 Sigma0
    { 7206, {0, -211, +211, 2212} }, // pi- Sigma+
    { 7207, {0, -321, 2212} }, // K p
    { 7208, {1, +211, -211} }  // k0 n
  };
  const std::vector<G4int> m_forward_seg_narrow{16, 17, 18, 19, 20, 21, 22};
  const std::vector<G4int> m_forward_seg_wide{10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29};
  const std::vector<G4int> m_forward_seg_all{6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33};
  const G4double m_refractive_index_kvc = 1.46;

  G4int m_trig_flag_int;
  G4int m_focus_parent_id;  
  G4bool m_kaon_beam_flag;
  // --------------------------------

  
public:
  void BeginOfRunAction(G4int runnum);
  void EndOfRunAction();
  void BeginOfEventAction();
  int  EndOfEventAction();
  void MakeBranch(const G4String& sd_name);
  void MakeHistogram(const G4String& sd_name);
  void SetNhits(const G4String& sd_name, G4int nhits);
  void SetHitData(const VHitInfo* hit);
  void SetTPCData(G4int tpctr, G4int tpcpid, G4int tpcparentid,
                  G4int tpcparentid_pid, G4double tpcpx, G4double tpcpy,
                  G4double tpcpz,G4double tpcpp,  G4int tpcqq, G4double tpcpm,
                  G4double tpcde, G4double tpclen, G4int tpclay,
                  G4double vtxpxtpc2,G4double vtxpytpc2,G4double vtxpztpc2,
                  G4double vtxxtpc2,G4double vtxytpc2,G4double vtxztpc2,
                  G4double vtxenetpc2);
  void SetCounterDataSimple(G4int ntrk, G4double time, G4ThreeVector pos,
                      G4ThreeVector mom, G4int track, G4int particle,
                      G4int iLay, G4int iRow, G4double beta, G4double edep,
			    G4int parentid, G4int parentpid, G4double tlength, G4double slength);
  void SetCounterDataExp(G4int ntrk, G4double time, G4ThreeVector pos,
                      G4ThreeVector mom, G4int track, G4int particle,
                      G4int iLay, G4int iRow, G4double beta, G4double edep,
			 G4int parentid, G4int parentpid, G4double tlength, G4double slength);
  void SetFermiMomentum(const G4ThreeVector& p);
  void SetGeneratorID(G4int generator);
  void SetModeID(G4int mode);
  void SetIncID(G4int inc);
  void SetPrimaryParticle(G4int id, G4int pdg,
                          const G4LorentzVector& p,
                          const G4LorentzVector& v,
                          G4bool is_virtual_beam=false);
  void SetSecondaryVertex(G4int pdg, G4int motherPdg,
			  const G4LorentzVector& p,
			   const G4LorentzVector& v);
  void SetBeamInfo(G4int pdg,
		   const G4LorentzVector& p,
		   const G4LorentzVector& v);
  void SetPrimaryVertex(G4int id, const G4ThreeVector& x);
  void SetPrimaryVertex(G4int id, G4double x, G4double y, G4double z);
  void SetEffectiveThickness(G4double effective_thickness);
  G4double GetEffectiveThickness();
  void SetMomKaonLab(G4double mom_kaon_lab);
  void SetCosTheta(G4double cos_theta);
  void SetCosThetaLambda(G4double cos_theta_lambda);
  void SetPreviousParticle(G4String particle_name, G4String process_name);
  
  // --------------------------------
  // combine beam and reaction generator
  void   SetDoHitTGT(G4bool do_hit_tgt);
  G4bool GetDoHitTGT();
  void   SetDoGenerateBeam(G4bool do_generate_beam);
  G4bool GetDoGenerateBeam();
  void   SetDoCombine(G4bool do_combine);
  G4bool GetDoCombine();
  void   SetThresholdCondition(G4bool threshold_con);
  G4bool GetThresholdCondition();
  void  SetEffectiveEvnum(G4int effective_evnum);
  G4int GetEffectiveEvnum();
  void  SetNextGenerator(G4int next_generator);
  G4int GetNextGenerator();
  void  SetFirstGenerator(G4int first_generator);
  G4int GetFirstGenerator();
  void  SetSecondGenerator(G4int second_generator);
  G4int GetSecondGenerator();
  void          SetNextPos(G4double vx, G4double vy, G4double vz);
  G4ThreeVector GetNextPos();
  void          SetNextMom(G4double px, G4double py, G4double pz);
  G4ThreeVector GetNextMom();
  void          SetVertexPos(G4double vx, G4double vy, G4double vz);
  G4ThreeVector GetVertexPos();
  void          SetDebugPos(G4double vx, G4double vy, G4double vz);
  G4ThreeVector GetDebugPos();
  // --------------------------------

  // --------------------------------
  // for checking decay particle
  std::pair<G4String, G4String> GetPreviousParticle();
  G4String GetFocusParticle(G4int generator_id);
  void  SetDecayParticleCode(G4int decay_particle_code);
  G4int GetDecayParticleCode();
  void SetDecayPosition(G4ThreeVector decay_position);
  G4ThreeVector GetDecayPosition();
  G4bool IsInsideHtof(G4ThreeVector position);
  // --------------------------------

  // --------------------------------
  // for trigger
  void SetFocusParentID(G4int focus_parent_id);
  // --------------------------------

  // === NEW (jaejin 2025-10-15): event-level truth flags for 2pi study ===
  // setters (called from stepping)
  void ClearForced2Pi();              // reset 2pi flags to (0, -1)
  void MarkForced2Pi(int channel);    // set (flag=1, channel=0/1)
  void MarkTargetTouch();             // set tgt_touch_flag=1

  // getters (used from analysis or other code)
  int  GetForced2PiFlag() const { return m_forced2pi_flag; }
  int  GetForced2PiChannel() const { return m_forced2pi_channel; }
  int  GetTargetTouchFlag() const { return m_tgt_touch_flag; }
  // =========================================================================

  
  int CircleIntersect(double x1, double y1, double r1, double x2, double y2, double r2,
		      double ca1, double cb1, double ct01, int qq1,
		      double ca2, double cb2, double ct02, int qq2,
		      double inter1[3], double inter2[3])
  {
    // function inputs: x1, y1, r1, x2, y2, r2
    // function output: inter1, inter2 = coordinates of intersections

    double d,e,f,g,a,b,c;
    double x,y,discrim;

    if(x1 == x2 && y1 == y2){
      G4cout << x1 << " " << y1 << " " << r1 << G4endl;
      return 0;
    }
    //    G4cout << x1 << " " << y1 << " " << r1 << G4endl;
    //    G4cout << x2 << " " << y2 << " " << r2 << G4endl;

    d = -0.5*(r1*r1 - r2*r2 - x1*x1 - y1*y1 + x2*x2 + y2*y2);
    e =  0.5*(r1*r1 + r2*r2 - x1*x1 - y1*y1 - x2*x2 - y2*y2);
    if(fabs(y1-y2) < 1.0e-20) {
      x = d/(x1 - x2);
      a = 1.0;
      b = 0.0;
      c = x*x - x*(x1 + x2) - e;
      discrim = -4*a*c;
      if(discrim < 0) return 0;
      y = sqrt(discrim) / (2*a);
      inter1[0] = x;
      inter1[1] = y;
      inter2[0] = x;
      inter2[1] = -y;
      return 1;
    }
    f = (x1 - x2) / (y1 - y2);
    g = d / (y1 - y2);
    // cout << "d=" << d << " e=" << e << " f=" << f << " g=" << g << endl;

    a = 1. + f*f;
    b = f*(y1 + y2) - 2*f*g - (x1 + x2);
    c = g*g - g*(y1 + y2) - e;
    //    G4cout << "a=" << a << " b=" << b << " c=" << c << G4endl;

    discrim = b*b - 4*a*c;
    //    G4cout << "discrim = " << discrim << G4endl;
    if(discrim < 0) return 0;
    inter1[0] = (-b + sqrt(discrim)) / (2*a);
    inter1[1] = g - f*inter1[0];
    inter2[0] = (-b - sqrt(discrim)) / (2*a);
    inter2[1] = g - f*inter2[0];


    double theta11 = atan2(inter1[1]-y1, inter1[0]-x1);
    double theta21 = atan2(inter1[1]-y2, inter1[0]-x2);
    double theta12 = atan2(inter2[1]-y1, inter2[0]-x1);
    double theta22 = atan2(inter2[1]-y2, inter2[0]-x2);

    double tmp_y11 = -1.*(double)qq1*ca1*r1*(theta11-ct01)+cb1;
    double tmp_y21 = -1.*(double)qq2*ca2*r2*(theta21-ct02)+cb2;
    double tmp_y12 = -1.*(double)qq1*ca1*r1*(theta12-ct01)+cb1;
    double tmp_y22 = -1.*(double)qq2*ca2*r2*(theta22-ct02)+cb2;

    // std::cout<<"theta11="<<theta11<<", theta21="<<theta21
    // 	     <<", theta12="<<theta12<<", theta22="<<theta22<<std::endl;

    // std::cout<<"tmp_y11="<<tmp_y11<<", tmp_y21="<<tmp_y21<<std::endl;
    // std::cout<<"tmp_y12="<<tmp_y12<<", tmp_y22="<<tmp_y22<<std::endl;
    //getchar();

    inter1[2] = (tmp_y11+tmp_y21)/2.;
    inter2[2] = (tmp_y12+tmp_y22)/2.;


    return 1;
  }


  double linearFitter(const int np,
		      const double *x,
		      const double *y, double *er,
		      double *a, double *b){

    int i;

    double alpha=0.;
    double beta=0.;
    double gamma=0.;
    double AA=0;
    double BB=0;

    for(i=0;i<np;i++){
      alpha+=x[i]*x[i]/er[i]/er[i];
      beta+=x[i]/er[i]/er[i];
      gamma+=1./er[i]/er[i];
      AA+=y[i]*x[i]/er[i]/er[i];
      BB+=y[i]/er[i]/er[i];

      //  G4cout<<"x test: "<<x[i]<<G4endl;
      //  G4cout<<"y test: "<<y[i]<<G4endl;
    }

    //  G4cout<<"beta test: "<<beta<<G4endl;
    //  G4cout<<"alpha test: "<<alpha<<G4endl;
    //  G4cout<<"gamma test: "<<gamma<<G4endl;

    *a=(gamma*AA -  beta*BB)/(alpha*gamma-beta*beta);
    *b=(-beta *AA + alpha*BB)/(alpha*gamma-beta*beta);

    //  G4cout<<"a test: "<<(*a)<<G4endl;
    //  G4cout<<"b test: "<<(*b)<<G4endl;

    return 1.;
  }
};

//_____________________________________________________________________________
inline G4String
AnaManager::ClassName()
{
  static G4String s_name("AnaManager");
  return s_name;
}

//_____________________________________________________________________________
inline AnaManager&
AnaManager::GetInstance()
{
  static AnaManager s_instance;
  return s_instance;
}

#endif

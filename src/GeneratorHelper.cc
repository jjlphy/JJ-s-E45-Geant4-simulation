// -*- C++ -*-

/**
 * for the EventGeneration for the E27 experiment
 */

#include "GeneratorHelper.hh"

#include <cmath>
#include <iomanip>

#include <CLHEP/Units/PhysicalConstants.h>
#include <Randomize.hh>
#include <G4LorentzVector.hh>

#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TGenPhaseSpace.h>

#include "AngDisGenerator.hh"

//_____________________________________________________________________________
G4ThreeVector
UniformDirectionInUV(double u0, double v0, double hu, double hv)
{
  double du=0.,dv=0.;
  if(hu!=0.0) du=(G4UniformRand()-0.5)*hu;
  if(hv!=0.0) dv=(G4UniformRand()-0.5)*hv;

  double u=u0+du, v=v0+dv;
  double ninv=1./sqrt(1.+u*u+v*v);

  return G4ThreeVector(u*ninv, v*ninv, ninv);
}

//_____________________________________________________________________________
G4ThreeVector
GaussDirectionInUV(double u0, double v0, double su, double sv)
{
  double du=0.,dv=0.;
  if(su!=0.0) du=G4RandGauss::shoot(0.0,su);
  if(sv!=0.0) dv=G4RandGauss::shoot(0.0,sv);

  double u=u0+du, v=v0+dv;
  double ninv=1./sqrt(1.+u*u+v*v);

  return G4ThreeVector(u*ninv, v*ninv, ninv);
}

//_____________________________________________________________________________
G4ThreeVector
UniformPosition(double hx, double hy, double hz)
{
  double x=0., y=0., z=0.;
  if(hx!=0.0) x+=(G4UniformRand()-0.5)*hx;
  if(hy!=0.0) y+=(G4UniformRand()-0.5)*hy;
  if(hz!=0.0) z+=(G4UniformRand()-0.5)*hz;

  return G4ThreeVector(x, y, z);
}

//_____________________________________________________________________________
G4ThreeVector
GaussPosition(double sx, double sy, double hz)
{
  double x=0., y=0., z=0.;
  if(sx!=0.0) x+=G4RandGauss::shoot(0.0,sx);;
  if(sy!=0.0) y+=G4RandGauss::shoot(0.0,sy);;
  if(hz!=0.0) z+=(G4UniformRand()-0.5)*hz;

  return G4ThreeVector(x, y, z);
}

//_____________________________________________________________________________
G4ThreeVector
GaussPosition_LqTarg(double x0, double y0, double z0,
                     double dx, double dy, double targ_r, double targ_height)
{
  double x=0., y=0., z=0.;
  int ntry=0;
  while(1){
    if(dx!=0.0) x+=G4RandGauss::shoot(0.,dx);;
    if(dy!=0.0) y+=G4RandGauss::shoot(0.,dy);;
    if(targ_r!=0.0) z+=(G4UniformRand()-0.5)*2.*targ_r;
    x+=x0;
    y+=y0;
    if((x*x+z*z) < targ_r*targ_r && fabs(y) < targ_height)
      break;
    x=0.;
    y=0.;
    z=0.;
    if(ntry>1000)
      G4Exception("GeneratorHelper::GaussPosition_LqTarg",
		  "Wrong target input",
		  RunMustBeAborted,
		  "GeneratorHelper::Wrong target input!!");
    else
      ++ntry;
  }
  z+=z0;
  return G4ThreeVector(x, y, z);
}

//_____________________________________________________________________________
double
BreitWigner(double mean, double gamma)
{
  double r=G4UniformRand()-0.5;
  return mean+0.5*gamma*tan(acos(-1.)*r);
}

//_____________________________________________________________________________
bool
Decay2Body(double Mini, double Mf1, double Mf2,
           const G4ThreeVector & Pini,
           G4ThreeVector & Pf1,  G4ThreeVector & Pf2,
           const AngDisGenerator & generator)
{
  if(Mini<Mf1+Mf2){
    std::cerr << "Mini < Mf1+Mf2 Mini=" << Mini/CLHEP::GeV << "GeV/c2 "
	      << "Mf1=" <<  Mf1/CLHEP::GeV << "GeV/c2 "
	      << "Mf2=" <<  Mf2/CLHEP::GeV << "GeV/c2 " << std::endl;
    return false;
  }

  G4ThreeVector beta(Pini/sqrt(Mini*Mini+Pini.mag2()));

  double Ecmf1=(Mini*Mini-Mf2*Mf2+Mf1*Mf1)/Mini*0.5;
  double Ecmf2=(Mini*Mini-Mf1*Mf1+Mf2*Mf2)/Mini*0.5;
  double Pcm=sqrt((Mini*Mini-(Mf1+Mf2)*(Mf1+Mf2))*
		  (Mini*Mini-(Mf1-Mf2)*(Mf1-Mf2)))/Mini*0.5;

  G4ThreeVector UnitDir=generator.GenerateDirection();

  UnitDir.rotateUz(Pini.unit());

  //  std::cout << "UnitDir=" << UnitDir << " Pini.unit()="
  //	    << Pini.unit() << std::endl;


  G4ThreeVector Pcmf1 =  Pcm*UnitDir;
  G4ThreeVector Pcmf2 = -Pcm*UnitDir;

  G4LorentzVector LVf1(Pcmf1, Ecmf1), LVf2(Pcmf2, Ecmf2);

  LVf1.boost(beta); Pf1=LVf1.vect();
  LVf2.boost(beta); Pf2=LVf2.vect();

  //   std::cout << "CosTCM=" << cost << " PhiCM=" << phi/degree
  //    << " degree" << std::endl;
  // std::cout << "PCM=" << Pcm/CLHEP::GeV << " " << Pcmf1/CLHEP::GeV << " -->"
  //    	    << Pf1/CLHEP::GeV << " " << Pf1.mag()/CLHEP::GeV << std::endl;
  return true;
}

//_____________________________________________________________________________
bool
Scattering2Body_theta(double Mi1, double Mi2, double Mf1, double Mf2,
                      const G4ThreeVector & Pini1,const G4ThreeVector & Pini2,
                      G4ThreeVector & Pf1,  G4ThreeVector & Pf2,
                      double & theta_CM,
                      const AngDisGenerator & generator)
{
  //  std::cout << "Pini=" << Pini/CLHEP::GeV << "GeV/c" << std::endl;
  TLorentzVector Pini1_lv ;
  TLorentzVector Pini2_lv ;
  TLorentzVector Pf1_lv  ;
  TLorentzVector Pf2_lv  ;

  TLorentzVector Total_lv;
  TLorentzVector Pf1_cm_lv;

  TVector3 Pini1_mom = TVector3(Pini1.x(),Pini1.y(),Pini1.z());
  TVector3 Pini2_mom = TVector3(Pini2.x(),Pini2.y(),Pini2.z());
  TVector3 Pf1_mom = TVector3(Pf1.x(),Pf1.y(),Pf1.z());
  TVector3 Pf2_mom = TVector3(Pf2.x(),Pf2.y(),Pf2.z());

  Pini1_lv.SetVectM(Pini1_mom,Mi1);
  Pini2_lv.SetVectM(Pini2_mom,Mi2);
  Pf1_lv.SetVectM(Pf1_mom,Mf1);
  Pf2_lv.SetVectM(Pf2_mom,Mf2);
  Total_lv = Pini1_lv + Pini2_lv;

  double masses[2]={Mf1,Mf2};
  // G4cout<<"Mi1:"<<Mi1<<G4endl;
  // G4cout<<"Mi2:"<<Mi2<<G4endl;
  // G4cout<<"Mf1:"<<Mf1<<G4endl;
  // G4cout<<"Mf2:"<<Mf2<<G4endl;


  if(Total_lv.Mag() < Mf1+Mf2){
    return false;
  }

  TGenPhaseSpace event1;
  event1.SetDecay(Total_lv,2,masses);
  //  G4cout<<"--------pbeb=="<<Pini1_lv.E() <<G4endl;

  while(1){

    // double weight1 = event1.Generate();
    Pf1_lv= *event1.GetDecay(0);
    Pf2_lv= *event1.GetDecay(1);

    TVector3 b = -Total_lv.BoostVector();

    Pf1_cm_lv=Pf1_lv;
    Pf1_cm_lv.Boost(b);
    double y=G4UniformRand();

    //G4cout<<"dist_func"<<generator.GetDfuncVal(Pf1_cm_lv.Vect().CosTheta())<<G4endl;
    //G4cout<<"ooooooooo"<<Mf2<<"oooooooooooo"<<G4endl;
    if(y<generator.GetDfuncVal(Pf1_cm_lv.Vect().CosTheta())){
      //G4cout<<" coscsos== "<<Pf1_cm_lv.Vect().CosTheta()<<G4endl;
      break;
    }
  };
  //G4cout<<"gen-----Pf1_lveb =="<<Pf1_lv.Vect().Mag() <<G4endl;
  //G4cout<<"gen-----pf1_x =="<< (Pf1_lv.Vect()).X() <<G4endl;
  //G4cout<<"gen-----pf1_y =="<< (Pf1_lv.Vect()).Y() <<G4endl;
  //G4cout<<"gen-----pf1_z =="<< (Pf1_lv.Vect()).Z() <<G4endl;
  //G4cout<<"gen-----pf1_e =="<< Pf1_lv.E() <<G4endl;
  //G4cout<<"gen-----missing =="<< (Pini1_lv + Pini2_lv + Pf1_lv *(-1) ).Mag() <<G4endl;

  theta_CM = Pf1_cm_lv.Vect().Theta();
  //std::cout<<"theta_CM _generator="<<theta_CM<<std::endl;
  //  double (*func)(double)=func_in;
  //

  Pf1 = G4ThreeVector(Pf1_lv.Vect().X(), Pf1_lv.Vect().Y() ,Pf1_lv.Vect().Z());
  Pf2 = G4ThreeVector(Pf2_lv.Vect().X(), Pf2_lv.Vect().Y() ,Pf2_lv.Vect().Z());

  return true;
}

//_____________________________________________________________________________
bool
Scattering3Body_theta(double Mi1, double Mi2,
                      double Mf1, double Mf2,double Mf3,
                      const G4ThreeVector & Pini1,const G4ThreeVector & Pini2,
                      G4ThreeVector & Pf1,
                      G4ThreeVector & Pf2, G4ThreeVector &Pf3,
                      double & theta_CM, const AngDisGenerator&)
{
  //  std::cout << "Pini=" << Pini/CLHEP::GeV << "GeV/c" << std::endl;
  TLorentzVector Pini1_lv ;
  TLorentzVector Pini2_lv ;
  TLorentzVector Pf1_lv  ;
  TLorentzVector Pf2_lv  ;
  TLorentzVector Pf3_lv  ;

  TLorentzVector Total_lv;
  TLorentzVector Pf1_cm_lv;

  TVector3 Pini1_mom = TVector3(Pini1.x(),Pini1.y(),Pini1.z());
  TVector3 Pini2_mom = TVector3(Pini2.x(),Pini2.y(),Pini2.z());
  TVector3 Pf1_mom = TVector3(Pf1.x(),Pf1.y(),Pf1.z());
  TVector3 Pf2_mom = TVector3(Pf2.x(),Pf2.y(),Pf2.z());
  TVector3 Pf3_mom = TVector3(Pf3.x(),Pf3.y(),Pf3.z());

  Pini1_lv.SetVectM(Pini1_mom,Mi1);
  Pini2_lv.SetVectM(Pini2_mom,Mi2);
  Pf1_lv.SetVectM(Pf1_mom,Mf1);
  Pf2_lv.SetVectM(Pf2_mom,Mf2);
  Pf3_lv.SetVectM(Pf3_mom,Mf3);
  Total_lv = Pini1_lv + Pini2_lv;

  double masses[3]={Mf1,Mf2,Mf3};
  // G4cout<<"Mi1:"<<Mi1<<G4endl;
  // G4cout<<"Mi2:"<<Mi2<<G4endl;
  // G4cout<<"Mf1:"<<Mf1<<G4endl;
  // G4cout<<"Mf2:"<<Mf2<<G4endl;
  // G4cout<<"Mf3:"<<Mf3<<G4endl;

  if(Total_lv.Mag() < Mf1+Mf2+Mf3){
    return false;
  }

  TGenPhaseSpace event1;
  event1.SetDecay(Total_lv,3,masses);
  //  G4cout<<"--------pbeb=="<<Pini1_lv.E() <<G4endl;

  while(1){

    double weight1 = event1.Generate();
    Pf1_lv= *event1.GetDecay(0);
    Pf2_lv= *event1.GetDecay(1);
    Pf3_lv= *event1.GetDecay(2);

    TVector3 b = -Total_lv.BoostVector();

    Pf1_cm_lv=Pf1_lv;
    Pf1_cm_lv.Boost(b);
    double y=G4UniformRand();

    // G4cout<<"weight=="<<weight1<<G4endl;
    if(y<weight1){
      break;
    }

    //G4cout<<"dist_func"<<generator.GetDfuncVal(Pf1_cm_lv.Vect().CosTheta())<<G4endl;
    //G4cout<<"ooooooooo"<<Mf2<<"oooooooooooo"<<G4endl;
    //if(y<generator.GetDfuncVal(Pf1_cm_lv.Vect().CosTheta())){
    //G4cout<<" coscsos== "<<Pf1_cm_lv.Vect().CosTheta()<<G4endl;
    //break;
    //}
  };
  //G4cout<<"gen-----Pf1_lveb =="<<Pf1_lv.Vect().Mag() <<G4endl;
  //G4cout<<"gen-----pf1_x =="<< (Pf1_lv.Vect()).X() <<G4endl;
  //G4cout<<"gen-----pf1_y =="<< (Pf1_lv.Vect()).Y() <<G4endl;
  //G4cout<<"gen-----pf1_z =="<< (Pf1_lv.Vect()).Z() <<G4endl;
  //G4cout<<"gen-----pf1_e =="<< Pf1_lv.E() <<G4endl;
  //G4cout<<"gen-----missing =="<< (Pini1_lv + Pini2_lv + Pf1_lv *(-1) ).Mag() <<G4endl;

  //  double (*func)(double)=func_in;
  //

  theta_CM = Pf1_cm_lv.Vect().Theta();

  Pf1 = G4ThreeVector(Pf1_lv.Vect().X(), Pf1_lv.Vect().Y() ,Pf1_lv.Vect().Z());
  Pf2 = G4ThreeVector(Pf2_lv.Vect().X(), Pf2_lv.Vect().Y() ,Pf2_lv.Vect().Z());
  Pf3 = G4ThreeVector(Pf3_lv.Vect().X(), Pf3_lv.Vect().Y() ,Pf3_lv.Vect().Z());

  // G4cout<<"Pf1= "<<Pf1.mag()<<G4endl;
  // G4cout<<"Pf2= "<<Pf2.mag()<<G4endl;
  // G4cout<<"Pf3= "<<Pf3.mag()<<G4endl;

  return true;
}

//_____________________________________________________________________________
bool
Decay3BodyPhaseSpace(double Mini, double Mf1, double Mf2, double Mf3,
                     const G4ThreeVector& Pini,
                     G4ThreeVector& Pf1,  G4ThreeVector& Pf2,
                     G4ThreeVector& Pf3)
{

  TLorentzVector Pini1_lv ;
  TLorentzVector Pf1_lv  ;
  TLorentzVector Pf2_lv  ;
  TLorentzVector Pf3_lv  ;

  TLorentzVector Total_lv;
  TLorentzVector Pf1_cm_lv;

  TVector3 Pini1_mom = TVector3(Pini.x(),Pini.y(),Pini.z());
  TVector3 Pf1_mom = TVector3(Pf1.x(),Pf1.y(),Pf1.z());
  TVector3 Pf2_mom = TVector3(Pf2.x(),Pf2.y(),Pf2.z());
  TVector3 Pf3_mom = TVector3(Pf3.x(),Pf3.y(),Pf3.z());

  Pini1_lv.SetVectM(Pini1_mom,Mini);
  Pf1_lv.SetVectM(Pf1_mom,Mf1);
  Pf2_lv.SetVectM(Pf2_mom,Mf2);
  Pf3_lv.SetVectM(Pf3_mom,Mf3);
  Total_lv = Pini1_lv;

  double masses[3]={Mf1,Mf2,Mf3};

  if(Total_lv.Mag() < Mf1+Mf2+Mf3){
    return false;
  }

  TGenPhaseSpace event1;
  event1.SetDecay(Total_lv,3,masses);

  while(1){

    double weight1 = event1.Generate();
    Pf1_lv= *event1.GetDecay(0);
    Pf2_lv= *event1.GetDecay(1);
    Pf3_lv= *event1.GetDecay(2);

    double y=G4UniformRand();

    if(y<weight1){
      break;
    }

  };


  Pf1 = G4ThreeVector(Pf1_lv.Vect().X(), Pf1_lv.Vect().Y() ,Pf1_lv.Vect().Z());
  Pf2 = G4ThreeVector(Pf2_lv.Vect().X(), Pf2_lv.Vect().Y() ,Pf2_lv.Vect().Z());
  Pf3 = G4ThreeVector(Pf3_lv.Vect().X(), Pf3_lv.Vect().Y() ,Pf3_lv.Vect().Z());

  return true;
}

// -*- C++ -*-

#include "ResHypTPC.hh"

#include <TH1.h>
#include <TF1.h>

#include <Randomize.hh>

#include "AnaManager.hh"

//_____________________________________________________________________________
ResHypTPC::ResHypTPC()
{
}


//_____________________________________________________________________________
ResHypTPC::ResHypTPC(G4double t_pad_size, G4double t_pad_length,
                     G4double t_threshold, G4double diff_T,
                     G4double smearing)
{
  //  //  const G4double ResHypTPC::y_resolution= 0.5;
  y_resolution= 0.5;
  sigma_amp = 13.3/53.9;
  neff = 26;//MIP
  neff_nmpv_correction = 0.7*1.2;// = 0.7*1.2;
  diffuse_GEM = 0.1;
  init();

  pad_size = t_pad_size;
  pad_length = t_pad_length;
  threshold = t_threshold;
  diffusion_T = diff_T;
  const_smearing=smearing;

  f_n_drift_electron = new TF1("f_n_drift_electron","landau(0)",0,ncutoff);
  f_n_drift_electron->SetParameter(0,1);
  f_n_drift_electron->SetParameter(1,nmpv);
  f_n_drift_electron->SetParameter(2,nsigma);

  //  f_diffusion = new TF1("f_diffusion","gaus(0)",-pad_size*10,pad_size*10);
  //  f_diffusion->SetParameter(0,1);

  f_GEM_avalanche = new TF1("f_GEM_avalanche","gaus(0)",-pad_size*10,pad_size*10);

  //    f_const_smearing = new TF1("f_const_smearing","gaus(0)",-pad_size*10,pad_size*10);
  //    f_const_smearing->SetParameter(0,1);
  //    f_const_smearing->SetParameter(1,0);
  //    f_const_smearing->SetParameter(2,const_smearing);
}

//_____________________________________________________________________________
ResHypTPC::~ResHypTPC()
{
  delete f_n_drift_electron;
  //    delete f_diffusion;
  delete f_GEM_avalanche;
  //    delete f_const_smearing;
}

//_____________________________________________________________________________
G4double
ResHypTPC::getXDeviation(G4int &n_electron, G4int &n_pad, G4double &x_rms,
                         G4double x_track, G4double y_track,
                         G4double dxdz_track=0, G4double dydz_track=0)
{
  if (debug) std::cout << "in getXDeviation" << std::endl;
  G4double path_length = pad_length*sqrt(1.+dxdz_track*dxdz_track+dydz_track*dydz_track);
  //  G4cout<<path_length<<G4endl;
  //RandLandau
  n_electron = (int)(f_n_drift_electron->GetRandom()*path_length/10.+0.5);
  //  n_electron = (int)(CLHEP::RandLandau::shoot(nmpv,nsigma)*path_length/10.+0.5);
  //  n_electron = (int)(f_n_drift_electron->CLHEP::RandFlat::shoot()*path_length/10.+0.5);
  //  G4cout<<"test random:"<<CLHEP::RandLandau::shoot()<<G4endl;
  //  G4cout<<"f_n_drift_electron:"<<f_n_drift_electron->GetRandom()<<G4endl;
  //  G4cout<<"n_electron:"<<n_electron<<G4endl;

  if (debug) std::cout << "n_electron = " << n_electron << std::endl;
  //  G4cout << "n_electron = " << n_electron << G4endl;
  //  G4cout << "y_track = " << y_track << G4endl;
  G4double dz = pad_length/(double)n_electron;

  TH1F *hpad = new TH1F("hpad","hpad",20,-pad_size*10,pad_size*10);
  for (int ie=0;ie<n_electron;ie++) {
    G4double z = dz*((double)ie+0.5);
    //    G4double x_init = x_track + dxdz_track*z;
    //    G4double y_init = y_track + dydz_track*z;
    //    x,y is at the center of the pad length
    G4double x_init = x_track + dxdz_track*(z-pad_length/2.);
    G4double y_init = y_track + dydz_track*(z-pad_length/2.)
      ;
    //    f_diffusion->SetParameter(1,x_init);
    //    f_diffusion->SetParameter(2,diffusion_T*sqrt(y_init/10.));
    G4double x_diffused = CLHEP::RandGauss::shoot(x_init,diffusion_T*sqrt(y_init/10.));
    //    G4double x_diffused = f_diffusion->GetRandom();
    f_GEM_avalanche->SetParameter(0,1./sqrt(2.*3.141592654)/diffuse_GEM);
    f_GEM_avalanche->SetParameter(1,x_diffused);
    f_GEM_avalanche->SetParameter(2,diffuse_GEM);
    for (int ipad=-10;ipad<10;ipad++) {
      G4double x0 = (double)ipad*pad_size;
      G4double x1 = ((double)ipad+1.)*pad_size;
      hpad->Fill((x0+x1)/2.,f_GEM_avalanche->Integral(x0,x1));
      //      G4cout<<x0<<":"<<x1<<G4endl;
    }
  }


  if (debug) std::cout << "end filling hpad" << std::endl;

  G4double npeak = hpad->GetMaximum();
  G4double nth = npeak*threshold;
  n_sum = 0;
  x_mean = 0;
  x2_mean = 0;
  n_pad = 0;
  for (int ipad=0;ipad<20;ipad++) {
    G4double n_in_pad = hpad->GetBinContent(ipad+1);
    G4double x_pad = hpad->GetBinCenter(ipad+1);
    if (n_in_pad>=nth) {
      n_sum   += n_in_pad;
      x_mean  += (n_in_pad*x_pad);
      x2_mean += (n_in_pad*x_pad*x_pad);
      n_pad ++;
    }
  }

  x_rms = 0;

  if (n_sum>0) {
    x_mean /= n_sum;
    x2_mean /= n_sum;
    if ((x2_mean - x_mean*x_mean)<0) {
      //      std::cerr << "x2_mean - x_mean*x_mean = " << x2_mean - x_mean*x_mean << std::endl;
      //      std::cerr << "x2_mean*n_sum - x_mean*n_sum = " << x2_mean*n_sum << " " << x_mean*n_sum << std::endl;
      x_rms = sqrt(x2_mean - x_mean*x_mean);
    } else {
      x_rms = sqrt(x2_mean - x_mean*x_mean);
    }
  }

  delete hpad;
  //  G4cout<<const_smearing<<G4endl;
  //  G4cout<<CLHEP::RandGauss::shoot(0.,const_smearing)<<G4endl;
  if (const_smearing>0) {
    //        x_mean += f_const_smearing->GetRandom();
    //    G4cout<<const_smearing<<G4endl;
    //    G4cout<<CLHEP::RandGauss::shoot(0.,const_smearing)<<G4endl;
    x_mean += CLHEP::RandGauss::shoot(0.,const_smearing);
    //    G4cout<<x_mean<<G4endl;
  }
  return x_mean-x_track;
}

//_____________________________________________________________________________
G4double
ResHypTPC::getYDeviation(G4double /* y_track */)
{
  return y_resolution;
}

//_____________________________________________________________________________
G4double
GetTransverseRes(G4double y_pos)
{
  //  double s0 = 0.204;// mm HIMAC result
  double s0 = 0.199;// mm HIMAC result
  double Dt = 0.18;//mm/sqrt(cm) at 1T
  double L_D = 30.+(y_pos*0.1);//cm
  // double N_eff = 42.8;
  double N_eff = 42.1;
  //double A = 0.0582*0.01;//m-1 -> cm-1
  double A = 0.055*0.01;//m-1 -> cm-1
  double e_ALD = exp(-1.*A*L_D);
  //  double sT2 = s0*s0 + (Dt*Dt*L_D/(N_eff*e_ALD));
  double sT2 = (Dt*Dt*L_D/(N_eff*e_ALD));
  double sT = sqrt(sT2);
  return sT;
}

G4double
GetVerticalRes( G4double y, G4double* par ){

	double p0 = par[0],p1 = par[1],p2 = par[2],p3=par[3];
	double dl = y + 300;
	double val = sqrt( 
			p0*p0 + p2*p2*dl/(p3*exp(-p1*dl))
			);
	return val;
}
G4double
GetHorizontalRes(G4double t,G4double y, G4double* par){
	double p0 = par[0],p1 = par[1],p2 = par[2],p3=par[3],p4=par[4],p5=par[5];
	double dl = y + 300;
	double val = 1.2 * sqrt(
			p0*p0 + p2*p2*dl/(p3*exp(-p1*dl))
			+p4*p4*tan(t)*tan(t)/(12*p5));
	return val;
}
G4ThreeVector
GetResVector(G4ThreeVector pos,G4ThreeVector mom, G4double* par_y, G4double* par_t){
	G4ThreeVector XZ(1,0,1);
	double PadAngle = atan2(-pos.x(),(pos.z()));
	double TrackAngle = PadAngle - atan2(-mom.x(),mom.z());
	double res_y = GetVerticalRes(pos.y(), par_y);
	double res_t = GetHorizontalRes(TrackAngle,pos.y(),par_t);
	double res_x =abs( cos(PadAngle)) *res_t; 	
	double res_z =abs( sin(PadAngle)) *res_t; 	
	return G4ThreeVector(res_x,res_y,res_z);
}
G4ThreeVector
GetSmearingVector(G4ThreeVector pos,G4ThreeVector mom, G4double* par_y, G4double* par_t){
	auto res_vect = GetResVector(pos,mom,par_y,par_t);
	double res_y = res_vect.y();
	double res_t = hypot(res_vect.x(),res_vect.z());
	double PadAngle = atan2(pos.x(),(pos.z()));
	
	double smear_y = CLHEP::RandGauss::shoot(0,res_y);
	double smear_t = CLHEP::RandGauss::shoot(0,res_t);
	double smear_x = cos(PadAngle) *smear_t; 	
	double smear_z = -sin(PadAngle) *smear_t; 	
	return G4ThreeVector(smear_x,smear_y,smear_z);
}

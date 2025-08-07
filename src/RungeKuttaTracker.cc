// -*- C++ -*-

/**
 *  2013.6.4. S.Hwang
 *  2013.6.17. revised the program as same as LEPS ana
 */

#include "RungeKuttaTracker.hh"

#include <globals.hh>

#include <TMath.h>

//#include "filter.h"
#include "MathTools.hh"
#include "track.hh"

//#include "nrutil.h"

#define NR_END 1
#define FREE_ARG char*
#define SWAP(a,b){temp=(a);(a)=(b);(b)=temp;}

//_____________________________________________________________________________
void
nrerror(G4String error_text)
{
  std::cout << "error in where?:" << error_text << std::endl;
  std::exit(1);
}

//_____________________________________________________________________________
int*
ivector(long nl, long nh)
{
  int *v;
  v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
  if(!v) nrerror("");
  return v-nl+NR_END;
}

//_____________________________________________________________________________
void
free_ivector(int *v, long nl, long /* nh */)
{
  free((FREE_ARG) (v+nl-NR_END));
}

//_____________________________________________________________________________
int
gaussj(double a[6][6], int n, double* /* b */, int /* m */)
{
  int *indxc,*indxr,*ipiv;
  int i,icol=0,irow=0,j,k,l,ll;
  double big,dum,pivinv;
  double temp;
  indxc=ivector(1,n);
  indxr=ivector(1,n);
  ipiv=ivector(1,n);
  for(j=1;j<=n;j++) ipiv[j]=0;
  for(i=1;i<=n;i++){
    big=0.0;
    for(j=1;j<=n;j++)
      if(ipiv[j]!=1)
	for(k=1;k<=n;k++){
	  if(ipiv[k]==0){
	    if(fabs(a[j][k])>=big){
	      big=fabs(a[j][k]);
	      irow=j;
	      icol=k;
	    }
	  } else if (ipiv[k]>1) nrerror("error");
	}
    ++(ipiv[icol]);
    if(irow!=icol){
      for(l=1;l<=n;l++) SWAP(a[irow][l],a[icol][l]);
      /*	  for(l=1;l<=m;l++) SWAP(b[irow][l],b[icol][l]);*/
    }
    indxr[i]=irow;
    indxc[i]=icol;
    if(a[icol][icol] ==0.0) nrerror("error");
    pivinv=1.0/a[icol][icol];
    a[icol][icol]=1.0;
    for(l=1;l<=n;l++){
      a[icol][l]*=pivinv;
    }
    /*	for(l=1;l<=m;l++) b[icol][l]*=pivinv;*/
    for(ll=1;ll<=n;ll++)
      if(ll != icol){
	dum=a[ll][icol];
	a[ll][icol]=0.0;
	for(l=1;l<=n;l++) a[ll][l] -= a[icol][l]*dum;
	/*		for(l=1;l<=m;l++) b[ll][l] -= b[icol][l]*dum;*/
      }
  }
  for(l=n;l>=1;l--){
    if(indxr[l]!=indxc[l])
      for(k=1;k<=n;k++)
	SWAP(a[k][indxr[l]],a[k][indxc[l]]);
  }

  free_ivector(ipiv,1,n);
  free_ivector(indxr,1,n);
  free_ivector(indxc,1,n);

  return 0;


}

//_____________________________________________________________________________
RungeKuttaTracker::RungeKuttaTracker(int c_use, Track* aTrack)
{
  RungeKuttaTracking(c_use, aTrack);
}

//_____________________________________________________________________________
RungeKuttaTracker::~RungeKuttaTracker()
{
}

//_____________________________________________________________________________
void
RungeKuttaTracker::RungeKuttaTracking(int /* c_use */, Track* aTrack)
{
  ///////////leps ana method
  // int ndf=2*(aTrack->numHits-aTrack->nout)-5;
  std::cout<<"start tracking-"<<std::endl;
  std::cout<<"numhits:"<<aTrack->numHits<<std::endl;
  /*  double dc_pos[22]=
      {0.,0.,0.,0.,
      0.,0.,0.,0.,
      850.,870.,890.,910.,930.,950.,
      2408.5,2416.3,2433.7,2441.5,
      2628.5,2660.5,2691.5,2721.5
      };
  */
  /*  double dc_res[22]=
      {0.050,0.05,0.05,0.05,//ssd1
      0.05,0.05,0.05,0.05,//ssd2
      0.200,0.200,0.200,0.200,0.200,0.200,//dc1
      0.200,0.200,0.200,0.200,//dc2
      0.200,0.200,0.200,0.200//dc3
      };
  */
  int iteration=0;
  double rk_par0[20];
  double rk_par1[20]={0};
  // double rk_hit[20][24];//[para name][number of layer]
  // double rk_hit1[20][24];//[para name][number of layer]
  // double chi2=0.;
  // double chi2_old=0.;
  // double chi2_new=0.;
  // int ierr=0.;

  //  1st Fit
  //  let's define rk_par0 for LEPSana
  //1: p
  //2:charge
  //3:z
  //4:x
  //5:y
  //6:dxdz
  //7:dydz
  //8:chi2
  //9:
  //10:nhit
  //
  for(int i=0;i<5;i++){
    rk_par0[i]=aTrack->rKInitPara[i];
    std::cout<<i<<":rk init:"<<rk_par0[i]<<std::endl;
  }

  iteration=1;
  ///include RKstep
  // ierr=RungeKuttaFit(c_use, iteration,aTrack,rk_par0,rk_par1,rk_hit);
  //  std::cout<<"number of step:"<<ierr<<std::endl;
  //prepare next step
  //  std::cout<<"after 1st rungekutta fit"<<std::endl;
  for(int i=0;i<5;i++){
    //    std::cout<<rk_par0[i]<<", "<<rk_par1[i]<<std::endl;
    rk_par0[i]=rk_par1[i];
  }


  //from 2nd fit to max itr(100)
  double dchi2[1000]={0.};
  double mom[1000]={0.};
  double prbchi2[1000]={0.};
  double calchi2[1000]={0.};
  double chi2old=-1;
  double chi2new=-1;
  for(int i=0;i<50;i++){
    iteration=iteration+1;
    chi2old=rk_par1[5]/(rk_par1[6]-5);///previous chi2
    // ierr=RungeKuttaFit(c_use, iteration,aTrack,rk_par0,rk_par1,rk_hit);
    chi2new=rk_par1[5]/(rk_par1[6]-5.);
    mom[iteration]=rk_par1[4];
    std::cout<<"km chi2/ndf:"<<rk_par1[5]/(rk_par1[6]-5.)<<std::endl;
    dchi2[iteration]=(chi2new-chi2old)/chi2old;
    prbchi2[iteration]=TMath::Prob(rk_par1[5],int(rk_par1[6]-5.));
    calchi2[iteration]=(rk_par1[5]/int(rk_par1[6]-5.));
    std::cout<<"chi2old:"<<chi2old<<std::endl;
    std::cout<<"chi2new:"<<chi2new<<std::endl;
    //    std::cout<<"dchi2:"<<dchi2[iteration]<<std::endl;
    std::cout<<"number of hits:"<<rk_par1[6]<<std::endl;
    std::cout<<"ndf:"<<int(rk_par1[6]-5.)<<std::endl;
  }
  std::cout<<"display dchi2"<<std::endl;
  std::cout<<"iteration:dchi2:prbchi2:chi2:mom"<<std::endl;

  for(int i=2;i < iteration;i++){
    std::cout<<i<<" : "<<dchi2[i]<<" : "<< prbchi2[i] <<" : "<< calchi2[i]<<" : " << 1./mom[i]<<std::endl;
  }
}

//_____________________________________________________________________________
void
RungeKuttaTracker::RungeKuttaFieldIntegral(int /* c_use */,
                                           double* /* init_par */,
                                           double* /* final_par */,
                                           G4ThreeVector& /* B */)
{
}

//_____________________________________________________________________________
void
RungeKuttaTracker::RungeKuttaFieldIntegral(int /* c_use */,
                                           double* /* init_par */,
                                           double* /* final_par */,
                                           G4ThreeVector& /* B */,
                                           G4ThreeVector& /* dBdY */,
                                           G4ThreeVector& /* dBdZ */)
{
}

//_____________________________________________________________________________
int
RungeKuttaTracker::RungeKuttaFit(int c_use, int iteration, Track* aTrack,
                                 double* rk_par0, double* rk_par1,
                                 double* /* rk_hit */)
{
  // double chi2=0.;
  //  double u0[2];
  //  double dudz0[2];
  //  double qp0;

  //  double u[2];
  //  double dudz[2];
  //  double qp;

  //  double dudw0[2][5];
  //  double ddudw0[2][5];
  //  double dudw[2][5];
  //  double ddudw[2][5];
  // lets
  if(iteration==1.){
    KMFinit1();
  }
  u0[0]=rk_par0[0];
  u0[1]=rk_par0[1];
  dudz0[0]=rk_par0[2];
  dudz0[1]=rk_par0[3];
  qp0=rk_par0[4];

  kmf_.KMFpar0[0][0]=u0[0];
  kmf_.KMFpar0[1][0]=u0[1];
  kmf_.KMFpar0[2][0]=dudz0[0];
  kmf_.KMFpar0[3][0]=dudz0[1];
  kmf_.KMFpar0[4][0]=qp0;

  ///
  for(int i=0;i<5;i++){
    for(int j=0;j<5;j++){
      kmf_.Fplane[i][j]=0.;
    }
  }

  z0=aTrack->x[0][2]-5.; //5 mm offset
  //  z=0.;
  /*
    double dc_pos[22]=
    {0.,0.,0.,0.,
    0.,0.,0.,0.,
    850.,870.,890.,910.,930.,950.,
    2408.5,2416.3,2433.7,2441.5,
    2628.5,2660.5,2691.5,2721.5
    };
  */
  /*
    double dc_res[22]=
    {0.050,0.05,0.05,0.05,//ssd1
    0.05,0.05,0.05,0.05,//ssd2
    0.200,0.200,0.200,0.200,0.200,0.200,//dc1
    0.200,0.200,0.200,0.200,//dc2
    0.200,0.200,0.200,0.200//dc3
    };
  */

  double h=1.;//first try
  int num_step=0;

  //start tracking
  num_step=0;
  int det_plane_hit=0;
  double path_length=0;

  for(int ihit=0;ihit<aTrack->numHits;ihit++){
    //init dudw0 and ddudw0
    for(int k=0;k<2;k++){
      for(int j=0;j<5;j++){
	dudw0[k][j]=0.;
	ddudw0[k][j]=0.;
      }
    }
    dudw0[0][0]=1.;
    dudw0[1][1]=1.;
    ddudw0[0][2]=1.;
    ddudw0[1][3]=1.;
    double z_hit=aTrack->x[ihit][2];


    while(1){
      num_step +=1;
      if((z0+h)<z_hit){
	h=1.;
	det_plane_hit=0;
      }else if((z0+h)>=z_hit){
	h=z_hit-z0;
	det_plane_hit=1;
      }else{
	std::cout<<"something is wrong"<<std::endl;
	det_plane_hit=0;
	exit(1);
      }

      //////////runge kutta step///////////////
      std::cout<<"number of step:"<<num_step<<std::endl;
      std::cout<<"run rungekuttastep:"<<num_step<<std::endl;
      RungeKuttaStep(c_use, qp0,h,
                     z0,u0,dudz0,dudw0,ddudw0);
      //      prepare next step
      z=z0+h;
      z0=z;
      for(int j=0;j<2;j++){
	u0[j]=u[j];
	dudz0[j]=dudz[j];
	for(int k=0;k<5;k++){
	  dudw0[j][k]=dudw[j][k];
	  ddudw0[j][k]=ddudw[j][k];
	}
      }

      path_length=path_length+sqrt(pow(u[0]-u0[0],2) + pow(u[1]-u0[1],2) + h*h);
      double x_hit=aTrack->x[ihit][0];
      double y_hit=aTrack->x[ihit][1];
      double res=aTrack->res[ihit];
      int idwi0=aTrack->idwi[ihit];
      /*
	std::cout<<"0000000000000000000000000"<<std::endl;
	std::cout<<"0000000000000000000000000"<<std::endl;
	std::cout<<"0000000000000000000000000"<<std::endl;
	std::cout<<"0000000000000000000000000"<<std::endl;
	std::cout<<"ihit-->"<<ihit<<std::endl;
	std::cout<<"idwi-->"<<idwi0<<std::endl;
	std::cout<<"0000000000000000000000000"<<std::endl;
	std::cout<<"0000000000000000000000000"<<std::endl;
	std::cout<<"0000000000000000000000000"<<std::endl;
	std::cout<<"0000000000000000000000000"<<std::endl;
      */
      //      if(z>=z_hit){
      //	chi2=chi2+(pow(x_hit-u0[0],2)+pow(y_hit-u0[1],2))/(res*res);
      //	std::cout<<"-----------------------------"<<std::endl;
      //	std::cout<<"res:"<<res<<std::endl;
      //	std::cout<<"chi2:"<<chi2<<std::endl;
      //	std::cout<<"-----------------------------"<<std::endl;
      //      }
      ///detector plane
      if(det_plane_hit==1.){
	for(int j=0;j<5;j++){
	  kmf_.Fplane[0][j] = dudw[0][j];
	  kmf_.Fplane[1][j] = dudw[1][j];
	  kmf_.Fplane[2][j] = ddudw[0][j];
	  kmf_.Fplane[3][j] = ddudw[1][j];
	  kmf_.Fplane[4][j] = 0.;
	}
	kmf_.Fplane[4][4]=1.0;
	//////////
	std::cout<<ihit<<"th Fplane"<<std::endl;
	for(int j=0;j<5;j++){
#if 0
	  std::cout<<kmf_.Fplane[j][0]<<" ,";
	  std::cout<<kmf_.Fplane[j][1]<<" ,";
	  std::cout<<kmf_.Fplane[j][2]<<" ,";
	  std::cout<<kmf_.Fplane[j][3]<<" ,";
	  std::cout<<kmf_.Fplane[j][4]<<std::endl;
#endif
	}

	//	std::cout<<"----------------------------"<<std::endl;

	if(ihit==0){
	  kmf_.KMFpar1[0][ihit]=u[0];
	  kmf_.KMFpar1[1][ihit]=u[1];
	  kmf_.KMFpar1[2][ihit]=dudz[0];
	  kmf_.KMFpar1[3][ihit]=dudz[1];
	  kmf_.KMFpar1[4][ihit]=qp0;
	  kmf_.reso[ihit]=res;
	  kmf_.idwi=idwi0;
	  if(idwi0==1){
	    kmf_.meas[ihit]=x_hit;//
	  }else if(idwi0==2){
	    kmf_.meas[ihit]=y_hit;//
	  }
	  kmf_.measz[ihit]=z_hit;
	  //	idwi=
	  RKHits[2][ihit]=z;
	  RKHits[3][ihit]=path_length;

	  ////////////////KMF init2
	  KMFinit2();

	}else{
	  kmf_.KMFpar0[0][ihit]=u[0];
	  kmf_.KMFpar0[1][ihit]=u[1];
	  kmf_.KMFpar0[2][ihit]=dudz[0];
	  kmf_.KMFpar0[3][ihit]=dudz[1];
	  kmf_.KMFpar0[4][ihit]=qp0;

	  pGeV=fabs(1/qp0);
	  kmf_.mulsth[ihit]=0.;
	  if(idwi0==1){
	    kmf_.meas[ihit]=x_hit;//
	  }else if(idwi0==2){
	    kmf_.meas[ihit]=y_hit;//
	  }
	  kmf_.reso[ihit]=res;
	  kmf_.idwi=idwi0;
	  kmf_.path[ihit]=1.;

	  kmf_.iplane=ihit;
	  //	  std::cout<<"0000000000000000000000000"<<std::endl;
	  //	  std::cout<<"kmfilter"<<std::endl;
	  //	  std::cout<<"0000000000000000000000000"<<std::endl;
	  KMFilter();

	  u0[0]=kmf_.KMFpar1[0][ihit];
	  u0[1]=kmf_.KMFpar1[1][ihit];
	  dudz0[0]=kmf_.KMFpar1[2][ihit];
	  dudz0[1]=kmf_.KMFpar1[3][ihit];
	  qp0=kmf_.KMFpar1[4][ihit];

	  /*
	    u0[0]=u[0];
	    u0[1]=u[1];
	    dudz0[0]=dudz[0];
	    dudz0[1]=dudz[1];
	  */
	  //	  qp0=qp;



	  std::cout<<"<<<<<<< check RK --> RF & Matrix >>>>>>>"<<std::endl;
	  std::cout<<u[0]<<" : "<< kmf_.KMFpar1[0][ihit]<<" : "<<aTrack->x[ihit][0]<<std::endl;
	  std::cout<<u[1]<<" : "<< kmf_.KMFpar1[1][ihit]<<" : "<<aTrack->x[ihit][1]<<std::endl;
	  std::cout<<dudz[0]<<" : "<< kmf_.KMFpar1[2][ihit]<<" : "<<aTrack->uv[ihit][0]<<std::endl;
	  std::cout<<dudz[1]<<" : "<< kmf_.KMFpar1[3][ihit]<<" : "<<aTrack->uv[ihit][1]<<std::endl;
	  std::cout<<1./qp0<<" : "<< 1./kmf_.KMFpar1[4][ihit]<<" : "<<1./aTrack->rKInitPara[4]<<std::endl;
	  std::cout<<"----------------------------------"<<std::endl;




	  //	  std::cout<<"0000000000000000000000000"<<std::endl;
	  //	  std::cout<<"end kmfilter"<<std::endl;
	  //	  std::cout<<"0000000000000000000000000"<<std::endl;
	}

      }



      //      std::cout<<"==========================:"<<std::endl;
      //      std::cout<<"z0:"<<z<<std::endl;
      //      std::cout<<"h:"<<h<<std::endl;

      //      std::cout<<"aTrack->x[ihit][2]:"<<aTrack->x[ihit][2]<<std::endl;
      //      std::cout<<"==========================:"<<std::endl;

      if(z>=z_hit){

	//	std::cout<<"==========================:"<<std::endl;
	//	std::cout<<"z0:"<<z<<std::endl;
	//	std::cout<<"h:"<<h<<std::endl;
	//	std::cout<<"aTrack->x[ihit][2]:"<<aTrack->x[ihit][2]<<std::endl;
	//	std::cout<<"X cal hit:mea hit -->"<< u[0] <<" : "<<aTrack->x[ihit][0]<<std::endl;
	//	std::cout<<"Y cal hit:mea hit -->"<< u[1] <<" : "<<aTrack->x[ihit][1]<<std::endl;
	//	std::cout<<"step -->"<< num_step <<std::endl;
	//	std::cout<<"==========================:"<<std::endl;

	break;
      }

      if(z>2721.5 || num_step>5000){
	//	std::cout<<"exceed num of step in RK integration:"<<ihit<<" : "<<aTrack->x[ihit][2]<<std::endl;
	break;
      }
      //      std::cout<<"end of while"<<std::endl;
    }
    //    std::cout<<"end of ihit"<<std::endl;
  }

  //  smooth
  kmf_.iplane=aTrack->numHits;
  fltosm();

  for(int ihit=0;ihit<aTrack->numHits;ihit++){
    kmf_.iplane=aTrack->numHits-ihit-1;
    kmf_.idwi=aTrack->idwi[kmf_.iplane];
    //    std::cout<<"before KMsmooth"<<std::endl;
    //    std::cout<<"iplane:"<<kmf_.iplane<<std::endl;
    //    std::cout<<"idwi:"<<kmf_.idwi<<std::endl;
    KMsmooth();
  }
  /*
    u0[0]=rk_par0[0];
    u0[1]=rk_par0[1];
    dudz0[0]=rk_par0[2];
    dudz0[1]=rk_par0[3];
    qp0=rk_par0[4];
  */

  std::cout<<"--------------------------------"<<std::endl;
  std::cout<<"initial & final"<<std::endl;
  std::cout<<rk_par0[0]<<", "<<kmf_.SMTpar[0][0]<<std::endl;
  std::cout<<rk_par0[1]<<", "<<kmf_.SMTpar[1][0]<<std::endl;
  std::cout<<rk_par0[2]<<", "<<kmf_.SMTpar[2][0]<<std::endl;
  std::cout<<rk_par0[3]<<", "<<kmf_.SMTpar[3][0]<<std::endl;
  std::cout<<rk_par0[4]<<", "<<kmf_.SMTpar[4][0]<<std::endl;
  std::cout<<"--------------------------------"<<std::endl;
  /*
    rk_par0[0]=kmf_.SMTpar[0][0];
    rk_par0[1]=kmf_.SMTpar[1][0];
    rk_par0[2]=kmf_.SMTpar[2][0];
    rk_par0[3]=kmf_.SMTpar[3][0];
    rk_par0[4]=kmf_.SMTpar[4][0];
  */

  rk_par1[0]=kmf_.SMTpar[0][0];
  rk_par1[1]=kmf_.SMTpar[1][0];
  rk_par1[2]=kmf_.SMTpar[2][0];
  rk_par1[3]=kmf_.SMTpar[3][0];
  rk_par1[4]=kmf_.SMTpar[4][0];


  std::cout<<"kmchi2:"<<kmf_.kmfchi2<<std::endl;
  rk_par1[5]=kmf_.kmfchi2;
  rk_par1[6]=aTrack->numHits;
  return num_step;
}

//_____________________________________________________________________________
void
RungeKuttaTracker::RungeKuttaStep(int c_use, double /* qp0 */, double h,
                                  double /* z0 */, double* /* u0 */,
                                  double* /* dudz0 */,
                                  double dudw0_[2][5], double ddudw0_[2][5])
{
  //  std::cout<<"rungekutta step"<<std::endl;

  double u1[3];
  // double b1[3];
  double dudz1[2];

  double u2[3];
  // double b2[3];
  double dudz2[2];

  double u3[3];
  // double b3[3];
  double dudz3[2];

  double u4[3];
  // double b4[3];
  double dudz4[2];

  double A1[2][2], A2[2][2], A3[2][2], A4[2][2];
  // double C1[2][2], C2[2][2], C3[2][2], C4[2][2];
  double K1[2],K2[2],K3[2],K4[2];
  double f1[2],f2[2],f3[2],f4[2];

  ///magnet position 1105+400
  double b0[3]={0};
  double dbdu0[3][2]={{0}};
  double h2=h*h;
  z = z0+h;

  ///////////////////////////////
  ///////// 1st point /////////
  ///////////////////////////////
  for(int i=0;i<2;i++){
    u1[i]=u0[i];
    dudz1[i]=dudz0[i];
  }
  u1[2]=z0;

#if 0
  std::cout<<"+++++++++++++++in the RKStep+++++++++++++++++"<<std::endl;
  std::cout<<"z0:"<<z0<<std::endl;
  std::cout<<"u0[0]:"<<u0[0]<<std::endl;
  std::cout<<"u0[1]:"<<u0[1]<<std::endl;
  std::cout<<"dudz0[0]:"<<dudz0[0]<<std::endl;
  std::cout<<"dudz0[1]:"<<dudz0[1]<<std::endl;
  std::cout<<"qp0:"<<qp0<<std::endl;
#endif
  ///////////////////////////////
  ///////// 2nd point /////////
  ///////////////////////////////
  if(u1[2]>1550.-400. && u1[2]<1550.+400.){
    b0[0]=b0[2]=0.;
    b0[1]=-0.875*0.299792458;
    //    b0[1]=-0.875;
  }else{
    b0[0]=b0[1]=b0[2]=0.;
  }

  fACK(c_use,dudz1,b0,dbdu0,qp0);
  for(int i=0;i<2;i++){
    f1[i]=f[i];
    K1[i]=K[i];
    for(int j=0;j<2;j++){
      A1[i][j]=A[i][j];
      C1[i][j]=C[i][j];
    }
  }
#if 0
  std::cout<<"+++++++++++++++in the RKStep+++++++++++++++++"<<std::endl;
  std::cout<<"f1[0]:"<<f1[0]<<std::endl;
  std::cout<<"f1[1]:"<<f1[1]<<std::endl;
  std::cout<<"K1[0]:"<<K1[0]<<std::endl;
  std::cout<<"K1[1]:"<<K1[1]<<std::endl;
  std::cout<<"A1:"<<std::endl;
  std::cout<<"          "<<A1[0][0]<<", "<<A1[0][1]<<std::endl;
  std::cout<<"          "<<A1[1][0]<<", "<<A1[1][1]<<std::endl;
  std::cout<<"C1:"<<std::endl;
  std::cout<<"          "<<C1[0][0]<<", "<<C1[0][1]<<std::endl;
  std::cout<<"          "<<C1[1][0]<<", "<<C1[1][1]<<std::endl;
#endif
  for(int i=0;i<2;i++){
    u2[i]=u0[i]+h*dudz0[i]/2.+h2*K1[i]/8.;
    dudz2[i]=dudz0[i]+h*K1[i]/2.;
  }
  u2[2]=z0+h/2.;

  ///////////////////////////////
  ///////// 3rd point /////////
  ///////////////////////////////
  if(u2[2]>1550.-400. && u2[2]<1550.+400.){
    b0[0]=b0[2]=0.;
    b0[1]=-0.875*0.299792458;
    //    b0[1]=-0.875;
  }else{
    b0[0]=b0[1]=b0[2]=0.;
  }

  fACK(c_use,dudz2,b0,dbdu0,qp0);
  for(int i=0;i<2;i++){
    f2[i]=f[i];
    K2[i]=K[i];
    for(int j=0;j<2;j++){
      A2[i][j]=A[i][j];
      C2[i][j]=C[i][j];
    }
  }
#if 0
  std::cout<<"+++++++++++++++in the RKStep+++++++++++++++++"<<std::endl;
  std::cout<<"f2[0]:"<<f2[0]<<std::endl;
  std::cout<<"f2[1]:"<<f2[1]<<std::endl;
  std::cout<<"K2[0]:"<<K2[0]<<std::endl;
  std::cout<<"K2[1]:"<<K2[1]<<std::endl;
  std::cout<<"A2:"<<std::endl;
  std::cout<<"          "<<A2[0][0]<<", "<<A2[0][1]<<std::endl;
  std::cout<<"          "<<A2[1][0]<<", "<<A2[1][1]<<std::endl;
  std::cout<<"C2:"<<std::endl;
  std::cout<<"          "<<C2[0][0]<<", "<<C2[0][1]<<std::endl;
  std::cout<<"          "<<C2[1][0]<<", "<<C2[1][1]<<std::endl;
#endif

  for(int i=0;i<2;i++){
    u3[i]=u2[i];
    dudz3[i]=dudz0[i]+h*K2[i]/2.;
  }
  u3[2]=z0+h/2.;

  ///////////////////////////////
  ///////// 4th point /////////
  ///////////////////////////////
  if(u3[2]>1550.-400. && u3[2]<1550.+400.){
    b0[0]=b0[2]=0.;
    b0[1]=-0.875*0.299792458;
    //    b0[1]=-0.875;
  }else{
    b0[0]=b0[1]=b0[2]=0.;
  }

  fACK(c_use,dudz3,b0,dbdu0,qp0);
  for(int i=0;i<2;i++){
    f3[i]=f[i];
    K3[i]=K[i];
    for(int j=0;j<2;j++){
      A3[i][j]=A[i][j];
      C3[i][j]=C[i][j];
    }
  }
#if 0
  std::cout<<"+++++++++++++++in the RKStep+++++++++++++++++"<<std::endl;
  std::cout<<"f3[0]:"<<f3[0]<<std::endl;
  std::cout<<"f3[1]:"<<f3[1]<<std::endl;
  std::cout<<"K3[0]:"<<K3[0]<<std::endl;
  std::cout<<"K3[1]:"<<K3[1]<<std::endl;
  std::cout<<"A3:"<<std::endl;
  std::cout<<"          "<<A3[0][0]<<", "<<A3[0][1]<<std::endl;
  std::cout<<"          "<<A3[1][0]<<", "<<A3[1][1]<<std::endl;
  std::cout<<"C3:"<<std::endl;
  std::cout<<"          "<<C3[0][0]<<", "<<C3[0][1]<<std::endl;
  std::cout<<"          "<<C3[1][0]<<", "<<C3[1][1]<<std::endl;
#endif
  for(int i=0;i<2;i++){
    u4[i]=u0[i]+h*dudz0[i]+h2*K3[i]/2;
    dudz4[i]=dudz0[i]+h*K3[i];
  }
  u4[2]=z0+h;

  if(u4[2]>1550.-400. && u4[2]<1550.+400.){
    b0[0]=b0[2]=0.;
    b0[1]=-0.875*0.299792458;
    //    b0[1]=-0.875;
  }else{
    b0[0]=b0[1]=b0[2]=0.;
  }
  fACK(c_use,dudz4,b0,dbdu0,qp0);
  for(int i=0;i<2;i++){
    f4[i]=f[i];
    K4[i]=K[i];
    for(int j=0;j<2;j++){
      A4[i][j]=A[i][j];
      C4[i][j]=C[i][j];
    }
  }
#if 0
  std::cout<<"+++++++++++++++in the RKStep+++++++++++++++++"<<std::endl;
  std::cout<<"f4[0]:"<<f4[0]<<std::endl;
  std::cout<<"f4[1]:"<<f4[1]<<std::endl;
  std::cout<<"K4[0]:"<<K4[0]<<std::endl;
  std::cout<<"K4[1]:"<<K4[1]<<std::endl;
  std::cout<<"A4:"<<std::endl;
  std::cout<<"          "<<A4[0][0]<<", "<<A4[0][1]<<std::endl;
  std::cout<<"          "<<A4[1][0]<<", "<<A4[1][1]<<std::endl;
  std::cout<<"C4:"<<std::endl;
  std::cout<<"          "<<C4[0][0]<<", "<<C4[0][1]<<std::endl;
  std::cout<<"          "<<C4[1][0]<<", "<<C4[1][1]<<std::endl;
#endif
  /////calculation of u, dudz
  for(int i=0;i<2;i++){
    u[i]=u0[i]+h*dudz0[i]+h2*(K1[i] + K2[i] + K3[i])/6.;
    dudz[i]=dudz0[i]+h*(K1[i] + 2.*K2[i] + 2.*K3[i] + K4[i])/6.;
  }

  ///calculation of dudw, ddudw
  double dudv0[2]={0};
  double ddudv0[2]={0};
  double dudp0[2]={0};
  double ddudp0[2]={0};

  for(int i=0;i<5;i++){
    if(i<4){

      for(int j=0;j<2;j++){
	dudv0[j]=dudw0_[j][i];
	ddudv0[j]=ddudw0_[j][i];
	dKdw0(c_use,h,dudv0,ddudv0,A1,A2,A3,A4);///out put dK1,dK2,dK3,dK4
      }
    }else{///qp0
      for(int j=0;j<2;j++){
	dudp0[j]=dudw0_[j][i];
	ddudp0[j]=ddudw0_[j][i];
	dKdp0(c_use, h,dudp0,ddudp0,f1,f2,f3,f4,A1,A2,A3,A4);///out put dK1,dK2,dK3,dK4
      }
    }

    for(int j=0;j<2;j++){
      dudw[j][i]=dudw0_[j][i]+h*ddudw0_[j][i]+h2*(dK1[j]+dK2[j]+dK3[j])/6.;
      ddudw[j][i]=ddudw0_[j][i]+h*(dK1[j]+2.*dK2[j]+2.*dK3[j]+dK4[j])/6.;
    }
    //    if(i==4){
    //      std::cout<<"dudw:dudw0-->"<<dudw[0][i]<<":"<<dudw0[0][i]<<std::endl;
    //      std::cout<<"ddudw:ddudw0-->"<<ddudw[0][i]<<":"<<ddudw0[0][i]<<std::endl;
    //    }
  }
  //    std::cout<<std::endl;

  //      std::cout<<"ddudw:ddudw0-->"<<ddudw[0][i]<<":"<<ddudw0[0][i]<<std::endl;
  //    }
}

//_____________________________________________________________________________
void
RungeKuttaTracker::fACK(int c_use, double *dudz_, double *b,
                        double /* dbdu */ [3][2], double qp0_)
{
  double dx, dy, bx, by, bz, dx2, dy2, dd;
  dx=dudz_[0];
  dy=dudz_[1];
  bx=b[0];
  by=b[1];
  bz=b[2];
  dx2=dx*dx;
  dy2=dy*dy;

  dd=sqrt(1.+dx2+dy2);
  f[0]=dd*(dx*dy*bx -(1.+dx2)*by +dy*bz);
  f[1]=dd*((1.+dy2)*bx -dx*dy*by -dx*bz);
  K[0]=qp0_*f[0];
  K[1]=qp0_*f[1];

  //  std::cout<<"in the fACK"<<std::endl;
  ///A matrix
  ///
  //  double A[2][2]={{0.}};
  //  double C[2][2]={{0.}};
  A[0][0]=  A[0][1]=  A[1][0]=  A[1][1]=  0.;
  ////// a[0][0]=dfx/dx'
  ////// a[0][1]=dfx/dy'
  ////// a[1][0]=dfy/dx'
  ////// a[1][1]=dfy/dy'

  //  A[0][0]=dd*(dy*()
  //((1.+2.*dx2+dy2)*dy*bx - dx*(3.+3.*dx2+2.*dy2)*by + dx*dy*bz)*qp0_/dd;
  A[0][0]=((1.+2.*dx2+dy2)*dy*bx - dx*(3.+3.*dx2+2.*dy2)*by + dx*dy*bz)*qp0_/dd;
  A[0][1]=(dx*(1.+dx2+2.*dy2)*bx - dy*(1.+dx2)*by + (1.+dx2+2.*dy2)*bz)*qp0_/dd;
  A[1][0]=(dx*(1.+dy2)*bx - dy*(1.+2.*dx2+dy2)*by - (1.+2.*dx2+dy2)*bz)*qp0_/dd;
  A[1][1]=(dy*(3.+2.*dx2+3.*dy2)*bx - dx*(1.+dx2+2.*dy2)*by - dx*dy*bz)*qp0_/dd;

#if 0
  std::cout<<"dx:"<<dx<<std::endl;
  std::cout<<"dy:"<<dy<<std::endl;
  std::cout<<"dd:"<<dd<<std::endl;
  std::cout<<"bx:"<<bx<<std::endl;
  std::cout<<"by:"<<by<<std::endl;
  std::cout<<"bz:"<<bz<<std::endl;
  std::cout<<"in the fACK"<<std::endl;
  std::cout<<"A"<<std::endl;
  std::cout<<A[0][0]<<" ,"<<A[0][1]<<std::endl;
  std::cout<<A[1][0]<<" ,"<<A[1][1]<<std::endl;
  std::cout<<"end cal of A"<<std::endl;
#endif

  ///C matrix
  if(c_use==1){
    C[0][0]=0.;
    C[0][1]=0.;
    C[1][0]=0.;
    C[1][1]=0.;
  }
}

//_____________________________________________________________________________
void
RungeKuttaTracker::dKdw0(int c_use, double h,double *dudv0, double *ddudv0,
                         double A1[2][2], double A2[2][2],
                         double A3[2][2], double A4[2][2])
{
  double dKA1[2],dKA2[2],dKA3[2],dKA4[2];
  double dKC1[2],dKC2[2],dKC3[2],dKC4[2];
  double h2=h*h;

#if 0
  std::cout<<"in the dkdw0"<<std::endl;
  std::cout<<"AAAAA"<<std::endl;
  std::cout<<A[0][0]<<" ,"<<A[0][1]<<std::endl;
  std::cout<<A[1][0]<<" ,"<<A[1][1]<<std::endl;
  std::cout<<std::endl;

  std::cout<<"A1"<<std::endl;
  std::cout<<A1[0][0]<<" ,"<<A1[0][1]<<std::endl;
  std::cout<<A1[1][0]<<" ,"<<A1[1][1]<<std::endl;
  std::cout<<std::endl;

  std::cout<<"A2"<<std::endl;
  std::cout<<A2[0][0]<<" ,"<<A2[0][1]<<std::endl;
  std::cout<<A2[1][0]<<" ,"<<A2[1][1]<<std::endl;
  std::cout<<std::endl;

  std::cout<<"A3"<<std::endl;
  std::cout<<A3[0][0]<<" ,"<<A3[0][1]<<std::endl;
  std::cout<<A3[1][0]<<" ,"<<A3[1][1]<<std::endl;
  std::cout<<std::endl;

  std::cout<<"A4"<<std::endl;
  std::cout<<A4[0][0]<<" ,"<<A4[0][1]<<std::endl;
  std::cout<<A4[1][0]<<" ,"<<A4[1][1]<<std::endl;
  std::cout<<std::endl;
#endif

  //dK1/dx0
  for(int i=0;i<2;i++){
    dKA1[i]=0.;
    dKC1[i]=0.;
    for(int j=0;j<2;j++){
      dKA1[i]=dKA1[i]+A1[i][j]*ddudv0[j];
      if(c_use==1) dKC1[i]=dKC1[i]+C1[i][j]*dudv0[j];
    }
    dK1[i]=dKA1[i]+dKC1[i];
  }

  //dK2/dx0
  for(int i=0;i<2;i++){
    dKA2[i]=0.;
    dKC2[i]=0.;
    for(int j=0;j<2;j++){
      dKA2[i]=dKA2[i]+A2[i][j]*(ddudv0[j] + h*dK1[j]/2.);
      if(c_use==1) dKC2[i]=dKC2[i]+C2[i][j]*(dudv0[j]+h*ddudv0[j]/2.+h2*dK1[j]/8.);
    }
    dK2[i]=dKA2[i]+dKC2[i];
  }

  //dK3/dx0
  for(int i=0;i<2;i++){
    dKA3[i]=0.;
    dKC3[i]=0.;
    for(int j=0;j<2;j++){
      dKA3[i]=dKA3[i]+A3[i][j]*(ddudv0[j] + h/2*dK2[j]);
      if(c_use==1) dKC3[i]=dKC3[i]+C3[i][j]*(dudv0[j]+h*ddudv0[j]/2.+h2*dK1[j]/8.);
    }
    dK3[i]=dKA3[i]+dKC3[i];
  }


  //dK4/dx0
  for(int i=0;i<2;i++){
    dKA4[i]=0.;
    dKC4[i]=0.;
    for(int j=0;j<2;j++){
      dKA4[i]=dKA4[i]+A4[i][j]*(ddudv0[j] + h*dK3[j]);
      if(c_use==1) dKC4[i]=dKC4[i]+C4[i][j]*(dudv0[j]+h*ddudv0[j]+h2*dK1[j]/2.);
    }
    dK4[i]=dKA4[i]+dKC4[i];
  }



}

void RungeKuttaTracker::dKdp0(int c_use, double h,double *dudp0, double *ddudp0,
			      double *f1,double *f2,double *f3,double *f4, double A1[2][2],double A2[2][2],double A3[2][2],double A4[2][2]){

  double dKA1[2],dKA2[2],dKA3[2],dKA4[2];
  double dKC1[2],dKC2[2],dKC3[2],dKC4[2];
  double h2=h*h;


  //dK1/dp0
  for(int i=0;i<2;i++){
    dKA1[i]=0.;
    dKC1[i]=0.;
    for(int j=0;j<2;j++){
      dKA1[i]=dKA1[i]+A1[i][j]*ddudp0[j];
      if(c_use==1) dKC1[i]=dKC1[i]+C1[i][j]*dudp0[j];
    }
    dK1[i]=f1[i]+dKA1[i]+dKC1[i];
  }

  //dK2/dp0
  for(int i=0;i<2;i++){
    dKA2[i]=0.;
    dKC2[i]=0.;
    for(int j=0;j<2;j++){
      dKA2[i]=dKA2[i]+A2[i][j]*(ddudp0[j] + h/2*dK1[j]);
      if(c_use==1) dKC2[i] = dKC2[i] + C2[i][j]*(dudp0[j]+h*ddudp0[j]/2.+h2*dK1[j]/8.);
    }
    dK2[i]=f2[i]+dKA2[i]+dKC2[i];
  }

  //dK3/dp0
  for(int i=0;i<2;i++){
    dKA3[i]=dKC3[i]=0.;
    for(int j=0;j<2;j++){
      dKA3[i]=dKA3[i]+A3[i][j]*(ddudp0[j] + h/2*dK2[j]);
      if(c_use==1) dKC3[i]=dKC3[i]+C3[i][j]*(dudp0[j]+h*ddudp0[j]/2.+h2*dK1[j]/8.);
    }
    dK3[i]=f3[i]+dKA3[i]+dKC3[i];
  }


  //dK4/dp0
  for(int i=0;i<2;i++){
    dKA4[i]=dKC4[i]=0.;
    for(int j=0;j<2;j++){
      dKA4[i]=dKA4[i]+A4[i][j]*(ddudp0[j] + h*dK3[j]);
      if(c_use==1) dKC4[i]=dKC4[i]+C4[i][j]*(dudp0[j]+h*ddudp0[j]+h2*dK1[j]/2.);
    }
    dK4[i]=f4[i]+dKA4[i]+dKC4[i];
  }


}

///KMF initiallization
int RungeKuttaTracker::KMFinit1(){
  std::cout<<"---------------------------"<<std::endl;
  std::cout<<"KMFinit1"<<std::endl;
  std::cout<<"---------------------------"<<std::endl;
  for(int i = 0; i < PARASIZE; i++){
    for(int j = 0; j < PARASIZE; j++){
      for(int k= 0; k < PLANESIZE; k++){
	Cmats[i][j][k] = 0.;
      }
    }
  }
  ///initial value of covarnance vector
  Cmats[0][0][0] = 2.0;//1
  Cmats[1][1][0] = 2.0;//1
  Cmats[2][2][0] = 0.05;//0.02
  Cmats[3][3][0] = 0.05;//0.02
  Cmats[4][4][0] = 0.05;//0.02

  return 0;

}

int RungeKuttaTracker::KMFinit2(){
  std::cout<<"---------------------------"<<std::endl;
  std::cout<<"------------KMFinit2------------"<<std::endl;
  std::cout<<"---------------------------"<<std::endl;

  reso  = kmf_.reso[0];
  Cms   = pow(kmf_.mulsth[0],2);
  meas  = kmf_.meas[0];
  measz = kmf_.measz[0];
  zcoor = kmf_.zcoor[0];

  for(int i = 0; i < PARASIZE; i++){
    for(int j = 0; j < PLANESIZE; j++){
      KMFpar0[i][j] = 0.;
      KMFpar1[i][j] = 0.;
      SMTpar[i][j]  = 0.;
    }
  }

  for(int i = 0; i < PARASIZE; i++){
    for(int j = 0; j < PARASIZE; j++){
      for(int k= 0; k < PLANESIZE; k++){
	Fmat[i][j][k]  = 0.;
	Cmat0[i][j][k] = 0.;
	Cinv0[i][j][k] = 0.;
	Cmat1[i][j][k] = 0.;
	Qmat[i][j][k]  = 0.;
      }
    }
  }

  for(int i = 0; i < PLANESIZE; i++){
    KMFchi2[i]    = 0.;
    KMFchi2new[i] = 0.;
    SMTchi2[i]    = 0.;
    SMTchi2new[i] = 0.;
  }

  for(int i = 0; i < PLANESIZE; i++){
    KMFresi[i] = 0.;
    SMTresi[i] = 0.;
  }
  /*
    KMFpar0[0][0] = kmf_.KMFpar0[0][0];
    KMFpar0[1][0] = kmf_.KMFpar0[1][0];
    KMFpar0[2][0] = kmf_.KMFpar0[2][0];
    KMFpar0[3][0] = kmf_.KMFpar0[3][0];
    KMFpar0[4][0] = kmf_.KMFpar0[4][0];
  */
  KMFpar1[0][0] = kmf_.KMFpar1[0][0];
  KMFpar1[1][0] = kmf_.KMFpar1[1][0];
  KMFpar1[2][0] = kmf_.KMFpar1[2][0];
  KMFpar1[3][0] = kmf_.KMFpar1[3][0];
  KMFpar1[4][0] = kmf_.KMFpar1[4][0];
  //  std::cout<<std::endl;
  //  std::cout<<"kmfpar0,par1 : "<<KMFpar0[0][0]<<":"<<KMFpar1[0][0]<<std::endl;
  //  std::cout<<std::endl;
  ///here initial vcalue
  for(int i=0;i<PARASIZE;i++){
    for(int j=0;j<PARASIZE;j++){
      Fmat[i][j][0] = kmf_.Fplane[i][j];
      Cmat1[i][j][0] = Cmats[i][j][0];
      std::cout<<Cmat1[i][j][0]<<", ";
    }
    std::cout<<std::endl;
  }



  if(kmf_.idwi == 1){
    KMFresi[0] =  meas - KMFpar1[0][0];
  }else if(kmf_.idwi == 2){
    KMFresi[0] =  meas - KMFpar1[1][0];
  }

  kmf_.kmfresi[0] =  KMFresi[0];
  Rerr = pow(0.035,2)- Cmat1[0][0][0];
  std::cout<<"Rerr in the init:"<<Rerr<<std::endl;
  if(Rerr < 0){
    kmf_.kmferr[0] = 1.;
  }else{
    kmf_.kmferr[0] = sqrt(Rerr);
  }


  KMFchi2new[0] = 0;
  KMFchi2[0] = 0;

  kmf_.kmfchi2       = KMFchi2[0];
  kmf_.kmfchi2new[0] = KMFchi2new[0];
  kmf_.smtchi2 = 0.;

  return 0;

}



int RungeKuttaTracker::KMFilter()
{

  //  int iplane;
  //  double c30,s30,c45,s45;

  iplane = kmf_.iplane;
  std::cout<<"--------------------------------"<<std::endl;
  std::cout<<"iplane-->"<<iplane<<std::endl;
  std::cout<<"--------------------------------"<<std::endl;

  c30 = cos(30.0*3.14159/180.0);
  s30 = sin(30.0*3.14159/180.0);
  c45 = cos(45.0*3.14159/180.0);
  s45 = sin(45.0*3.14159/180.0);


  // track parameter = 5 vector KMFpar0 or KMFpar1  //
  // x,y,tx,ty,lam //
  // KMFpar0 = par{k|k-1} //
  KMFpar0[0][iplane] = kmf_.KMFpar0[0][iplane];
  KMFpar0[1][iplane] = kmf_.KMFpar0[1][iplane];
  KMFpar0[2][iplane] = kmf_.KMFpar0[2][iplane];
  KMFpar0[3][iplane] = kmf_.KMFpar0[3][iplane];
  KMFpar0[4][iplane] = kmf_.KMFpar0[4][iplane];

  Cms   = pow(kmf_.mulsth[iplane],2);
  meas  = kmf_.meas[iplane];
  measz = kmf_.measz[iplane];
  zcoor = kmf_.zcoor[iplane];
  path  = kmf_.path[iplane];
  idwi  = kmf_.idwi;
  reso  = kmf_.reso[iplane];

  // track model F{k} //
  for(int i=0;i<PARASIZE;i++){
    for(int j=0;j<PARASIZE;j++){
      Fmat[i][j][iplane] = kmf_.Fplane[i][j];
      //	  std::cout<<iplane<<"fmat:-->"<<Fmat[i][j][iplane]<<std::endl;
    }
  }

  ///////////////check Fmat/////
  std::cout<<"7777777777777777777777777777777777"<<std::endl;
  std::cout<<"7777777777777777777777777777777777"<<std::endl;
  double state0[5]={KMFpar0[0][iplane],
		    KMFpar0[1][iplane],
		    KMFpar0[2][iplane],
		    KMFpar0[3][iplane],
		    KMFpar0[4][iplane]};
  double state1[5]={0.};
  for(int i=0;i<PARASIZE;i++){
    for(int j=0;j<PARASIZE;j++){
      state1[i]=state1[i]+Fmat[i][j][iplane]*state0[j];
    }
  }
  std::cout<<"state vector0:state vector1"<<std::endl;
  std::cout<<state0[0]<<" : "<<state1[0]<<std::endl;
  std::cout<<state0[1]<<" : "<<state1[1]<<std::endl;
  std::cout<<state0[2]<<" : "<<state1[2]<<std::endl;
  std::cout<<state0[3]<<" : "<<state1[3]<<std::endl;
  std::cout<<state0[4]<<" : "<<state1[4]<<std::endl;
  std::cout<<"99999999999999999999999999999999999"<<std::endl;

  std::cout<<"cal filter_Q"<<std::endl;
  filter_Q(iplane-1);    // MS error matrix //

  std::cout<<"cal filter_C0"<<std::endl;
  filter_C0(iplane-1);   //` cov matrix C{k|k-1} //

  std::cout<<"cal filter_C0_1"<<std::endl;
  filter_C0_1(iplane-1);   // inverse of cov matrix C{k|k-1} //

  std::cout<<"cal filter_C1"<<std::endl;
  filter_C1(iplane-1);  // cov matrix C{k|k} //

  std::cout<<"cal filter_par1"<<std::endl;
  filter_par1(iplane-1);  //  KMFpar1 = par{k|k} //

  std::cout<<"cal filter_resi"<<std::endl;
  filter_resi(iplane-1); // residual //

  std::cout<<"cal filter_chi2"<<std::endl;
  filter_chi2(iplane-1);  // chi2 //


  return 0;
}

//_____________________________________________________________________________
int RungeKuttaTracker::filter_Q(int ip)
{
  double tx,ty,lam;
  double acont,bcont,ccont;
  double Cxx,Cyy,Cxy;
  iplane1 = ip + 1;
  for(int i = 0; i < PARASIZE; i++){
    for(int j = 0; j < PARASIZE; j++){
      Qmat[i][j][iplane1] = 0.0;
    }
  }

  tx  = KMFpar1[2][iplane0];
  ty  = KMFpar1[3][iplane0];
  lam = KMFpar1[4][iplane0];
  //  std::cout<<tx<<","<<ty<<std::endl;
  acont = 1+pow(tx,2)+pow(ty,2);
  bcont = 1+pow(tx,2);
  ccont = 1+pow(ty,2);

  Cxx = bcont*acont*Cms;
  Cyy = ccont*acont*Cms;
  Cxy = (tx*ty)*acont*Cms;

  // thick //
  Qmat[0][0][iplane1] = Cxx*pow(path,2)/3;
  Qmat[0][1][iplane1] = Cxy*pow(path,2)/3;
  Qmat[0][2][iplane1] = Cxx*pow(path,1)/2;
  Qmat[0][3][iplane1] = Cxy*pow(path,1)/2;

  Qmat[1][0][iplane1] = Qmat[0][1][iplane1];
  Qmat[1][1][iplane1] = Cyy*pow(path,2)/3;
  Qmat[1][2][iplane1] = Cxy*pow(path,1)/2;
  Qmat[1][3][iplane1] = Cyy*pow(path,1)/2;

  Qmat[2][0][iplane1] = Qmat[0][2][iplane1];
  Qmat[2][1][iplane1] = Qmat[1][2][iplane1];
  Qmat[2][2][iplane1] = Cxx;
  Qmat[2][3][iplane1] = Cxy;

  Qmat[3][0][iplane1] = Qmat[0][3][iplane1];
  Qmat[3][1][iplane1] = Qmat[1][3][iplane1];
  Qmat[3][2][iplane1] = Qmat[2][3][iplane1];
  Qmat[3][3][iplane1] = Cyy;

  // thin
  Qmat[2][2][iplane1] =  bcont*acont*Cms;
  Qmat[3][3][iplane1] =  ccont*acont*Cms;
  Qmat[2][3][iplane1] =(tx*ty)*acont*Cms;
  Qmat[3][2][iplane1] = Qmat[2][3][iplane1];
  Qmat[4][4][iplane1] = pow(ty,4)*pow(tx,2)/(acont*bcont)*Cms/pow(lam,2);
  //
  //
  Qmat[4][4][iplane1] = kmf_.eloss[iplane1];
  //
  return 0;

}

//_____________________________________________________________________________
int
RungeKuttaTracker::filter_C0(int ip)
{
  double Emat[PARASIZE][PARASIZE][PLANESIZE];
  iplane1 = ip + 1;
  std::cout<<"filter_C0"<<std::endl;
  std::cout<<"planes:"<<iplane1<<","<<iplane0<<std::endl;
  aba_t(Fmat,Cmat1,Emat,iplane1,iplane0);

  std::cout<<"In the filter_C0"<<std::endl;
  std::cout<<"Fmat"<<std::endl;
  for(int i = 0; i < PARASIZE; i++){
    for(int j = 0; j < PARASIZE; j++){
      std::cout<<Fmat[i][j][iplane1]<<",";
    }
    std::cout<<std::endl;
  }

  std::cout<<"Cmat1"<<std::endl;
  for(int i = 0; i < PARASIZE; i++){
    for(int j = 0; j < PARASIZE; j++){
      std::cout<<Cmat1[i][j][iplane0]<<",";
    }
    std::cout<<std::endl;
  }

  std::cout<<"cmat0"<<std::endl;
  for(int i = 0; i < PARASIZE; i++){
    for(int j = 0; j < PARASIZE; j++){
      Cmat0[i][j][iplane1] = Emat[i][j][iplane1];
      // + Qmat[i][j][iplane1];
      std::cout<<Cmat0[i][j][iplane1]<<",";
    }
    std::cout<<std::endl;
  }

  return 0;

}

//_____________________________________________________________________________
int
RungeKuttaTracker::filter_C0_1(int ip)
{
  //  int gaussj();
  double vect[PARASIZE];
  double Matrix[6][6];
  // double B[PARASIZE][PARASIZE];
  iplane1 = ip + 1;
  //  std::cout<<"cmat0"<<std::endl;
  for(int i = 0; i < PARASIZE;i ++){
    for(int j = 0; j < PARASIZE; j++){
      Matrix[i+1][j+1] = Cmat0[i][j][iplane1];
      //	  std::cout<<Matrix[i+1][j+1]<<",";
    }
    //	std::cout<<std::endl;
    vect[i] = 0.0;
  }
  //  std::cout<<"filter_c0_1, gauss_t"<<std::endl;
  //  MathTools::GaussJordan(Matrix,PARASIZE,vect,PARASIZE)
  gaussj(Matrix,PARASIZE,vect,PARASIZE);
  //  std::cout<<"end filter_c0_1, gauss t"<<std::endl;

  for(int i = 0; i < PARASIZE; i++){
    for(int j = 0; j < PARASIZE; j++){
      Cinv0[i][j][iplane1] = Matrix[i+1][j+1];
    }
  }

  // for test //
  //
  //  for(i = 0; i < PARASIZE; i++){
  //	for(j = 0; j < PARASIZE; j++){
  //	  B[i][j] = 0.0;
  //	  for(k = 0; k < PARASIZE; k++){
  //		B[i][j] += Cinv0[i][k][iplane1]*Cmat0[k][j][iplane1];
  //	  }
  //	}
  //  }
  //
  return 0;
}

//_____________________________________________________________________________
int
RungeKuttaTracker::filter_C1(int ip)
{
  //  int gaussj();
  double vect[PARASIZE];
  double Matrix[6][6];
  // double B[PARASIZE][PARASIZE];
  // double C_[PARASIZE][PARASIZE];

  iplane1 = ip + 1;

  for(int i = 0; i < PARASIZE; i++){
    for(int j = 0; j < PARASIZE; j++){
      Matrix[i+1][j+1] = Cinv0[i][j][iplane1];
    }
  }

  if(idwi == 1){     // x //
    Matrix[1][1] += 1/pow(reso,2);
  }else if(idwi == 2){     // y //
    Matrix[2][2] += 1/pow(reso,2);
  }else if(idwi == 3){     // u //
    Matrix[1][1] += pow(c30,2)/pow(reso,2);
    Matrix[1][2] +=   -c30*s30/pow(reso,2);
    Matrix[2][1] +=   -c30*s30/pow(reso,2);
    Matrix[2][2] += pow(s30,2)/pow(reso,2);
  }else if(idwi == 4){     // v //
    Matrix[1][1] += pow(c30,2)/pow(reso,2);
    Matrix[1][2] +=    c30*s30/pow(reso,2);
    Matrix[2][1] +=    c30*s30/pow(reso,2);
    Matrix[2][2] += pow(s30,2)/pow(reso,2);
  }else if(idwi == 5){     // u //
    Matrix[1][1] += pow(c45,2)/pow(reso,2);
    Matrix[1][2] +=   -c45*s45/pow(reso,2);
    Matrix[2][1] +=   -c45*s45/pow(reso,2);
    Matrix[2][2] += pow(s45,2)/pow(reso,2);
  }else if(idwi == 6){     // v //
    Matrix[1][1] += pow(c45,2)/pow(reso,2);
    Matrix[1][2] +=    c45*s45/pow(reso,2);
    Matrix[2][1] +=    c45*s45/pow(reso,2);
    Matrix[2][2] += pow(s45,2)/pow(reso,2);
  }else{
    printf("1 error in odd or even plane\n");
  }

  // for(int i = 0; i < PARASIZE; i++){
  //   for(int j = 0; j < PARASIZE; j++){
  //     C_[i][j] = Matrix[i+1][j+1];
  //   }
  // }

  gaussj(Matrix,PARASIZE,vect,PARASIZE);

  for(int i = 0; i < PARASIZE; i++){
    for(int j = 0; j < PARASIZE; j++){
      Cmat1[i][j][iplane1] = Matrix[i+1][j+1];
    }
  }

  // for(int i = 0; i < PARASIZE; i++){
  //   for(int j = 0; j < PARASIZE; j++){
  //     B[i][j] = 0.0;
  //     for(int k = 0; k < PARASIZE; k++){
  // 	B[i][j] += C[i][k]*Cmat1[k][j][iplane1];
  //     }
  //   }
  // }

  return 0;
}

//_____________________________________________________________________________
int
RungeKuttaTracker::filter_par1(int ip)
{
  double par[PARASIZE][PLANESIZE];
  double par1[PARASIZE][PLANESIZE];
  //  int MaxV();
  iplane1 = ip + 1;
  for(int i = 0; i < PARASIZE; i++){
    for(int j = 0; j < PLANESIZE; j++){
      par[i][j] = 0.0;
      par1[i][j] = 0.0;
    }
  }

  MaxV(Cinv0,KMFpar0,par,iplane1);

  if(idwi == 1){     // x //
    par[0][iplane1] += meas/pow(reso,2);
  }else if(idwi == 2){     // y //
    par[1][iplane1] += meas/pow(reso,2);
  }else if(idwi == 3){     // u //
    par[0][iplane1] +=  c30*meas/pow(reso,2);
    par[1][iplane1] += -s30*meas/pow(reso,2);
  }else if(idwi == 4){     // v //
    par[0][iplane1] +=  c30*meas/pow(reso,2);
    par[1][iplane1] +=  s30*meas/pow(reso,2);
  }else if(idwi == 5){     // u //
    par[0][iplane1] +=  c45*meas/pow(reso,2);
    par[1][iplane1] += -s45*meas/pow(reso,2);
  }else if(idwi == 6){     // v //
    par[0][iplane1] +=  c45*meas/pow(reso,2);
    par[1][iplane1] +=  s45*meas/pow(reso,2);
  }else{
    printf("2 error in odd or even plane\n");
  }

  MaxV(Cmat1,par,par1,iplane1);

  for(int i = 0; i < PARASIZE; i++){
    KMFpar1[i][iplane1] = par1[i][iplane1];
    kmf_.KMFpar1[i][iplane1] = par1[i][iplane1];
  }

  // if(idwi == 1){
  //   double temp = (meas - KMFpar0[0][iplane1])*Cmat0[0][0][iplane1]
  //     /(Cmat0[0][0][iplane1]+reso)+KMFpar0[0][iplane1];
  // }

  return 0;

}

//_____________________________________________________________________________
int
RungeKuttaTracker::filter_resi(int i)
{
  iplane1 = i + 1;
  if(idwi == 1){     // x //
    KMFresi[iplane1] = meas - KMFpar1[0][iplane1];
  }else if(idwi == 2){     // y //
    KMFresi[iplane1] = meas - KMFpar1[1][iplane1];
  }else if(idwi == 3){     // u //
    KMFresi[iplane1] = meas -
      (c30*KMFpar1[0][iplane1] - s30*KMFpar1[1][iplane1]);
  }else if(idwi == 4){     // v //
    KMFresi[iplane1] = meas -
      (c30*KMFpar1[0][iplane1] + s30*KMFpar1[1][iplane1]);
  }else if(idwi == 5){     // u //
    KMFresi[iplane1] = meas -
      (c45*KMFpar1[0][iplane1] - s45*KMFpar1[1][iplane1]);
  }else if(idwi == 6){     // v //
    KMFresi[iplane1] = meas -
      (c45*KMFpar1[0][iplane1] + s45*KMFpar1[1][iplane1]);
  }else{
    printf("4 error in odd or even plane\n");
  }

  kmf_.kmfresi[iplane1] = KMFresi[iplane1];
  return 0;
}

//_____________________________________________________________________________
int
RungeKuttaTracker::filter_chi2(int i)
{
  iplane1 = i + 1;
  if(idwi == 1){     // x //
    Rerr = pow(reso,2) - Cmat1[0][0][iplane1];
  }else if(idwi == 2){     // y //
    Rerr = pow(reso,2) - Cmat1[1][1][iplane1];
  }else if(idwi == 3){     // u //
    Rerr = pow(reso,2) -
      (pow(c30,2)*Cmat1[0][0][iplane1] +
       pow(s30,2)*Cmat1[1][1][iplane1] -
       c30*s30*(Cmat1[0][1][iplane1] + Cmat1[1][0][iplane1]));
  }else if(idwi == 4){     // v //
    Rerr = pow(reso,2) -
      (pow(c30,2)*Cmat1[0][0][iplane1] +
       pow(s30,2)*Cmat1[1][1][iplane1] +
       c30*s30*(Cmat1[0][1][iplane1] + Cmat1[1][0][iplane1]));
  }else if(idwi == 5){     // u //
    Rerr = pow(reso,2) -
      (pow(c45,2)*Cmat1[0][0][iplane1] +
       pow(s45,2)*Cmat1[1][1][iplane1] -
       c45*s45*(Cmat1[0][1][iplane1] + Cmat1[1][0][iplane1]));
  }else if(idwi == 6){     // v //
    Rerr = pow(reso,2) -
      (pow(c45,2)*Cmat1[0][0][iplane1] +
       pow(s45,2)*Cmat1[1][1][iplane1] +
       c45*s45*(Cmat1[0][1][iplane1] + Cmat1[1][0][iplane1]));
  }else{
    printf("error in odd or even plane\n");
    return -1;
  }

  KMFchi2new[iplane1] = pow(KMFresi[iplane1],2)/Rerr;
  KMFchi2[iplane1]    = KMFchi2[iplane0] + KMFchi2new[iplane1];

  kmf_.kmfchi2             = KMFchi2[iplane1];
  kmf_.kmfchi2new[iplane1] = KMFchi2new[iplane1];
  kmf_.kmferr[iplane1]     = sqrt(Rerr);

  return 0;
}


//INmat1,INmat2,INmat3,OUTmat,plane0,plane1)
int RungeKuttaTracker::ab_tc(
			     double INmat1[PARASIZE][PARASIZE][PLANESIZE],
			     double INmat2[PARASIZE][PARASIZE][PLANESIZE],
			     double INmat3[PARASIZE][PARASIZE][PLANESIZE],
			     double OUTmat[PARASIZE][PARASIZE][PLANESIZE],
			     int plane0,int plane1)
{
  double Xmat[PARASIZE][PARASIZE];

  for(int i=0;i<PARASIZE;i++){
    for(int j=0;j<PARASIZE;j++){
      for(int k=0;k<PLANESIZE;k++){
	OUTmat[i][j][k]=0.;
      }
      Xmat[i][j]=0.;
    }
  }

  for(int i=0;i<PARASIZE;i++){
    for(int j=0;j<PARASIZE;j++){
      for(int k=0;k<PARASIZE;k++){
	Xmat[i][j]=Xmat[i][j]+INmat1[i][k][plane0]*INmat2[j][k][plane1];
      }
    }
  }

  for(int i=0;i<PARASIZE;i++){
    for(int j=0;j<PARASIZE;j++){
      for(int k=0;k<PARASIZE;k++){
	OUTmat[i][j][plane0]=OUTmat[i][j][plane0]+
	  Xmat[i][k]*INmat3[k][j][plane1];
      }
    }
  }
  return 0;


}
//INmat1,INmat2,OUTmat,plane1,plane0)
int RungeKuttaTracker::aba_t(
			     double INmat1[PARASIZE][PARASIZE][PLANESIZE],
			     double INmat2[PARASIZE][PARASIZE][PLANESIZE],
			     double OUTmat[PARASIZE][PARASIZE][PLANESIZE],
			     int plane1,int plane0)
{

  double Xmat[PARASIZE][PARASIZE];

  for(int i=0;i<PARASIZE;i++){
    for(int j=0;j<PARASIZE;j++){
      for(int k=0;k<PLANESIZE;k++){
	OUTmat[i][j][k]=0.;
      }
      Xmat[i][j]=0.;
    }
  }

  for(int i=0;i<PARASIZE;i++){
    for(int j=0;j<PARASIZE;j++){
      for(int k=0;k<PARASIZE;k++){
	Xmat[i][j]=Xmat[i][j]+INmat1[i][k][plane1]*INmat2[k][j][plane0];
	/*	std::cout<<"xmat:"<<Xmat[i][j]<<std::endl;
		std::cout<<"inmat1*inmat2:"<<INmat1[i][k][plane1]*INmat2[k][j][plane0]<<std::endl;
		std::cout<<"inmat1:inmat2:"<<INmat1[i][k][plane1]<<":"<<INmat2[k][j][plane0]<<std::endl;
		std::cout<<"Cmat1:"<<Cmat1[k][j][0]<<std::endl;
		std::cout<<"plane0:"<<plane0<<std::endl;
	*/

	//		std::cout<<"inmat1:"<<INmat1[i][k][plane1]<<std::endl;
	//		std::cout<<"inmat2:"<<INmat2[k][j][plane1]<<std::endl;
      }
    }
  }

  for(int i=0;i<PARASIZE;i++){
    for(int j=0;j<PARASIZE;j++){
      for(int k=0;k<PARASIZE;k++){
	OUTmat[i][j][plane1]=OUTmat[i][j][plane1]
	  +Xmat[i][k]*INmat1[j][k][plane1];
	//		std::cout<<"OUTmat:"<<OUTmat[i][j][plane1]<<std::endl;
      }
    }
  }
  return 0;
}

//_____________________________________________________________________________
//INmat,INvec,OUTvec,plane)
int
RungeKuttaTracker::MaxV(double INmat[PARASIZE][PARASIZE][PLANESIZE],
                        double INvec[PARASIZE][PLANESIZE],
                        double OUTvec[PARASIZE][PLANESIZE],
                        int plane)
{
  for(int i=0;i<PARASIZE;i++){
    for(int k=0;k<PLANESIZE;k++){
      OUTvec[i][k]=0.;
    }
  }

  for(int i=0;i<PARASIZE;i++){
    for(int k=0;k<PARASIZE;k++){
      OUTvec[i][plane]=OUTvec[i][plane]
	+INmat[i][k][plane]*INvec[k][plane];
    }
  }
  return 0;
}

//_____________________________________________________________________________
int
RungeKuttaTracker::fltosm()
{
  iplane=kmf_.iplane;
  //    std::cout<<"in the fltosm iplane:"<<iplane<<std::endl;
  for(int i = 0; i < PARASIZE; i++){
    for(int j = 0; j < PARASIZE; j++){
      Cmats[i][j][iplane] = Cmat1[i][j][iplane];
    }
    SMTpar[i][iplane] = kmf_.KMFpar1[i][iplane];
    kmf_.SMTpar[i][iplane] = kmf_.KMFpar1[i][iplane];
    //        std::cout<<"int fltosm:"<<kmf_.SMTpar[i][iplane]<<", "<<KMFpar1[i][iplane]<<std::endl;
  }

  kmf_.smtchi2new[iplane]  = KMFchi2new[iplane];
  kmf_.smtresi[iplane]     = KMFresi[iplane];
  kmf_.smterr[iplane]      = kmf_.kmferr[iplane];

  SMTchi2[iplane]          = kmf_.smtchi2 ;
  SMTchi2new[iplane]       = kmf_.smtchi2new[iplane];
  SMTresi[iplane]          = kmf_.smtresi[iplane];

  return 0;
}

//_____________________________________________________________________________
int
RungeKuttaTracker::KMsmooth()
{
  iplane = kmf_.iplane;

  reso  = kmf_.reso[iplane];
  meas  = kmf_.meas[iplane];
  idwi  = kmf_.idwi;
  iplane0 = iplane;
  iplane1 = iplane + 1;

  smoother_A();
  smoother_pars();
  smoother_Cs();
  smoother_resi();
  smoother_chi2();

  return 0;
}

//_____________________________________________________________________________
int
RungeKuttaTracker::smoother_A()
{
  ab_tc(Cmat1,Fmat,Cinv0,Amat,iplane0,iplane1);
  return 0;
}

//_____________________________________________________________________________
int
RungeKuttaTracker::smoother_pars()
{
  double par0[PARASIZE][PLANESIZE];
  double par1[PARASIZE][PLANESIZE];

  for(int i = 0; i < PARASIZE; i++){
    for(int j = 0; j < PLANESIZE; j++){
      par0[i][j] = 0.;
      par1[i][j] = 0.;
    }
  }

  for(int i = 0; i < PARASIZE; i++){
    par0[i][iplane0] = SMTpar[i][iplane1] - KMFpar0[i][iplane1];

  }

  MaxV(Amat,par0,par1,iplane0);
  //  std::cout<<"in the smoother_pars:"<<iplane0<<std::endl;
  for(int i = 0; i < PARASIZE; i++){
    SMTpar[i][iplane0] = KMFpar1[i][iplane0] + par1[i][iplane0];
    kmf_.SMTpar[i][iplane0] = SMTpar[i][iplane0];
    //    std::cout<<i<<":"<<kmf_.SMTpar[i][iplane0]<<;
  }
  //  std::cout<<std::endl;
  return 0;
}

//_____________________________________________________________________________
int
RungeKuttaTracker::smoother_Cs()
{
  double Dmat[PARASIZE][PARASIZE][PLANESIZE];
  double Emat[PARASIZE][PARASIZE][PLANESIZE];
  for(int i=0;i<PARASIZE;i++){
    for(int j=0;j<PARASIZE;j++){
      for(int k=0;k<PLANESIZE;k++){
	Dmat[i][j][k]=0.;
	Emat[i][j][k]=0.;
      }
    }
  }

  for(int i=0;i<PARASIZE;i++){
    for(int j=0;j<PARASIZE;j++){
      Dmat[i][j][iplane0]=Cmats[i][j][iplane1]-Cmat0[i][j][iplane1];
    }
  }

  aba_t(Amat,Dmat,Emat,iplane0,iplane0);

  for(int i=0;i<PARASIZE;i++){
    for(int j=0;j<PARASIZE;j++){
      Cmats[i][j][iplane0]=Cmat1[i][j][iplane0]+Emat[i][j][iplane0];
    }
  }

  return(0);
}

//_____________________________________________________________________________
int
RungeKuttaTracker::smoother_resi()
{
  if(idwi == 1){     // x //
    SMTresi[iplane0] = meas - SMTpar[0][iplane0];
  }else if(idwi == 2){     /* y */
    SMTresi[iplane0] = meas - SMTpar[1][iplane0];
  }else if(idwi == 3){     /* u */
    SMTresi[iplane0] = meas -
      (c30*SMTpar[0][iplane0] - s30*SMTpar[1][iplane0]);
  }else if(idwi == 4){     /* v */
    SMTresi[iplane0] = meas -
      (c30*SMTpar[0][iplane0] + s30*SMTpar[1][iplane0]);
  }else if(idwi == 5){     /* u */
    SMTresi[iplane0] = meas -
      (c45*SMTpar[0][iplane0] - s45*SMTpar[1][iplane0]);
  }else if(idwi == 6){     /* v */
    SMTresi[iplane0] = meas -
      (c45*SMTpar[0][iplane0] + s45*SMTpar[1][iplane0]);
  }else{
    printf("4 error in odd or even plane\n");
  }

  kmf_.smtresi[iplane0] = SMTresi[iplane0];
  return 0;
}

//_____________________________________________________________________________
int
RungeKuttaTracker::smoother_chi2()
{
  if(idwi == 1){     /* x */
    Rerr = pow(reso,2) - Cmats[0][0][iplane0];
  }else if(idwi == 2){     /* y */
    Rerr = pow(reso,2) - Cmats[1][1][iplane0];
  }else if(idwi == 3){     /* u */
    Rerr = pow(reso,2) -
      (pow(c30,2)*Cmats[0][0][iplane0] +
       pow(s30,2)*Cmats[1][1][iplane0] -
       c30*s30*(Cmats[0][1][iplane0] + Cmats[1][0][iplane0]));
  }else if(idwi == 4){     /* v */
    Rerr = pow(reso,2) -
      (pow(c30,2)*Cmats[0][0][iplane0] +
       pow(s30,2)*Cmats[1][1][iplane0] +
       c30*s30*(Cmats[0][1][iplane0] + Cmats[1][0][iplane0]));
  }else if(idwi == 5){     /* u */
    Rerr = pow(reso,2) -
      (pow(c45,2)*Cmats[0][0][iplane0] +
       pow(s45,2)*Cmats[1][1][iplane0] -
       c45*s45*(Cmats[0][1][iplane0] + Cmats[1][0][iplane0]));
  }else if(idwi == 6){     /* v */
    Rerr = pow(reso,2) -
      (pow(c45,2)*Cmats[0][0][iplane0] +
       pow(s45,2)*Cmats[1][1][iplane0] +
       c45*s45*(Cmats[0][1][iplane0] + Cmats[1][0][iplane0]));
  }else{
    printf("error in odd or even plane\n");
    return -1;
  }

  SMTchi2new[iplane0] = pow(SMTresi[iplane0],2)/Rerr;

  kmf_.smtchi2new[iplane0] = SMTchi2new[iplane0];
  kmf_.smterr[iplane0]     = sqrt(Rerr);

  return 0;

}

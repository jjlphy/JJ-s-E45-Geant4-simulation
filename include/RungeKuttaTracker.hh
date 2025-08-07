// -*- C++ -*-

#ifndef RUNGE_KUTTA_TRACKER_HH
#define RUNGE_KUTTA_TRACKER_HH

#include <G4ThreeVector.hh>

#include "track.hh"
#include "kmf.h"
//#include "filter.h"

//_____________________________________________________________________________
class RungeKuttaTracker
{
public:
  RungeKuttaTracker(int c_use, Track* aTrack);
  ~RungeKuttaTracker();

private:
  void RungeKuttaTracking(int c_use, Track* aTrack);
  void RungeKuttaFieldIntegral(int c_use, double *init_par,double *final_par,
			       G4ThreeVector &B);
  void RungeKuttaFieldIntegral(int c_use, double *init_par,double *final_par,
			       G4ThreeVector &B,
			       G4ThreeVector &dBdY, G4ThreeVector &dBdZ);
  int RungeKuttaFit(int c_use, int iteration, Track* aTrack,double *rk_par0,
		    double *rk_par, double* rk_hit);
  void RungeKuttaStep(int c_use, double qp0,double h,
		      double z0,double *u0,double *dudz0,
		      double dudw0[2][5],double ddudw0[2][5]);
  double z;
  double u[2];
  double dudz[2];
  double dudw[2][5];
  double ddudw[2][5];
  // double qp;

  double z0;
  double u0[2];
  double dudz0[2];
  double dudw0[2][5];
  double ddudw0[2][5];
  double qp0;

  void fACK(int c_use, double *dudz,double *b, double dbdu[3][2], double qp0);
  //matrix ACK
  double A[2][2];
  double C[2][2];
  double K[2];
  double f[2];

  // double A1[2][2],A2[2][2],A3[2][2],A4[2][2];
  double C1[2][2],C2[2][2],C3[2][2],C4[2][2];

  ///calculation of deviation
  void dKdw0(int c_use, double h,double *dudv0, double *ddudv0,
             double A1[2][2], double A2[2][2],
             double A3[2][2], double A4[2][2]);
  void dKdp0(int c_use, double h,double *dudv0, double *ddudv0,
             double *f1, double *f2,double *f3,double *f4,
             double A1[2][2], double A2[2][2],
             double A3[2][2], double A4[2][2]);

  double dK1[2],dK2[2],dK3[2],dK4[2];

  ///rungekutta hits
  // double rk_hits[50][30];//[par][nhit]
  double RKHits[50][22];

  //  F-matrix
  //  double F[5][5];

  // double kmf_par[30][5];

  /////kmf function
  int KMFinit1();
  int KMFinit2();



  //  double Fplane[5][5];

  KMF kmf_;
  double reso,zcoor,Cms,meas,measz, Rerr, path;
  int idwi;
  double pGeV;
  double c30,s30,c45,s45;

  int filter_Q(int ip);
  int filter_C0(int ip);//FC?? strange
  // int filter_C0_v2(int i);//FCFt
  int filter_C0_1(int ip);
  int filter_C1(int ip);
  int filter_par1(int ip);
  int filter_resi(int i);
  int filter_chi2(int i);

  int KMFilter();

  ///libs
#define PARASIZE 5
#define PLANESIZE 20
#define RESISIZE 2

  double KMFpar0[PARASIZE][PLANESIZE];
  double KMFpar1[PARASIZE][PLANESIZE];
  double Fmat[PARASIZE][PARASIZE][PLANESIZE];
  double Cmat0[PARASIZE][PARASIZE][PLANESIZE];
  double Cinv0[PARASIZE][PARASIZE][PLANESIZE];
  double Cmat1[PARASIZE][PARASIZE][PLANESIZE];
  double Qmat[PARASIZE][PARASIZE][PLANESIZE];
  double KMFresi[PLANESIZE];
  double KMFchi2[PLANESIZE];
  double KMFchi2new[PLANESIZE];

  double SMTpar[PARASIZE][PLANESIZE];
  double Amat[PARASIZE][PARASIZE][PLANESIZE];
  double Cmats[PARASIZE][PARASIZE][PLANESIZE];
  double SMTresi[PLANESIZE];
  double SMTchi2[PLANESIZE];
  double SMTchi2new[PLANESIZE];

  int iplane;


  int ab_tc(double INmat1[PARASIZE][PARASIZE][PLANESIZE],
            double INmat2[PARASIZE][PARASIZE][PLANESIZE],
            double INmat3[PARASIZE][PARASIZE][PLANESIZE],
            double OUTmat[PARASIZE][PARASIZE][PLANESIZE],
            int plane0, int plane1);
  int aba_t(double INmat1[PARASIZE][PARASIZE][PLANESIZE],
            double INmat2[PARASIZE][PARASIZE][PLANESIZE],
            double OUTmat[PARASIZE][PARASIZE][PLANESIZE],
            int plane1,int plane0);
  int MaxV(double INmat[PARASIZE][PARASIZE][PLANESIZE],
           double INvec[PARASIZE][PLANESIZE],
           double OUTvec[PARASIZE][PLANESIZE],
           int plane);

  //  int gaussj(double a[6][6],int n, double *b, int m);

  int iplane0;
  int iplane1;
  int fltosm();
  int KMsmooth();
  int smoother_A();
  int smoother_pars();
  int smoother_Cs();
  int smoother_resi();
  int smoother_chi2();

};
#endif

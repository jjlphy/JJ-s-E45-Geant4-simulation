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
/*
int filter_Q();
int filter_C0();
int filter_C0_1();
int filter_C1();
int filter_par1();
int filter_resi();
int filter_chi2();

int KMFilter();

///libs
int ab_tc(double INmat1[PARASIZE][PARASIZE][PLANESIZE],
          double INmat2[PARASIZE][PARASIZE][PLANESIZE],
          double INmat3[PARASIZE][PARASIZE][PLANESIZE],
          double OUTmat[PARASIZE][PARASIZE][PLANESIZE],
          int plane0, int plane1);
int aba_t(double INmat1[PARASIZE][PARASIZE][PLANESIZE],
          double INmat2[PARASIZE][PARASIZE][PLANESIZE],
          double OUTmat[PARASIZE][PARASIZE][PLANESIZE],
          int plane0,int plane1);
int MaxV(double INmat[PARASIZE][PARASIZE][PLANESIZE],
         double INvec[PARASIZE][PLANESIZE],
         double OUTvec[PARASIZE][PLANESIZE],
         int plane);
*/

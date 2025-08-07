/* MINUIT functions
 *						by tsugu
 */

#ifndef MINUIT_H
#define MINUIT_H

#ifndef MAX_CH_LEN
#define MAX_CH_LEN 80
#endif /* MAX_CH_LEN */

#ifndef FORTRAN_DEFAULT_CHAR
#define FORTRAN_DEFAULT_CHAR ' '
#endif /* FORTRAN_DEFAULT_CHAR */
#define MINUIT_VAR_LEN 10

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* function type defined by users */
typedef void ( * MinuitFUtil )();
typedef void ( * MinuitFCN )( int*,
	double*, double*, double*, int*, MinuitFUtil );

/* Functions call original CERN Library functions */
void MINTIO( int iread, int iwrite, int isave );
void MNCOMD( MinuitFCN fcn, char* ichstr, int* icondn, MinuitFUtil futil );
void MNCONT( MinuitFCN fcn, int num1, int num2, int npt,
	double* xpt, double* ypt, int* nfound, MinuitFUtil futil );
void MNEMAT( double* emat, int ndim );
void MNERRS( int num, double* eplus,
	double* eminus, double* eparab, double* globcc );
void MNEXCM( MinuitFCN fcn, char* ichcom,
	double* arglis, int narg, int* ierflg, MinuitFUtil futil );
void MNINIT( int ird, int iwr, int isav );
void MNINPU( int nunit, int* ierr );
void MNINTR( MinuitFCN fcn, MinuitFUtil futil );
void MNPARM( int num, char* ichnam,
	double stval, double step, double bnd1, double bnd2, int* ierflg );
void MNPARS( char* ichstr, int* icondn );
void MNPOUT( int num, char* chnam, double* val,
	double* error, double* bnd1, double* bnd2, int* ivarbl );
void MNSETI( char* ictitle );
void MNSTAT( double* fmin, double* fedm,
	double* errdef, int* npari, int* nparx, int* istat );

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* MINUIT_H */

/* End of file */


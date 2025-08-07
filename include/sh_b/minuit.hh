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


#endif /* __cplusplus */

#endif /* MINUIT_H */

/* End of file */


#ifndef _DRIVERS_H_
#define _DRIVERS_H_ 
/*
   ------------------------------------------------------------- 
   File drivers.h of ADOL-C version 1.8.7      as of Feb/28/2000
   ------------------------------------------------------------- 
   Easy to use drivers for optimization and nonlinear equations
   (with C and C++ callable interfaces including Fortran 
    callable versions).

   Last changes:
    20000228  olvo:    corrected comment at lagra_hess_vec
      981130  olvo:    newly created from adutilsc.h

   ------------------------------------------------------------- 
*/

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/
/*                                                      PUBLIC EXPORTS ONLY */

#ifdef __cplusplus
/****************************************************************************/
/****************************************************************************/
/*                                                            No C++ THINGS */


/****************************************************************************/
/****************************************************************************/
/*                                                         Now the C THINGS */
extern "C" {
#endif


/****************************************************************************/
/*                         DRIVERS FOR OPTIMIZATION AND NONLINEAR EQUATIONS */

/*--------------------------------------------------------------------------*/
/*                                                                 function */
/* function(tag, m, n, x[n], y[m])                                          */

int function(short,int,int,double*,double*);

fint function_(fint*,fint*,fint*,fdouble*,fdouble*);


/*--------------------------------------------------------------------------*/
/*                                                                 gradient */
/* gradient(tag, n, x[n], g[n])                                             */

int gradient(short,int,double*,double*);

fint gradient_(fint*,fint*,fdouble*,fdouble*);


/*--------------------------------------------------------------------------*/
/*                                                                 jacobian */
/* jacobian(tag, m, n, x[n], J[m][n])                                       */

int jacobian(short,int,int,double*,double**);

fint jacobian_(fint*,fint*,fint*,fdouble*,fdouble*);


/*--------------------------------------------------------------------------*/
/*                                                         vector_jacobian  */
/* vec_jac(tag, m, n, repeat, x[n], u[m], v[n])                             */

int vec_jac(short,int,int,int,double*,double*,double*);

fint vec_jac_(fint*,fint*,fint*,fint*,fdouble*,fdouble*,fdouble*);


/*--------------------------------------------------------------------------*/
/*                                                          jacobian_vector */
/* jac_vec(tag, m, n, x[n], v[n], u[m]);                                    */

int jac_vec(short,int,int,double*,double*,double*);

fint jac_vec_(fint*,fint*,fint*,fdouble*,fdouble*,fdouble*);


/*--------------------------------------------------------------------------*/
/*                                                                  hessian */
/* hessian(tag, n, x[n], lower triangle of H[n][n])                         */

int hessian(short,int,double*,double**);

fint hessian_(fint*,fint*,fdouble*,fdouble*);


/*--------------------------------------------------------------------------*/
/*                                                           hessian_vector */
/* hess_vec(tag, n, x[n], v[n], w[n])                                       */

int hess_vec(short,int,double*,double*,double*);

fint hess_vec_(fint*,fint*,fdouble*,fdouble*,fdouble*);


/*--------------------------------------------------------------------------*/
/*                                                  lagrange_hessian_vector */
/* lagra_hess_vec(tag, m, n, x[n], v[n], u[m], w[n])                        */

int lagra_hess_vec(short,int,int,double*,double*,double*,double*);

fint lagra_hess_vec_(fint*,fint*,fint*,fdouble*,fdouble*,fdouble*,fdouble*);


/****************************************************************************/
/*                                                               THAT'S ALL */
#ifdef __cplusplus
}
#endif

#endif


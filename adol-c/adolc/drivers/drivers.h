/*---------------------------------------------------------------------------- 
 ADOL-C -- Automatic Differentiation by Overloading in C++
 File:     drivers/drivers.h
 Revision: $Id: drivers.h,v 1.2 2004/10/14 13:29:48 e_arnold Exp $
 Contents: Easy to use drivers for optimization and nonlinear equations
           (with C and C++ callable interfaces including Fortran 
            callable versions).

 Copyright (c) 2004
               Technical University Dresden
               Department of Mathematics
               Institute of Scientific Computing
  
 This file is part of ADOL-C. This software is provided under the terms of
 the Common Public License. Any use, reproduction, or distribution of the
 software constitutes recipient's acceptance of the terms of this license.
 See the accompanying copy of the Common Public License for more details.

 History:
          20040416 kowarz: adapted to configure - make - make install
          20030303 andrea: new hess_mat(..), new hessian2(..)
          20000228 olvo:   corrected comment at lagra_hess_vec
          19981130 olvo:   newly created from adutilsc.h

----------------------------------------------------------------------------*/
#if !defined(ADOLC_DRIVERS_DRIVERS_H)
#define ADOLC_DRIVERS_DRIVERS_H 1

#include "../common.h"

BEGIN_C_DECLS

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
/*                                                           jacobian_partx */
/* jacobian_partx(tag, m, n, x[n][], J[m][n][])                             */

int jacobian_partx(short,int,int,int*,double**,double***);
fint jacobian_partx_(fint*,fint*,fint*,fint*,fdouble**,fdouble***);

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
/* uses Hessian-vector product                                              */
int hessian(short,int,double*,double**);
fint hessian_(fint*,fint*,fdouble*,fdouble*);

/*--------------------------------------------------------------------------*/
/*                                                                 hessian2 */
/* hessian2(tag, n, x[n], lower triangle of H[n][n])                        */
/* uses Hessian-matrix product                                              */
int hessian2(short,int,double*,double**);
fint hessian2_(fint*,fint*,fdouble*,fdouble*);

/*--------------------------------------------------------------------------*/
/*                                                           hessian_vector */
/* hess_vec(tag, n, x[n], v[n], w[n])                                       */
int hess_vec(short,int,double*,double*,double*);
fint hess_vec_(fint*,fint*,fdouble*,fdouble*,fdouble*);

/*--------------------------------------------------------------------------*/
/*                                                           hessian_matrix */
/* hess_mat(tag, n, q, x[n], V[n][q], W[q][n])                              */
int hess_mat(short,int,int,double*,double**,double**);
fint hess_mat_(fint*,fint*,fint*,fdouble*,fdouble**,fdouble**);

/*--------------------------------------------------------------------------*/
/*                                                  lagrange_hessian_vector */
/* lagra_hess_vec(tag, m, n, x[n], v[n], u[m], w[n])                        */
int lagra_hess_vec(short,int,int,double*,double*,double*,double*);
fint lagra_hess_vec_(fint*,fint*,fint*,fdouble*,fdouble*,fdouble*,fdouble*);

END_C_DECLS

/****************************************************************************/
#endif


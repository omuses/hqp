/*---------------------------------------------------------------------------- 
 ADOL-C -- Automatic Differentiation by Overloading in C++
 File:     drivers/odedrivers.h
 Revision: $Id: odedrivers.h,v 1.1 2004/10/13 14:18:12 e_arnold Exp $
 Contents: Easy to use drivers for ordinary differential equations (ODE)
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
          19981201 olvo:   newly created from adutils.h & adutilsc.h

----------------------------------------------------------------------------*/
#if !defined(ADOLC_DRIVERS_ODEDRIVERS_H)
#define ADOLC_DRIVERS_ODEDRIVERS_H 1

#include "../common.h"
#include "../interfaces.h"

BEGIN_C_DECLS

/****************************************************************************/
/*                                                         DRIVERS FOR ODEs */

/*--------------------------------------------------------------------------*/
/*                                                                  forodec */
/* forodec(tag, n, tau, dold, dnew, X[n][d+1])                              */
int forodec(short,int,double,int,int,double**);
fint forodec_(fint*,fint*,fdouble*,fint*,fint*,fdouble*);

/*--------------------------------------------------------------------------*/
/*                                                                  accodec */
/* accodec(n, tau, d, Z[n][n][d+1], B[n][n][d+1], nz[n][n])                 */
void accodec(int,double,int,double***,double***,short**);
fint accodec_(fint*,fdouble*,fint*,fdouble*,fdouble*);

END_C_DECLS

/****************************************************************************/
/****************************************************************************/
/*                                                       Now the C++ THINGS */
#if defined(__cplusplus)

/****************************************************************************/
/*                                       DRIVERS FOR ODEs, overloaded calls */

/*--------------------------------------------------------------------------*/
/*                                                                   forode */
/* forode(tag, n, tau, dold, dnew, X[n][d+1])                               */
inline int forode(
            short  tag,            // tape identifier
            int    n,              // space dimension
            double tau,            // scaling
            int    dold,           // previous degree defaults to zero
            int    dnew,           // New degree of consistency
            double **X)            // Taylor series
{
  return forodec(tag,n,tau,dold,dnew,X);
}

/*--------------------------------------------------------------------------*/
/*                                                                   forode */
/*        the scaling tau defaults to 1                                     */
/*                                                                          */
/*  forode(tag, n, dold, dnew, X[n][d+1])                                   */
inline int forode(short tag, int n, int dold, int dnew, double** X)
{
  return forodec(tag,n,1.0,dold,dnew,X);
}

/*--------------------------------------------------------------------------*/
/*                                                                   forode */
/*        previous order defaults to 0                                      */
/*                                                                          */
/* forode(tag, n, tau, dnew, X[n][d+1])                                     */
inline int forode( short tag, int n, double tau, int deg, double **X)
{
  return  forodec(tag,n,tau,0,deg, X);
}

/*--------------------------------------------------------------------------*/
/*                                                                   forode */
/*        both tau and dold default                                         */
/*                                                                          */
/* forode(tag, n, dnew, X[n][d+1])                                          */
inline int forode(short tag, int n, int deg, double** X)
{
  return  forode(tag,n,1.0,0,deg,X);
}

/*--------------------------------------------------------------------------*/
/*                                                                   accode */
/* accode(n, tau, d, Z[n][n][d+1], B[n][n][d+1], nz[n][n])                  */
inline void accode(
            int    n,               // space dimension
            double tau,             // scaling defaults to 1.0
            int    deg,             // highest degree
            double ***A,            // input tensor of "partial" Jacobians
            double ***B,            // output tensor of "total" Jacobians
            short  **nonzero = 0)   // optional sparsity characterization 
{
  accodec(n,tau,deg,A,B,nonzero);
}

/*--------------------------------------------------------------------------*/
/*                                                                   accode */
/*       scaling defaults to 1                                              */
/*                                                                          */
/* accode(n, d, Z[n][n][d+1], B[n][n][d+1], nz[n][n])                       */
inline void accode(
            int    n,             // space dimension
            int    deg,           // highest degree
            double ***A,          // input tensor of "partial" Jacobians
            double ***B,          // output tensor of "total" Jacobians
            short  **nonzero = 0) // optional sparsity characterization 
{
  accodec(n,1.0,deg,A,B,nonzero);
}

#endif

/****************************************************************************/
#endif


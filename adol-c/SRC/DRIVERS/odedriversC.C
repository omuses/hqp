#define _ODEDRIVERSC_CPP_
#define _ADOLC_SRC_
/*
   ------------------------------------------------------------------------
   File odedriversC.C of ADOL-C version 1.8.0               as of Nov/30/98
   ------------------------------------------------------------------------
   Easy to use drivers for ordinary differential equations (ODE)
   (Implementation of C++ callable interfaces).

   Last changes:
      981130  olvo:   newly created from odedrivers.c

  ------------------------------------------------------------------------
*/


/****************************************************************************/
/*                                                                 INCLUDES */
#include "dvlparms.h" /* Developers Parameters */
#include "usrparms.h" /* Users Parameters      */
#include "odedrivers.h"


/****************************************************************************/
/*                                       DRIVERS FOR ODEs, overloaded calls */

/*--------------------------------------------------------------------------*/
/*                                                                   forode */
/*                                                                          */
int forode( short  tag,            // tape identifier
            int    n,              // space dimension
            double tau,            // scaling defaults to 1.0
            int    dol,            // previous degree defaults to zero
            int    deg,            // New degree of consistency
            double **y)            // Taylor series
/* forode(tag, n, tau, dold, dnew, X[n][d+1])                               */
{
  return forodec(tag,n,tau,dol,deg,y);
}


/*--------------------------------------------------------------------------*/
/*                                                                   forode */
/*        the scaling tau defaults to 1                                     */
/*                                                                          */
int forode( short tag, int n, double tau, int deg, double **y)
/* forode(tag, n, tau, deg, X[n][d+1])                                      */
{
  /***   Default for previous degree is zero, do things from scratch ***/
  int zero = 0;
  return  forodec(tag, n, tau, zero, deg, y);
}


/*--------------------------------------------------------------------------*/
/*                                                                   forode */
/*        previous order defaults to 0                                      */
/*                                                                          */
int forode(short tag, int n, int dol, int deg, double** y)
/* forode(tag, n, dold, dnew, X[n][d+1])                                    */
{
  /***   Default for scaling is 1.0       *****/
  double tau = 1.0;
  return forodec(tag, n, tau, dol, deg, y);
}


/*--------------------------------------------------------------------------*/
/*                                                                   forode */
/*        both tau and dold default                                         */
/*                                                                          */
int forode(short tag, int n, int deg, double** y)
/* forode(tag, n, deg, X[n][d+1])                                           */
{
  /***    Combination of both previous defaults   ******/
  double tau = 1.0;
  return  forode(tag, n, tau, deg, y);
}


/*--------------------------------------------------------------------------*/
/*                                                                   accode */
/*                                                                          */
void accode(int    n,           // space dimension
            double tau,         // scaling defaults to 1.0
            int    deg,         // highest degree
	    double ***A,        // input tensor of "partial" Jacobians
            double ***B,        // output tensor of "total" Jacobians
	    short  **nonzero )  // optional sparsity characterization 
/* accode(n, tau, d, Z[n][n][d+1], B[n][n][d+1], nz[n][n])                  */
{
  accodec(n,tau,deg,A,B,nonzero);
}


/*--------------------------------------------------------------------------*/
/*                                                                   accode */
/*       scaling defaults to 1                                              */
void accode(int    n,            // space dimension
            int    deg,          // highest degree
	    double ***A,         // input tensor of "partial" Jacobians
            double ***B,         // output tensor of "total" Jacobians
	    short  **nonzero )   // optional sparsity characterization 
/* accode(n, d, Z[n][n][d+1], B[n][n][d+1], nz[n][n])                       */
{

  double tau = 1.0;
  accodec(n,tau,deg,A,B,nonzero);
}


/****************************************************************************/
/*                                                               THAT'S ALL */
#undef _ADOLC_SRC_
#undef _ODEDRIVERSC_CPP_

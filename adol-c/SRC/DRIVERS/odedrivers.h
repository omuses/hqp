#ifndef _ODEDRIVERS_H_
#define _ODEDRIVERS_H_ 
/*
   ------------------------------------------------------------- 
   File odedrivers.h of ADOL-C version 1.8.0     as of Dec/01/98
   ------------------------------------------------------------- 
   Easy to use drivers for ordinary differential equations (ODE)
   (with C and C++ callable interfaces including Fortran 
    callable versions).

   Last changes:
      981201  olvo:   newly created from adutils.h & adutilsc.h

   ------------------------------------------------------------- 
*/

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/
/*                                                      PUBLIC EXPORTS ONLY */

#ifdef __cplusplus
/****************************************************************************/
/****************************************************************************/
/*                                                       Now the C++ THINGS */


/****************************************************************************/
/*                                       DRIVERS FOR ODEs, overloaded calls */

/*--------------------------------------------------------------------------*/
/*                                                                   forode */
/* forode(tag, n, tau, dold, dnew, X[n][d+1])                               */

int forode(short, int, double, int, int, double**);

/*--------------------------------------------------------------------------*/
/*                                                                   forode */
/*        the scaling tau defaults to 1                                     */
/*                                                                          */
/*  forode(tag, n, dold, dnew, X[n][d+1])                                   */

int forode(short, int, /*1.0,*/ int, int, double**);

/*--------------------------------------------------------------------------*/
/*                                                                   forode */
/*        previous order defaults to 0                                      */
/*                                                                          */
/* forode(tag, n, tau, dnew, X[n][d+1])                                     */

int forode(short, int, double, /*0,*/ int, double**);

/*--------------------------------------------------------------------------*/
/*                                                                   forode */
/*        both tau and dold default                                         */
/*                                                                          */
/* forode(tag, n, dnew, X[n][d+1])                                          */

int forode(short, int, /*1.0 , 0,*/ int, double**);

/*--------------------------------------------------------------------------*/
/*                                                                   accode */
/* accode(n, tau, d, Z[n][n][d+1], B[n][n][d+1], nz[n][n])                  */

void accode(int, double, int, double***, double***, short** = 0 );

/*--------------------------------------------------------------------------*/
/*                                                                   accode */
/*       scaling defaults to 1                                              */
/*                                                                          */
/* accode(n, d, Z[n][n][d+1], B[n][n][d+1], nz[n][n])                       */

void accode(int, /*1.0,*/ int, double***, double***, short** = 0 );


/****************************************************************************/
/****************************************************************************/
/*                                                         Now the C THINGS */
extern "C" {
#endif


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


/****************************************************************************/
/*                                                               THAT'S ALL */
#ifdef __cplusplus
}
#endif

/*--------------------------------------------------------------------------*/
/*                                                        automatic include */
#include "interfaces.h"

#endif


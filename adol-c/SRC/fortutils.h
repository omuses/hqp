#ifndef _FORTUTILS_H_
#define _FORTUTILS_H_
/*
   ----------------------------------------------------------------
   File fortutils.h of ADOL-C version 1.8.0         as of Nov/30/98
   ----------------------------------------------------------------
   Internal tools to handle Fortran arrays

   Last changed : 
        981130 olvo   newly created from driversc.c

   ----------------------------------------------------------------
*/


/****************************************************************************/
/****************************************************************************/
/****************************************************************************/
/*                                              ADOLC INTERNAL EXPORTS ONLY */
#ifdef _ADOLC_SRC_

#include "usrparms.h"

#ifdef __cplusplus
/****************************************************************************/
/****************************************************************************/
/*                                                            No C++ THINGS */


/****************************************************************************/
/****************************************************************************/
/*                                                         Now the C THINGS */
extern "C" {
#endif

void spread1(int m, fdouble* x, double* X);
void pack1(int m, double* X, fdouble* x);

void spread2(int m, int n, fdouble* x, double** X);
void pack2(int m, int n, double** X, fdouble* x);

void spread3(int m, int n, int p, fdouble* x, double*** X);
void pack3(int m, int n, int p, double*** X, fdouble* x);


/****************************************************************************/
/*                                                               THAT'S ALL */
#ifdef __cplusplus
}
#endif

#endif /* _ADOLC_SRC_ */

#endif

#ifndef _LU_H
#define _LU_H
/*
   --------------------------------------------------------------
   File LU.C of ADOL-C version 1.8.5              as of Nov/22/99
   --------------------------------------------------------------

   Example:  'active' LU-decomposition and according solver

   Last changes:
     990922 olvo    first version 

   --------------------------------------------------------------
*/

/****************************************************************************/
/*                                                                 INCLUDES */
#include "../../../adolc/adolc.h"           // use of ALL ADOL-C interfaces


/****************************************************************************/
/* Simple LU-factorization according to Crout's algorithm without pivoting */
void LUfact(int n, adouble **A);


/****************************************************************************/
/* Solution of A*x=b by forward and backward substitution */
void LUsolve(int n, adouble **A, adouble *bx);


/****************************************************************************/
/*                                                              END OF FILE */
#endif

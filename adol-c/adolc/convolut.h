/*---------------------------------------------------------------------------- 
 ADOL-C -- Automatic Differentiation by Overloading in C++
 File:     convolut.h
 Revision: $Id: convolut.h,v 1.2 2004/10/14 13:29:47 e_arnold Exp $
 Contents: Convolution routines (used by ho_rev.mc)

 Copyright (c) 2004
               Technical University Dresden
               Department of Mathematics
               Institute of Scientific Computing
  
 This file is part of ADOL-C. This software is provided under the terms of
 the Common Public License. Any use, reproduction, or distribution of the
 software constitutes recipient's acceptance of the terms of this license.
 See the accompanying copy of the Common Public License for more details.

 History:
          20040423 kowarz: adapted to configure - make - make install
          19981130 olvo:   last check
          19980707 olvo:   created this file from parts of adutilsc.h
          19980616 olvo:   (1) void copyAndZeroset(..)
                           (2) void inconv0(..)
                               void deconv0(..)
 
----------------------------------------------------------------------------*/

#if !defined(ADOLC_CONVOLUT_H)
#define ADOLC_CONVOLUT_H 1

#include "common.h"

BEGIN_C_DECLS

/****************************************************************************/
/*                                                              CONVOLUTION */

/*--------------------------------------------------------------------------*/
/* Evaluates convolution of a and b to c */
void conv( int dim, double *a, double *b, double *c );

/****************************************************************************/
/*                                                  INCREMENTAL CONVOLUTION */

/*--------------------------------------------------------------------------*/
/* Increments truncated convolution of a and b to c */
void inconv ( int dim, double *a, double *b, double* c );

/*--------------------------------------------------------------------------*/
/* Increments truncated convolution of a and b to c and sets a to zero */
void inconv0( int dim, double *a, double *b, double* c );


/****************************************************************************/
/*                                                  DECREMENTAL CONVOLUTION */

/*--------------------------------------------------------------------------*/
/* Decrements truncated convolution of a and b to c */
void deconv ( int dim, double* a, double *b, double* c );

/*--------------------------------------------------------------------------*/
/* Decrements truncated convolution of a and b to c and sets a to zero */
void deconv0( int dim, double* a, double *b, double* c );


/****************************************************************************/
/*                                                    OTHER USEFUL ROUTINES */

/*--------------------------------------------------------------------------*/
void divide(int dim, double* a, double *b, double* c);

/*--------------------------------------------------------------------------*/
void recipr(int dim, double  a, double *b, double* c);


/****************************************************************************/
/*                                                                  ZEROING */

/*--------------------------------------------------------------------------*/
/* Set a to zero */
void zeroset(int dim, double* a);

/*--------------------------------------------------------------------------*/
/* Copies a to tmp and initializes a to zero */
void copyAndZeroset( int dim, double *a, double* tmp);


/****************************************************************************/
END_C_DECLS

#endif

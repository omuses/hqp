/*---------------------------------------------------------------------------- 
 ADOL-C -- Automatic Differentiation by Overloading in C++
 File:     drivers/taylor.h
 Revision: $Id: taylor.h,v 1.1 2004/10/13 14:18:12 e_arnold Exp $
 Contents: Easy to use drivers for the evaluation of higher order derivative
           tensors and inverse/impicit function differentiation
 
 Copyright (c) 2004
               Technical University Dresden
               Department of Mathematics
               Institute of Scientific Computing
  
 This file is part of ADOL-C. This software is provided under the terms of
 the Common Public License. Any use, reproduction, or distribution of the
 software constitutes recipient's acceptance of the terms of this license.
 See the accompanying copy of the Common Public License for more details.

 History:
          20040416 kowarz:       adapted to configure - make - make install
          19981130 olvo:         last check (includes ...)
          19981120 olvo/walther: return values
          19980914 olvo:         Jac_solv --> jac_solv
          19980806 walther:      cleanup for tensors
          19980804 walther:      new access to tensors
 
----------------------------------------------------------------------------*/
#if !defined(ADOLC_DRIVERS_TAYLOR_H)
#define ADOLC_DRIVERS_TAYLOR_H 1

#include "../common.h"

BEGIN_C_DECLS

/****************************************************************************/
/*                                                       TENSOR EVALUATIONS */

/*--------------------------------------------------------------------------*/
/* tensor_eval(tag,m,n,d,p,x[n],tensor[m][dim],S[n][p])
      with dim = ((p+d) over d) */
int tensor_eval(int TAG, int m, int n, int d, int p, double *x, 
                 double **tensor, double **S);

/*--------------------------------------------------------------------------*/
/* inverse_tensor_eval(tag,n,d,p,x,tensor[n][dim],S[n][p])
      with dim = ((p+d) over d) */
int inverse_tensor_eval(int tag, int n, int d, int p, double *x, 
                         double **tensor, double **S);

/*--------------------------------------------------------------------------*/
/*  inverse_Taylor_prop(tag,n,d,Y[n][d+1],X[n][d+1]) */
int inverse_Taylor_prop(unsigned short tag, int n, int d, 
                         double** Y, double** X);

/****************************************************************************/
/*                                                  ACCESS TO TENSOR VALUES */

/*--------------------------------------------------------------------------*/
/* tensor_value(d,m,y[m],tensori[m][dim],multi[d]) 
      with dim = ((p+d) over d) */
void tensor_value(int d, int m, double *y, double **tensor, int *multi);

/*--------------------------------------------------------------------------*/
/* void** tensorsetup(m,p,d,tensorig) */
void** tensorsetup(int m, int p, int d, double** tensorig);

/*--------------------------------------------------------------------------*/
/* void freetensor(m,p,d,tensor) */
void freetensor(int m, int p, int d, double** tensor);

/*--------------------------------------------------------------------------*/
/* int address(d, im[d]) */
int address(int d, int* im);

/****************************************************************************/
/*                                                                    UTILS */

/*--------------------------------------------------------------------------*/
/* int binomi(a,b)  ---> binomial coefficient to compute tensor dimension */
long binomi(int a, int b);

/*--------------------------------------------------------------------------*/
/* jac_solv(tag,n,x,b,sparse,mode) */
int jac_solv(unsigned short tag, int n, double* x, double* b, 
             unsigned short sparse, unsigned short mode);

END_C_DECLS

/****************************************************************************/
#endif

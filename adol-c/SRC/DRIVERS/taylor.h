#ifndef _TAYLOR_H_
#define _TAYLOR_H_ 
/*
   ------------------------------------------------------------- 
   File taylor.h of ADOL-C version 1.8.0         as of Nov/30/98
   ------------------------------------------------------------- 
   "Easy to use" drivers for evaluation of higher order 
   derivative tensors and inverse/implicit function 
   differentiation
    
   Last changes:
      981130  olvo:     last check (includes ...)
      981120  olvo/walther:  return values
      980914  olvo:     Jac_solv --> jac_solv
      980806  walther:  cleanup for tensors
      980804  walther:  new access to tensors 

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
int binomi(int a, int b);

/*--------------------------------------------------------------------------*/
/* jac_solv(tag,n,x,b,sparse,mode) */
int jac_solv(unsigned short tag, int n, double* x, double* b, 
             unsigned short sparse, unsigned short mode);


/****************************************************************************/
/*                                                               THAT'S ALL */
#ifdef __cplusplus
}
#endif

#endif













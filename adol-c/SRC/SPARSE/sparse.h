#ifndef _SPARSE_H_
#define _SPARSE_H_
/*
   --------------------------------------------------------------
   File sparse.h of ADOL-C version 1.8.2          as of Mar/09/99
   --------------------------------------------------------------
   All "Easy To Use" interfaces of SPARSE package

   Last changes:
        990308 christo: bit patterns : unsigned int -> unsigned long int
        990308 christo: mode : short -> char
        990302 christo: new interface of jac_pat(...)
        981203 olvo: untransposing reverse
        981130 olvo: newly created from adutils.h & adutilsc.h

   --------------------------------------------------------------
*/


/****************************************************************************/
/****************************************************************************/
/****************************************************************************/
/*                                                        PUBLIC INTERFACES */

#ifdef __cplusplus
/****************************************************************************/
/****************************************************************************/
/*                                                       Now the C++ THINGS */

/****************************************************************************/
/*                                           FORWARD MODE, overloaded calls */

/* nBV = number of Boolean Vectors to be packed 
                    (see Chapter Dependence Analysis, ADOL-C Documentation)
   bits_per_long = 8*sizeof(unsigned long int)
   p = nBV / bits_per_long + ( (nBV % bits_per_long) != 0 )

   For the full Jacobian matrix set 
     p = indep / bits_per_long + ((indep % bits_per_long) != 0)
   and pass a bit pattern version of the identity matrix as an argument    */

/*--------------------------------------------------------------------------*/
/*  Bit pattern propagation call, d = 1, tight version                      */
/*                                                                          */
/* forward(tag, m, n, p, x[n], X[n][p], y[m], Y[m][p], mode) : intfov       */

int forward(short, int, int, int, double*, unsigned long int**, 
                                  double*, unsigned long int**, char =0);

/*--------------------------------------------------------------------------*/
/*  Bit pattern propagation call, d = 1, safe version (no x[] and y[])      */
/*                                                                          */
/* forward(tag, m, n, p, X[n][p], Y[m][p], mode) : intfov                   */

int forward(short, int, int, int, unsigned long int**, unsigned long int**, char =0);


/****************************************************************************/
/*                                           REVERSE MODE, overloaded calls */

/* nBV = number of Boolean Vectors to be packed
                     (see Chapter Dependence Analysis, ADOL-C Documentation)
   bits_per_long = 8*sizeof(unsigned long int)
   q = nBV / bits_per_long + ( (nBV % bits_per_long) != 0 )

   For the full Jacobian matrix set 
     q = depen / bits_per_long + ((depen % bits_per_long) != 0)
   and pass a bit pattern version of the identity matrix as an argument    */

/*--------------------------------------------------------------------------*/
/*                                                                          */
/*  Bit pattern propagation call, d = 0, tight & safe version               */ 
/*                                                                          */
/* reverse(tag, m, n, q, U[q][m], Z[q][n], mode) : intfov                   */

int reverse(short, int, int, int, unsigned long int**, unsigned long int**, char =0);


/****************************************************************************/
/****************************************************************************/
/*                                                         Now the C THINGS */
extern "C" {
#endif


/*--------------------------------------------------------------------------*/
/*                                                         jacobian pattern */
/* jac_pat(tag, m, n, argument, 
           rb[1+m], cb[1+n],  
           crs[rb[0]][ crs[][0] = non-zero independent blocks per row ],
           options[2])                                                      

   if (rb == NULL) each dependent variable will be considered 
                   as a block of dependent variables.
   if (cb == NULL) each independent variable will be considered 
                   as a block of independent variables.                     */
   
               
int jac_pat(short,int,int,double*,unsigned int*,unsigned int*,
            unsigned int**,int*);


/****************************************************************************/
/*                                                               THAT'S ALL */
#ifdef __cplusplus
}
#endif

#endif

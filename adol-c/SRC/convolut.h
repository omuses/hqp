#ifndef _CONVOLUT_H_
#define _CONVOLUT_H_
/*
   ----------------------------------------------------------------
   File convolut.h of ADOL-C version 1.8.0          as of Nov/30/98
   ----------------------------------------------------------------
   Contains convolution routines used in ho_reverse.c (declaration)

   Last changed : 
        981130 olvo   last check
        980707 olvo   created this file from parts of adutilsc.h
        980616 olvo   (1) void copyAndZeroset(..)
                      (2) void inconv0(..)
                          void deconv0(..)

   ----------------------------------------------------------------
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
/*                                                               THAT'S ALL */
#ifdef __cplusplus
}
#endif

#endif

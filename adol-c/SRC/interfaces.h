#ifndef _INTERFACES_H_
#define _INTERFACES_H_
/*
   --------------------------------------------------------------
   File interfaces.h of ADOL-C version 1.8.0      as of Dec/01/98
   --------------------------------------------------------------
   Contains the declaration of the standard interfaces to 
   ADOL-C forward and reverse calls (C++, C and Fortran callable
   C functions).

   Functions prototyped here are defined in the files
            ---> uni5_for.c for
                                 zos_forward.c  (zos_forward)
                                 fos_forward.c  (fos_forward)
                                 hos_forward.c  (hos_forward)
                                 fov_forward.c  (fov_forward)
                                 hov_forward.c  (hov_forward)
                 fo_rev.c for
		                 fos_reverse.c  (fos_reverse)
		                 fov_reverse.c  (fov_reverse)
                 ho_rev.c for
                                 hos_reverse.c  (hos_reverse)
		                 hov_reverse.c  (hov_reverse)
                 interfacesC.C
                 interfacesf.c
              
  ADOL-C Abreviations :
     zos : zero-order-scalar mode   
     fos : first-order-scalar mode
     hos : higher-order-scalar mode
     fov : first-order-vector mode
     hov : higher-order-vector mode

   Last changes:
        981201 olvo: automatic include of "SPARSE/sparse.h"
        981130 olvo: newly created by unification of ADOLC-kernel
                     routines of adutils?.h
         
   History of adutils.h:
     981126 olvo     last check (p's & q's)
     980727 olvo     ec U[m][p] ---> U[p][m]   

   History of adutilsc.h:
     981126 olvo:    last check

   --------------------------------------------------------------
*/


/****************************************************************************/
/****************************************************************************/
/****************************************************************************/
/*                                                      PUBLIC EXPORTS ONLY */
#include "usrparms.h" 

#ifdef __cplusplus
/****************************************************************************/
/****************************************************************************/
/*                                                       Now the C++ THINGS */

/****************************************************************************/
/*                                           FORWARD MODE, overloaded calls */

/*--------------------------------------------------------------------------*/
/*  General scalar call. For d=0 or d=1 done by specialized code            */
/*                                                                          */
/* forward(tag, m, n, d, keep, X[n][d+1], Y[m][d+1]) : hos || fos || zos    */

int forward(short, int, int, int, int, double**, double**);

/*--------------------------------------------------------------------------*/
/*    Y can be one dimensional if m=1. d=0 or d=1 done by specialized code  */
/*                                                                          */
/* forward(tag, m, n, d, keep, X[n][d+1], Y[d+1]) : hos || fos || zos       */

int forward(short, int, int, int, int, double**, double*);

/*--------------------------------------------------------------------------*/
/*    X and Y can be one dimensional if d = 0; done by specialized code     */
/*                                                                          */
/* forward(tag, m, n, d, keep, X[n], Y[m]) : zos                            */

int forward(short, int, int, int, int, double*, double*); 

/*--------------------------------------------------------------------------*/
/*    X and Y can be one dimensional if d omitted; done by specialized code */
/*                                                                          */
/* forward(tag, m, n, keep, X[n], Y[m]) : zos                               */

int forward(short, int, int, int, double*, double*);

/*--------------------------------------------------------------------------*/
/*  General vector call                                                     */
/*                                                                          */
/* forward(tag, m, n, d, p, x[n], X[n][p][d], y[m], Y[m][p][d]) : hov       */

int forward(short, int, int, int, int, double*, double***, 
                                       double*, double***);
 
/*--------------------------------------------------------------------------*/
/*  d = 1 may be omitted. General vector call, done by specialized code     */
/*                                                                          */
/* forward(tag, m, n, p, x[n], X[n][p], y[m], Y[m][p]) : fov                */

int forward(short, int, int, int, double*, double**, double*, double**);


/****************************************************************************/
/*                                           REVERSE MODE, overloaded calls */

/*--------------------------------------------------------------------------*/
/*  General call                                                            */
/*                                                                          */
/* reverse(tag, m, n, d, u[m], Z[n][d+1]) : hos                             */

int reverse(short, int, int, int, double*, double**);

/*--------------------------------------------------------------------------*/
/*    u can be a scalar if m=1                                              */
/*                                                                          */
/* reverse(tag, m, n, d, u, Z[n][d+1]) : hos                                */

int reverse(short, int, int, int, double, double**);

/*--------------------------------------------------------------------------*/
/*    Z can be vector if d = 0; done by specialized code                    */
/*                                                                          */
/* reverse(tag, m, n, d, u[m], Z[n]) : fos                                  */

int reverse(short, int, int, int, double*, double*);

/*--------------------------------------------------------------------------*/
/*    u can be a scalar if m=1 and d=0; done by specialized code            */
/*                                                                          */
/* reverse(tag, m, n, d, u, Z[n]) : fos                                     */

int reverse(short, int, int, int, double, double*);

/*--------------------------------------------------------------------------*/
/*  General vector call                                                     */
/*                                                                          */
/* reverse(tag, m, n, d, q, U[q][m], Z[q][n][d+1], nz[q][n]) : hov          */

int reverse(short, int, int, int, int, double**, double***, short** =0);

/*--------------------------------------------------------------------------*/
/*    U can be a vector if m=1                                              */
/*                                                                          */
/* reverse(tag, m, n, d, q, U[q], Z[q][n][d+1], nz[q][n]) : hov             */

int reverse(short, int, int, int, int, double*, double***, short** = 0);

/*--------------------------------------------------------------------------*/
/*                                                                          */
/*    If d=0 then Z may be a matrix, no nz; done by specialized cod         */
/*                                                                          */
/* reverse(tag, m, n, d, q, U[q][m], Z[q][n]) : fov                         */

int reverse(short, int, int, int, int, double**, double**);

/*--------------------------------------------------------------------------*/
/*                                                                          */
/*    d=0 may be omitted, Z is a matrix, no nz; done by specialized code    */
/*                                                                          */
/* reverse(tag, m, n, q, U[q][m], Z[q][n]) : fov                            */

int reverse(short, int, int, int, double**, double**);

/*--------------------------------------------------------------------------*/
/*                                                                          */
/*    If m=1 and d=0 then U can be vector and Z a matrix but no nz.         */
/*    Done by specialized code                                              */
/*                                                                          */
/* reverse(tag, m, n, d, q, U[q], Z[q][n]) : fov                            */

int reverse(short, int, int, int, int, double*, double**);

/*--------------------------------------------------------------------------*/
/*                                                                          */
/*    If q and U are omitted they default to m and I so that as above       */
/*                                                                          */
/* reverse(tag, m, n, d, Z[q][n][d+1], nz[q][n]) : hov                      */

int reverse(short, int, int, int, double***, short** =0);


/****************************************************************************/
/****************************************************************************/
/*                                                         Now the C THINGS */
extern "C" {
#endif

/****************************************************************************/
/*                                                             FORWARD MODE */

/*--------------------------------------------------------------------------*/
/*                                                                      ZOS */
/* zos_forward(tag, m, n, keep, x[n], y[m])                                 */
/* (defined in uni5_for.c)                                                  */

int zos_forward(short,int,int,int,double*,double*);

/* zos_forward_nk(tag, m, n, x[n], y[m])                                    */ 
/* (no keep, defined in uni5_for.c, but not supported in ADOL-C 1.8)        */

int zos_forward_nk(short,int,int,double*,double*);


/*--------------------------------------------------------------------------*/
/*                                                                      FOS */
/* fos_forward(tag, m, n, keep, x[n], X[n], y[m], Y[m])                     */
/* (defined in uni5_for.c)                                                  */

int fos_forward(short,int,int,int,double*,double*,double*,double*);

/* fos_forward_nk(tag,m,n,x[n],X[n],y[m],Y[m])                              */
/* (no keep, defined in uni5_for.c, but not supported in ADOL-C 1.8)        */

int fos_forward_nk(short,int,int,double*,double*,double*,double*);


/*--------------------------------------------------------------------------*/
/*                                                                      HOS */
/* hos_forward(tag, m, n, d, keep, x[n], X[n][d], y[m], Y[m][d])            */ 
/* (defined in uni5_for.c)                                                  */

int hos_forward(short,int,int,int,int,double*,double**,double*,double**);

/* hos_forward_nk(tag, m, n, d, x[n], X[n][d], y[m], Y[m][d])               */ 
/* (no keep, defined in uni5_for.c, but not supported in ADOL-C 1.8)        */

int hos_forward_nk(short,int,int,int,double*,double**,double*,double**);

/* now pack the arrays into vectors for Fortran calling                     */
fint hos_forward_(fint*,fint*,fint*,fint*,fint*,fdouble*,fdouble*,
                                                fdouble*,fdouble*);


/*--------------------------------------------------------------------------*/
/*                                                                      FOV */
/* fov_forward(tag, m, n, p, x[n], X[n][p], y[m], Y[m][p])                  */
/* (defined in uni5_for.c)                                                  */
 
int fov_forward(short, int,int,int,double*,double**,double*,double**);
 
/* now pack the arrays into vectors for Fortran calling                     */
fint fov_forward_(fint*,fint*,fint*,fint*,fdouble*,fdouble*,
                                          fdouble*,fdouble*);


/*--------------------------------------------------------------------------*/
/*                                                                      HOV */
/* hov_forward(tag, m, n, d, p, x[n], X[n][p][d], y[m], Y[m][p][d])         */ 
/* (defined in uni5_for.c)                                                  */

int hov_forward(short,int,int,int,int,double*,double***,double*,double***);

/* now pack the arrays into vectors for Fortran calling                     */
fint hov_forward_(fint*,fint*,fint*,fint*,fint*,fdouble*,fdouble*,
                                                fdouble*,fdouble*);



/****************************************************************************/
/*                                                             REVERSE MODE */

/*--------------------------------------------------------------------------*/
/*                                                                      FOS */
/* fos_reverse(tag, m, n, u[m], z[n])                                       */
/* (defined  in fo_rev.c)                                                   */

int fos_reverse(short,int,int,double*,double*);

/* now pack the arrays into vectors for Fortran calling                     */
fint fos_reverse_(fint*,fint*,fint*,fdouble*,fdouble*);


/*--------------------------------------------------------------------------*/
/*                                                                      HOS */
/*  hos_reverse(tag, m, n, d, u[m], Z[n][d+1])                              */
/* (defined  in ho_rev.c)                                                   */

int hos_reverse(short,int,int,int,double*,double**);

/* now pack the arrays into vectors for Fortran calling                     */
fint hos_reverse_(fint*,fint*,fint*,fint*,fdouble*,fdouble*);


/*--------------------------------------------------------------------------*/
/*                                                                      FOV */
/* fov_reverse(tag, m, n, d, p, U[p][m], Z[p][n])                           */
/* (defined  in fo_rev.c)                                                   */

int fov_reverse(short,int,int,int,double**,double**);

/* now pack the arrays into vectors for Fortran calling                     */
fint fov_reverse_(fint*,fint*,fint*,fint*,fdouble*,fdouble*);


/*--------------------------------------------------------------------------*/
/*                                                                      HOV */
/* hov_reverse(tag, m, n, d, p, U[p][m], Z[p][n][d+1], nz[p][n])            */ 
/* (defined  in ho_rev.c)                                                   */

int hov_reverse(short,int,int,int,int,double**,double***,short**);

/* now pack the arrays into vectors for Fortran calling      */
fint hov_reverse_(fint*,fint*,fint*,fint*,fint*,fdouble*,fdouble*);


/****************************************************************************/
/*                                                               THAT'S ALL */
#ifdef __cplusplus
}
#endif

/*--------------------------------------------------------------------------*/
/*                                                        Automatic include */
#include "SPARSE/sparse.h"

#endif

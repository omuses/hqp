/*---------------------------------------------------------------------------- 
 ADOL-C -- Automatic Differentiation by Overloading in C++
 File:     interfaces.h
 Revision: $Id: interfaces.h,v 1.1 2004/10/13 14:18:12 e_arnold Exp $
 Contents: Declaration of the standard interfaces to ADOL-C forward and
           reverse calls (C++, C and Fortran callable C functions).

           Functions prototyped here are defined in the files
                 uni5_for.mc for
                                 zos_forward.c
                                 fos_forward.c
                                 hos_forward.c
                                 fov_forward.c
                                 hov_forward.c
                                 hov_wk_forward.c
                 fo_rev.mc for
		                 fos_reverse.c
		                 fov_reverse.c
                 ho_rev.mc for
                                 hos_reverse.c
                                 hos_ov_reverse.c
		                 hov_reverse.c
                                 hos_ti_reverse.c
		                 hov_ti_reverse.c
                 interfacesC.C
                 interfacesf.c
              
           ADOL-C Abreviations :
                 zos : zero-order-scalar mode   
                 fos : first-order-scalar mode
                 hos : higher-order-scalar mode
                 fov : first-order-vector mode
                 hov : higher-order-vector mode
                 wk  : with keep
                 ov  : over vector (forward)
                 ti  : Taylor input

 Copyright (c) 2004
               Technical University Dresden
               Department of Mathematics
               Institute of Scientific Computing

 This file is part of ADOL-C. This software is provided under the terms of
 the Common Public License. Any use, reproduction, or distribution of the
 software constitutes recipient's acceptance of the terms of this license.
 See the accompanying copy of the Common Public License for more details.

 History:
          20040414 kowarz: adapted to configure - make - make install
          20031218 andrea: hos_forward_partx
          20020730 olvo:   allowing input of higher order adjoints
                           (hos_ti_reverse & hov_ti_reverse)
          20010719 andrea: forward vector modes with keep
                           (hov_wk_forward & hos_ov_reverse)  
----------------------------------------------------------------------------*/
#if !defined(ADOLC_INTERFACES_H)
#define ADOLC_INTERFACES_H 1

#include "common.h"
#include "sparse/sparse.h"

/****************************************************************************/
/****************************************************************************/
/*                                                       Now the C++ THINGS */
#if defined(__cplusplus)

/****************************************************************************/
/*                                           FORWARD MODE, overloaded calls */

/*--------------------------------------------------------------------------*/
/*    General scalar call. For d=0 or d=1 done by specialized code          */
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
/*    If d=0 then Z may be a matrix, no nz; done by specialized code        */
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

#endif

/****************************************************************************/
/****************************************************************************/
/*                                                         Now the C THINGS */
BEGIN_C_DECLS

/****************************************************************************/
/*                                                             FORWARD MODE */

/*--------------------------------------------------------------------------*/
/*                                                                      ZOS */
/* zos_forward(tag, m, n, keep, x[n], y[m])                                 */
/* (defined in uni5_for.mc)                                                 */
int zos_forward(short,int,int,int,double*,double*);

/* zos_forward_nk(tag, m, n, x[n], y[m])                                    */ 
/* (no keep, defined in uni5_for.c, but not supported in ADOL-C 1.8)        */
int zos_forward_nk(short,int,int,double*,double*);

/* zos_forward_partx(tag, m, n, ndim[n], x[n][d], y[m])                     */ 
/* (based on zos_forward)                                                   */

int zos_forward_partx(short,int,int,int*,double**,double*);

/*--------------------------------------------------------------------------*/
/*                                                                      FOS */
/* fos_forward(tag, m, n, keep, x[n], X[n], y[m], Y[m])                     */
/* (defined in uni5_for.mc)                                                 */
int fos_forward(short,int,int,int,double*,double*,double*,double*);

/* fos_forward_nk(tag,m,n,x[n],X[n],y[m],Y[m])                              */
/* (no keep, defined in uni5_for.c, but not supported in ADOL-C 1.8)        */
int fos_forward_nk(short,int,int,double*,double*,double*,double*);

/* fos_forward_partx(tag, m, n, ndim[n], x[n][][2], y[m][2])                */ 
/* (based on fos_forward)                                                   */
int fos_forward_partx(short,int,int,int*,double***,double**);

/*--------------------------------------------------------------------------*/
/*                                                                      HOS */
/* hos_forward(tag, m, n, d, keep, x[n], X[n][d], y[m], Y[m][d])            */ 
/* (defined in uni5_for.mc)                                                 */
int hos_forward(short,int,int,int,int,double*,double**,double*,double**);

/* hos_forward_nk(tag, m, n, d, x[n], X[n][d], y[m], Y[m][d])               */ 
/* (no keep, defined in uni5_for.c, but not supported in ADOL-C 1.8)        */
int hos_forward_nk(short,int,int,int,double*,double**,double*,double**);

/* hos_forward_partx(tag, m, n, ndim[n], d, X[n][d+1], Y[m][d+1])           */ 
/* (defined in forward_partx.c)                                             */
int hos_forward_partx(short,int,int,int*,int,double***,double**);

/* now pack the arrays into vectors for Fortran calling                     */
fint hos_forward_(fint*,fint*,fint*,fint*,fint*,fdouble*,fdouble*,
                                                fdouble*,fdouble*);

/*--------------------------------------------------------------------------*/
/*                                                                      FOV */
/* fov_forward(tag, m, n, p, x[n], X[n][p], y[m], Y[m][p])                  */
/* (defined in uni5_for.mc)                                                 */
int fov_forward(short, int,int,int,double*,double**,double*,double**);
 
/* now pack the arrays into vectors for Fortran calling                     */
fint fov_forward_(fint*,fint*,fint*,fint*,fdouble*,fdouble*,
                                          fdouble*,fdouble*);

/*  fov_forward_partx(tag, m, n, ndim[n], p,                                */
/*                    x[n][], X[n][][p],y[m], Y[m][p])                      */ 
int fov_forward_partx(short, int, int, int*, int,
                      double**, double***, double*, double**);

/*--------------------------------------------------------------------------*/
/*                                                                      HOV */
/* hov_forward(tag, m, n, d, p, x[n], X[n][p][d], y[m], Y[m][p][d])         */
/* (defined in uni5_for.mc)                                                 */
int hov_forward(short,int,int,int,int,double*,double***,double*,double***);

/* now pack the arrays into vectors for Fortran calling                     */
fint hov_forward_(fint*,fint*,fint*,fint*,fint*,fdouble*,fdouble*,
                                                fdouble*,fdouble*);

/*  hov_forward_partx(tag, m, n, ndim[n], d, p,                             */
/*                    x[n][], X[n][][p][d], y[m], Y[m][p][d])               */
int hov_forward_partx(short, int, int, int*, int, int,
                      double**, double****, double*, double***);

/*--------------------------------------------------------------------------*/
/*                                                                   HOV_WK */
/* hov_wk_forward(tag, m, n, d, p, x[n], X[n][p][d], y[m], Y[m][p][d])      */ 
/* (defined in uni5_for.mc)                                                 */
int hov_wk_forward(short,int,int,int,int,int,double*,double***,double*,
                                                               double***);

/* now pack the arrays into vectors for Fortran calling                     */
fint hov_wk_forward_(fint*,fint*,fint*,fint*,fint*,fint*,fdouble*,fdouble*,
                                                fdouble*,fdouble*);

/****************************************************************************/
/*                                                             REVERSE MODE */

/*--------------------------------------------------------------------------*/
/*                                                                      FOS */
/* fos_reverse(tag, m, n, u[m], z[n])                                       */
/* (defined  in fo_rev.mc)                                                  */
int fos_reverse(short,int,int,double*,double*);

/* now pack the arrays into vectors for Fortran calling                     */
fint fos_reverse_(fint*,fint*,fint*,fdouble*,fdouble*);

/*--------------------------------------------------------------------------*/
/*                                                                      HOS */
/*  hos_reverse(tag, m, n, d, u[m], Z[n][d+1])                              */
/* (defined  in ho_rev.mc)                                                  */
int hos_reverse(short,int,int,int,double*,double**);

/* now pack the arrays into vectors for Fortran calling                     */
fint hos_reverse_(fint*,fint*,fint*,fint*,fdouble*,fdouble*);

/*--------------------------------------------------------------------------*/
/*                                                                   HOS_TI */
/*  hos_ti_reverse(tag, m, n, d, U[m][d+1], Z[n][d+1])                      */
/* (defined  in ho_rev.mc)                                                  */
int hos_ti_reverse(short,int,int,int,double**,double**);

/* now pack the arrays into vectors for Fortran calling                     */
fint hos_ti_reverse_(fint*,fint*,fint*,fint*,fdouble*,fdouble*);

/*--------------------------------------------------------------------------*/
/*                                                                   HOS_OV */
/*  hos_ov_reverse(tag, m, n, d, u[m], Z[n][d+1])                           */
/* (defined  in ho_rev.mc)                                                  */
int hos_ov_reverse(short,int,int,int,int,double**,double***);

/* now pack the arrays into vectors for Fortran calling                     */
fint hos_ov_reverse_(fint*,fint*,fint*,fint*,fint*,fdouble*,fdouble***);

/*--------------------------------------------------------------------------*/
/*                                                                      FOV */
/* fov_reverse(tag, m, n, d, p, U[p][m], Z[p][n])                           */
/* (defined  in fo_rev.mc)                                                  */
int fov_reverse(short,int,int,int,double**,double**);

/* now pack the arrays into vectors for Fortran calling                     */
fint fov_reverse_(fint*,fint*,fint*,fint*,fdouble*,fdouble*);

/*--------------------------------------------------------------------------*/
/*                                                                      HOV */
/* hov_reverse(tag, m, n, d, p, U[p][m], Z[p][n][d+1], nz[p][n])            */ 
/* (defined  in ho_rev.mc)                                                  */
int hov_reverse(short,int,int,int,int,double**,double***,short**);

/* now pack the arrays into vectors for Fortran calling      */
fint hov_reverse_(fint*,fint*,fint*,fint*,fint*,fdouble*,fdouble*);

/*--------------------------------------------------------------------------*/
/*                                                                   HOV_TI */
/* hov_ti_reverse(tag, m, n, d, p, U[p][m][d+1], Z[p][n][d+1], nz[p][n])    */
/* (defined  in ho_rev.mc)                                                  */
int hov_ti_reverse(short,int,int,int,int,double***,double***,short**);

/* now pack the arrays into vectors for Fortran calling      */
fint hov_ti_reverse_(fint*,fint*,fint*,fint*,fint*,fdouble*,fdouble*);

END_C_DECLS

#endif

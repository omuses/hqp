
/**************************************************************************
**
** Copyright (C) 1993 David E. Steward & Zbigniew Leyk, all rights reserved.
**
**			     Meschach Library
** 
** This Meschach Library is provided "as is" without any express 
** or implied warranty of any kind with respect to this software. 
** In particular the authors shall not be liable for any direct, 
** indirect, special, incidental or consequential damages arising 
** in any way from use of the software.
** 
** Everyone is granted permission to copy, modify and redistribute this
** Meschach Library, provided:
**  1.  All copies contain this copyright notice.
**  2.  All modified copies shall carry a notice stating who
**      made the last modification and the date of such modification.
**  3.  No charge is made for this software or works derived from it.  
**      This clause shall not be construed as constraining other software
**      distributed on the same medium as this software, nor is a
**      distribution fee considered a charge.
**
***************************************************************************/


/*
	2nd header file for Meschach's complex routines.
	This file contains declarations for complex factorisation/solve
	routines.

*/


#ifndef ZMATRIX2H
#define ZMATRIX2H

#include "zmatrix.h"

MESCH__BEGIN_DECLS

MESCH_API ZVEC	*zUsolve(ZMAT *matrix, ZVEC *b, ZVEC *out, double diag);
MESCH_API ZVEC	*zLsolve(ZMAT *matrix, ZVEC *b, ZVEC *out, double diag);
MESCH_API ZVEC	*zUAsolve(ZMAT *U, ZVEC *b, ZVEC *out, double diag);
MESCH_API ZVEC	*zDsolve(ZMAT *A, ZVEC *b, ZVEC *x);
MESCH_API ZVEC	*zLAsolve(ZMAT *L, ZVEC *b, ZVEC *out, double diag);

MESCH_API ZVEC	*zhhvec(ZVEC *,int,Real *,ZVEC *,complex *);
MESCH_API ZVEC	*zhhtrvec(ZVEC *,double,int,ZVEC *,ZVEC *);
MESCH_API ZMAT	*zhhtrrows(ZMAT *,int,int,ZVEC *,double);
MESCH_API ZMAT	*zhhtrcols(ZMAT *,int,int,ZVEC *,double);
MESCH_API ZMAT     *zHfactor(ZMAT *,ZVEC *);
MESCH_API ZMAT     *zHQunpack(ZMAT *,ZVEC *,ZMAT *,ZMAT *);

MESCH_API ZMAT	*zQRfactor(ZMAT *A, ZVEC *diag);
MESCH_API ZMAT	*zQRCPfactor(ZMAT *A, ZVEC *diag, PERM *px);
MESCH_API ZVEC	*_zQsolve(ZMAT *QR, ZVEC *diag, ZVEC *b, ZVEC *x, ZVEC *tmp);
MESCH_API ZMAT	*zmakeQ(ZMAT *QR, ZVEC *diag, ZMAT *Qout);
MESCH_API ZMAT	*zmakeR(ZMAT *QR, ZMAT *Rout);
MESCH_API ZVEC	*zQRsolve(ZMAT *QR, ZVEC *diag, ZVEC *b, ZVEC *x);
MESCH_API ZVEC	*zQRAsolve(ZMAT *QR, ZVEC *diag, ZVEC *b, ZVEC *x);
MESCH_API ZVEC	*zQRCPsolve(ZMAT *QR,ZVEC *diag,PERM *pivot,ZVEC *b,ZVEC *x);
MESCH_API ZVEC	*zUmlt(ZMAT *U, ZVEC *x, ZVEC *out);
MESCH_API ZVEC	*zUAmlt(ZMAT *U, ZVEC *x, ZVEC *out);
MESCH_API double	zQRcondest(ZMAT *QR);

MESCH_API ZVEC	*zLsolve(ZMAT *, ZVEC *, ZVEC *, double);
MESCH_API ZMAT	*zset_col(ZMAT *, int, ZVEC *);

MESCH_API ZMAT	*zLUfactor(ZMAT *A, PERM *pivot);
MESCH_API ZVEC	*zLUsolve(ZMAT *A, PERM *pivot, ZVEC *b, ZVEC *x);
MESCH_API ZVEC	*zLUAsolve(ZMAT *LU, PERM *pivot, ZVEC *b, ZVEC *x);
MESCH_API ZMAT	*zm_inverse(ZMAT *A, ZMAT *out);
MESCH_API double	zLUcondest(ZMAT *LU, PERM *pivot);

MESCH_API void	zgivens(complex, complex, Real *, complex *);
MESCH_API ZMAT	*zrot_rows(ZMAT *A, int i, int k, double c, complex s,
			   ZMAT *out);
MESCH_API ZMAT	*zrot_cols(ZMAT *A, int i, int k, double c, complex s,
			   ZMAT *out);
MESCH_API ZVEC	*rot_zvec(ZVEC *x, int i, int k, double c, complex s,
			  ZVEC *out);
MESCH_API ZMAT	*zschur(ZMAT *A,ZMAT *Q);
/* MESCH_API ZMAT	*schur_vecs(ZMAT *T,ZMAT *Q,X_re,X_im) */

MESCH__END_DECLS

#endif


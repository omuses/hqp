
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


/* Sparse matrix factorise/solve header */
/* RCS id: $Id: sparse2.h,v 1.3 2002/12/09 10:57:47 e_arnold Exp $ */



#ifndef SPARSE2H

#define SPARSE2H

#include "sparse.h"

MESCH__BEGIN_DECLS

MESCH_API SPMAT	*spCHfactor(SPMAT *), *spICHfactor(SPMAT *), 
    *spCHsymb(SPMAT *);
MESCH_API VEC	*spCHsolve(SPMAT *,VEC *,VEC *);

MESCH_API SPMAT	*spLUfactor(SPMAT *,PERM *,double);
MESCH_API SPMAT	*spILUfactor(SPMAT *,double);
MESCH_API VEC	*spLUsolve(SPMAT *,PERM *,VEC *,VEC *),
    *spLUTsolve(SPMAT *,PERM *,VEC *,VEC *);

MESCH_API SPMAT	*spBKPfactor(SPMAT *, PERM *, PERM *, double);
MESCH_API VEC	*spBKPsolve(SPMAT *, PERM *, PERM *, VEC *, VEC *);

MESCH_API VEC	*pccg(VEC *(*A)(),void *A_par,VEC *(*M_inv)(),void *M_par,
		      VEC *b, double tol,VEC *x);
MESCH_API VEC	*sp_pccg(SPMAT *,SPMAT *,VEC *,double,VEC *);
MESCH_API VEC	*cgs(VEC *(*A)(),void *A_par,VEC *b,VEC *r0,double tol,VEC *x);
MESCH_API VEC	*sp_cgs(SPMAT *,VEC *,VEC *,double,VEC *);
MESCH_API VEC	*lsqr(VEC *(*A)(),VEC *(*AT)(),void *A_par,VEC *b,
		      double tol,VEC *x);
MESCH_API VEC	*sp_lsqr(SPMAT *,VEC *,double,VEC *);
MESCH_API int	cg_set_maxiter(int);

MESCH_API void	lanczos(VEC *(*A)(),void *A_par,int m,VEC *x0,VEC *a,VEC *b,
			Real *beta_m1,MAT *Q);
MESCH_API void	sp_lanczos(SPMAT *,int,VEC *,VEC *,VEC *,Real *,MAT *);
MESCH_API VEC	*lanczos2(VEC *(*A)(),void *A_par,int m,VEC *x0,VEC *evals,
			  VEC *err_est);
MESCH_API VEC	*sp_lanczos2(SPMAT *,int,VEC *,VEC *,VEC *);
MESCH_API void    scan_to(SPMAT *,IVEC *,IVEC *,IVEC *,int);
MESCH_API row_elt  *chase_col(SPMAT *,int,int *,int *,int);
MESCH_API row_elt  *chase_past(SPMAT *,int,int *,int *,int);
MESCH_API row_elt  *bump_col(SPMAT *,int,int *,int *);

MESCH__END_DECLS

#endif

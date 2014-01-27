
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
	Header file for ``matrix2.a'' library file
*/


#ifndef MATRIX2H
#define MATRIX2H

#include "matrix.h"

MESCH__BEGIN_DECLS

/* Unless otherwise specified, factorisation routines overwrite the
   matrix that is being factorised */

                 /* forms Bunch-Kaufman-Parlett factorisation for
                        symmetric indefinite matrices */
MESCH_API MAT	*BKPfactor(MAT *A,PERM *pivot,PERM *blocks),
                 /* Cholesky factorisation of A
                        (symmetric, positive definite) */
		*CHfactor(MAT *A),
                /* LU factorisation of A (with partial pivoting) */ 
                *LUfactor(MAT *A,PERM *pivot),
                /* QR factorisation of A; need dim(diag) >= # rows of A */
		*QRfactor(MAT *A,VEC *diag),
                /* QR factorisation of A with column pivoting */
		*QRCPfactor(MAT *A,VEC *diag,PERM *pivot),
                /* L.D.L^T factorisation of A */
		*LDLfactor(MAT *A), 
                /* Hessenberg factorisation of A -- for schur() */
                *Hfactor(MAT *A,VEC *diag1,VEC *diag2),
                /* modified Cholesky factorisation of A;
                        actually factors A+D, D diagonal with no
                        diagonal entry in the factor < sqrt(tol) */
                *MCHfactor(MAT *A,Real tol),
		*m_inverse(MAT *A,MAT *out);

                /* returns condition estimate for A after LUfactor() */
MESCH_API Real	LUcondest(MAT *A,PERM *pivot),
                /* returns condition estimate for Q after QRfactor() */
                QRcondest(MAT *A);

/* Note: The make..() and ..update() routines assume that the factorisation
        has already been carried out */

     /* Qout is the "Q" (orthongonal) matrix from QR factorisation */
MESCH_API MAT	*makeQ(MAT *A,VEC *diag,MAT *Qout),
                /* Rout is the "R" (upper triangular) matrix
                        from QR factorisation */
		*makeR(MAT *A,MAT *Rout),
                /* Qout is orthogonal matrix in Hessenberg factorisation */
		*makeHQ(MAT *A,VEC *diag1,VEC *diag2,MAT *Qout),
                /* Hout is the Hessenberg matrix in Hessenberg factorisation */
		*makeH(MAT *A,MAT *Hout);

                /* updates L.D.L^T factorisation for A <- A + alpha.u.u^T */
MESCH_API MAT	*LDLupdate(MAT *A,VEC *u,Real alpha),
                /* updates QR factorisation for QR <- Q.(R+u.v^T)
		   Note: we need explicit Q & R matrices,
                        from makeQ() and makeR() */
		*QRupdate(MAT *Q,MAT *R,VEC *u,VEC *v);

/* Solve routines assume that the corresponding factorisation routine
        has already been applied to the matrix along with auxiliary
        objects (such as pivot permutations)

        These solve the system A.x = b,
        except for LUTsolve and QRTsolve which solve the transposed system
                                A^T.x. = b.
        If x is NULL on entry, then it is created.
*/

MESCH_API VEC	*BKPsolve(MAT *A,PERM *pivot,PERM *blocks,VEC *b,VEC *x),
		*CHsolve(MAT *A,VEC *b,VEC *x),
		*LDLsolve(MAT *A,VEC *b,VEC *x),
		*LUsolve(MAT *A,PERM *pivot,VEC *b,VEC *x),
		*_Qsolve(MAT *A,VEC *,VEC *,VEC *, VEC *),
		*QRsolve(MAT *A,VEC *,VEC *b,VEC *x),
    		*QRTsolve(MAT *A,VEC *,VEC *b,VEC *x),


     /* Triangular equations solve routines;
        U for upper triangular, L for lower traingular, D for diagonal
        if diag_val == 0.0 use that values in the matrix */

		*Usolve(MAT *A,VEC *b,VEC *x,Real diag_val),
		*Lsolve(MAT *A,VEC *b,VEC *x,Real diag_val),
		*Dsolve(MAT *A,VEC *b,VEC *x),
		*LTsolve(MAT *A,VEC *b,VEC *x,Real diag_val),
		*UTsolve(MAT *A,VEC *b,VEC *x,Real diag_val),
                *LUTsolve(MAT *A,PERM *,VEC *,VEC *),
                *QRCPsolve(MAT *QR,VEC *diag,PERM *pivot,VEC *b,VEC *x);

MESCH_API BAND    *bdLUfactor(BAND *A,PERM *pivot),
                *bdLDLfactor(BAND *A);
MESCH_API VEC     *bdLUsolve(BAND *A,PERM *pivot,VEC *b,VEC *x),
                *bdLDLsolve(BAND *A,VEC *b,VEC *x);



MESCH_API VEC	*hhvec(VEC *,u_int,Real *,VEC *,Real *);
MESCH_API VEC	*hhtrvec(VEC *,Real,u_int,VEC *,VEC *);
MESCH_API MAT	*hhtrrows(MAT *,u_int,u_int,VEC *,Real);
MESCH_API MAT	*hhtrcols(MAT *,u_int,u_int,VEC *,Real);

MESCH_API void	givens(Real,Real,Real *,Real *);
MESCH_API VEC	*rot_vec(VEC *,u_int,u_int,Real,Real,VEC *); /* in situ */
MESCH_API MAT	*rot_rows(MAT *,u_int,u_int,Real,Real,MAT *); /* in situ */
MESCH_API MAT	*rot_cols(MAT *,u_int,u_int,Real,Real,MAT *); /* in situ */


/* eigenvalue routines */

               /* compute eigenvalues of tridiagonal matrix
                  with diagonal entries a[i], super & sub diagonal entries
                  b[i]; eigenvectors stored in Q (if not NULL) */
MESCH_API VEC	*trieig(VEC *a,VEC *b,MAT *Q),
                 /* sets out to be vector of eigenvectors; eigenvectors
                   stored in Q (if not NULL). A is unchanged */
		*symmeig(MAT *A,MAT *Q,VEC *out);

               /* computes real Schur form = Q^T.A.Q */
MESCH_API MAT	*schur(MAT *A,MAT *Q);
         /* computes real and imaginary parts of the eigenvalues
                        of A after schur() */
MESCH_API void	schur_evals(MAT *A,VEC *re_part,VEC *im_part);
          /* computes real and imaginary parts of the eigenvectors
                        of A after schur() */
MESCH_API MAT	*schur_vecs(MAT *T,MAT *Q,MAT *X_re,MAT *X_im);


/* singular value decomposition */

        /* computes singular values of bi-diagonal matrix with
                   diagonal entries a[i] and superdiagonal entries b[i];
                   singular vectors stored in U and V (if not NULL) */
MESCH_API VEC	*bisvd(VEC *a,VEC *b,MAT *U,MAT *V),
               /* sets out to be vector of singular values;
                   singular vectors stored in U and V */
	*svd(MAT *A,MAT *U,MAT *V,VEC *out);

/* matrix powers and exponent */
MESCH_API MAT  *_m_pow(MAT *,int,MAT *,MAT *);
MESCH_API MAT  *m_pow(MAT *,int, MAT *);
MESCH_API MAT  *m_exp(MAT *,Real,MAT *);
MESCH_API MAT  *_m_exp(MAT *,Real,MAT *,int *,int *);
MESCH_API MAT  *m_poly(MAT *,VEC *,MAT *);

/* FFT */
MESCH_API void fft(VEC *,VEC *);
MESCH_API void ifft(VEC *,VEC *);

MESCH__END_DECLS

#endif

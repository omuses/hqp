/*
 * addon_hqp.h
 *   - some Meschach dense matrix add-ons and extensions
 *     developed for the Hqp project
 *
 * E. Arnold  03/07/97
 *
 */

/*
    Copyright (C) 1997--2002  Eckhard Arnold

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Library General Public
    License as published by the Free Software Foundation; 
    version 2 of the License.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Library General Public License for more details.

    You should have received a copy of the GNU Library General Public
    License along with this library (file COPYING.LIB);
    if not, write to the Free Software Foundation, Inc.,
    59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef ADDON_HQP_H
#define ADDON_HQP_H

#include "matrix.h"
#include "sparse.h"

MESCH__BEGIN_DECLS

/* Entry level access to data structures                  */
/* Redefine some macros from matrix.h with typecast (int) */
/* Define access macros for integer vectors               */
#ifdef DEBUG
/* returns x[i] */
#undef  v_entry
#define	v_entry(x,i)	(((int) (i) < 0 || (int) (i) >= (int) (x)->dim) ? \
			 m_error(E_BOUNDS,"v_entry"), 0.0 : (x)->ve[i] )

/* x[i] <- val */
#undef	v_set_val
#define	v_set_val(x,i,val) ((x)->ve[i] = ((int) (i) < 0 || \
					  (int) (i) >= (int) (x)->dim) ? \
			    m_error(E_BOUNDS,"v_set_val"), 0.0 : (val))

/* x[i] <- x[i] + val */
#undef	v_add_val
#define	v_add_val(x,i,val) ((x)->ve[i] += ((int) (i) < 0 || \
					   (int) (i) >= (int) (x)->dim) ? \
			    m_error(E_BOUNDS,"v_add_val"), 0.0 : (val))

/* x[i] <- x[i] - val */
#undef  v_sub_val
#define	v_sub_val(x,i,val) ((x)->ve[i] -= ((int) (i) < 0 || \
					   (int) (i) >= (int) (x)->dim) ? \
			    m_error(E_BOUNDS,"v_sub_val"), 0.0 : (val))

/* returns A[i][j] */
#undef  m_entry
#define	m_entry(A,i,j)	(((int) (i) < 0 || (int) (i) >= (int) (A)->m || \
			  (int) (j) < 0 || (int) (j) >= (int) (A)->n) ? \
			 m_error(E_BOUNDS,"m_entry"), 0.0 : (A)->me[i][j] )

/* A[i][j] <- val */
#undef  m_set_val
#define	m_set_val(A,i,j,val) ((A)->me[i][j] = ((int) (i) < 0 || \
					       (int) (i) >= (int) (A)->m || \
					       (int) (j) < 0 || \
					       (int) (j) >= (int) (A)->n) ? \
			      m_error(E_BOUNDS,"m_set_val"), 0.0 : (val) )

/* A[i][j] <- A[i][j] + val */
#undef  m_add_val
#define	m_add_val(A,i,j,val) ((A)->me[i][j] += ((int) (i) < 0 || \
						(int) (i) >= (int) (A)->m || \
						(int) (j) < 0 || \
						(int) (j) >= (int) (A)->n) ? \
			      m_error(E_BOUNDS,"m_add_val"), 0.0 : (val) )

/* A[i][j] <- A[i][j] - val */
#undef  m_sub_val
#define	m_sub_val(A,i,j,val) ((A)->me[i][j] -= ((int) (i) < 0 || \
						(int) (i) >= (int) (A)->m || \
						(int) (j) < 0 || \
						(int) (j) >= (int) (A)->n) ? \
			      m_error(E_BOUNDS,"m_sub_val"), 0.0 : (val) )

/* returns x[i] */
#define	iv_entry(x,i)	(((int) (i) < 0 || (int) (i) >= (int) (x)->dim) ? \
			 m_error(E_BOUNDS,"iv_entry"), 0 : (x)->ive[i] )

/* x[i] <- val */
#define	iv_set_val(x,i,val) ((x)->ive[i] = ((int) (i) < 0 || \
					    (int) (i) >= (int) (x)->dim) ? \
			    m_error(E_BOUNDS,"iv_set_val"), 0 : (val))

/* x[i] <- x[i] + val */
#define	iv_add_val(x,i,val) ((x)->ive[i] += ((int) (i) < 0 || \
					     (int) (i) >= (int) (x)->dim) ? \
			    m_error(E_BOUNDS,"iv_add_val"), 0 : (val))

/* x[i] <- x[i] - val */
#define	iv_sub_val(x,i,val) ((x)->ive[i] -= ((int) (i) < 0 || \
					     (int) (i) >= (int) (x)->dim) ? \
			    m_error(E_BOUNDS,"iv_sub_val"), 0 : (val))

#else

/* returns x[i] */
#define	iv_entry(x,i)		((x)->ive[i])

/* x[i] <- val */
#define	iv_set_val(x,i,val)	((x)->ive[i]  = (val))

/* x[i] <- x[i] + val */
#define	iv_add_val(x,i,val)	((x)->ive[i] += (val))

 /* x[i] <- x[i] - val */
#define	iv_sub_val(x,i,val)	((x)->ive[i] -= (val))

#endif

/*----------   Meschach extensions   ----------*/

/*
 * m_output_g -- Prints internal parameters of MAT * data structure
 */
MESCH_API void m_output_g(const MAT *A);

/*
 * v_output_g -- Prints internal parameters of VEC * data structure
 */
MESCH_API void v_output_g(const VEC *A);

/*
 * m_copy1 -- copies matrix into new area, resizes out to correct size
 */
MESCH_API MAT *m_copy1(const MAT *in, MAT *out);

/*
 * v_copy1 -- copies vector into new area, resizes out to correct size
 */
MESCH_API VEC *v_copy1(const VEC *in, VEC *out);

/*
 * v_diag -- get diagonal entries of a matrix
 */
MESCH_API VEC *v_diag(MAT *A, VEC *C);

/*
 * dm_mlt -- diagonal matrix matrix multiplication
 */
MESCH_API MAT *dm_mlt(MAT *A, VEC *B, MAT *C);

/*
 * md_mlt -- matrix diagonal matrix multiplication
 */
MESCH_API MAT *md_mlt(MAT *A, VEC *B, MAT *C);

/*
 * m_symm -- make square matrix symmetric
 */
MESCH_API MAT *m_symm(const MAT *in, MAT *out);

/*
 * rel_symm -- check square matrix for symmetry
 */
MESCH_API double rel_symm(const MAT *in);

/*
 * CHsolveM -- CHsolve with matrix right hand side
 */
MESCH_API MAT *CHsolveM(MAT *A, MAT *B, MAT *C);

/*
 * CHsolveMT -- CHsolve with transposed matrix right hand side
 */
MESCH_API MAT *CHsolveMT(MAT *A, MAT *B, MAT *C);

/*
 * BKPsolveM -- BKPsolve with matrix right hand side
 */
MESCH_API MAT *BKPsolveM(MAT *A, PERM *pivot, PERM *blocks, MAT *B, MAT *C);

/*
 * BKPsolveMT -- BKPsolve with transposed matrix right hand side
 */
MESCH_API MAT *BKPsolveMT(MAT *A, PERM *pivot, PERM *blocks, MAT *B, MAT *C);

/*
 * LUsolveM -- LUsolve with matrix right hand side
 */
MESCH_API MAT *LUsolveM(MAT *A, PERM *pivot, MAT *B, MAT *C);

/*
 * LUsolveMT -- LUsolve with transposed matrix right hand side
 */
MESCH_API MAT *LUsolveMT(MAT *A, PERM *pivot, MAT *B, MAT *C);

/*
 * GE_QP -- Generalized elimination for quadratic programming
 */
MESCH_API void GE_QP(MAT *A, MAT *S, MAT *Z, MAT *PQ, double eps);

/*
 * m_concatc -- Matrix concatenation by rows C = [A; B]
 */
MESCH_API MAT *m_concatc(MAT *A, MAT *B, MAT *C);

/*
 * v_concat -- Vector concatenation C = [A; B]
 */
MESCH_API VEC *v_concat(VEC *A, VEC *B, VEC *C);

/*
 * v_move_iv -- C(i) = A(iv(i))
 */
MESCH_API VEC *v_move_iv(const VEC *A, const IVEC *iv, VEC *C);

/*
 * v_set_iv -- C(iv(i)) = A(i)
 */
MESCH_API VEC *v_set_iv(VEC *A, IVEC *iv, VEC *C);

/* 
 * v_dist2 -- 2-norm (Euclidean norm) of vector difference 
 */
double v_dist2(const VEC *x, const VEC *y);

/* 
 * v_dist_inf -- infinity-norm (supremum norm) of vector difference 
 */
double v_dist_inf(const VEC *x, const VEC *y);

MESCH__END_DECLS

#endif

/*
 * meschext_ea.h
 *   - some Meschach add-ons and extensions
 *
 * E. Arnold  03/07/97
 *            2001-08-16 sprow_norm1
 *
 */

/*
    Copyright (C) 1997--2001  Eckhard Arnold

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

#ifndef MESCHEXT_EA_H
#define MESCHEXT_EA_H

extern "C" {
#include <meschach/matrix.h>
#include <meschach/sparse.h>
}

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

//----------   Meschach extensions   ----------

//   Prints internal parameters of MAT * data structure.
extern void m_output_g(MAT *A);

//   Prints internal parameters of VEC * data structure.
extern void v_output_g(const VEC *A);

//   m_copy1 -- copies matrix into new area, resizes out to correct size.
extern MAT *m_copy1(const MAT *in, MAT *out);

//   v_copy1 -- copies vector into new area, resizes out to correct size.
extern VEC *v_copy1(const VEC *in, VEC *out);

//   v_diag -- get diagonal entries of a matrix
extern VEC *v_diag(MAT *A, VEC *C);

//   dm_mlt -- diagonal matrix matrix multiplication
extern MAT *dm_mlt(MAT *A, VEC *B, MAT *C);

//   md_mlt -- matrix diagonal matrix multiplication
extern MAT *md_mlt(MAT *A, VEC *B, MAT *C);

//   m_symm -- make square matrix symmetric,
extern MAT *m_symm(const MAT *in, MAT *out);

//   rel_symm -- check square matrix for symmetry.
extern double rel_symm(const MAT *in);

//   CHsolve with matrix right hand side.
extern MAT *CHsolveM(MAT *A, MAT *B, MAT *C);

//   CHsolve with transposed matrix right hand side.
extern MAT *CHsolveMT(MAT *A, MAT *B, MAT *C);

//   BKPsolve with matrix right hand side.
extern MAT *BKPsolveM(MAT *A, PERM *pivot, PERM *blocks, MAT *B, MAT *C);

//   BKPsolve with transposed matrix right hand side.
extern MAT *BKPsolveMT(MAT *A, PERM *pivot, PERM *blocks, MAT *B, MAT *C);

//   LUsolve with matrix right hand side.
extern MAT *LUsolveM(MAT *A, PERM *pivot, MAT *B, MAT *C);

//   LUsolve with transposed matrix right hand side.
extern MAT *LUsolveMT(MAT *A, PERM *pivot, MAT *B, MAT *C);

//   GE_QP -- Generalized elimination for quadratic programming
extern void GE_QP(MAT *A, MAT *S, MAT *Z, MAT *PQ, double eps);

//   Matrix concatenation by rows C = [A; B].
extern MAT *m_concatc(MAT *A, MAT *B, MAT *C);

//   Vector concatenation C = [A; B]
extern VEC *v_concat(VEC *A, VEC *B, VEC *C);

//   C(i) = A(iv(i))
extern VEC *v_move_iv(const VEC *A, const IVEC *iv, VEC *C);

//   C(iv(i)) = A(i)
extern VEC *v_set_iv(VEC *A, IVEC *iv, VEC *C);

//   Euclidean norm of a SPROW.
extern Real sprow_norm2(const SPROW *r);

//   L1 norm of a SPROW.
extern Real sprow_norm1(const SPROW *r);

//   Copy a block from sparse matrix src to dense matrix dst
extern void sp_extract_mat_iv(const SPMAT *src, const IVEC *iv, int j_offs, 
			      MAT *dst);

//   Block from sparse matrix - dense matrix multiplication
extern MAT *bspm_mlt(const SPMAT *A, int i0, int j0, int m0, 
		     const MAT *B, MAT *C);

//   Dense matrix - block from transposed sparse matrix multiplication
extern MAT *mbsptr_mlt(const MAT *B, const SPMAT *A, int i0, int j0, int n0, 
		       MAT *C);

//   Block from sparse matrix - vector multiplication.
extern VEC *bspv_mlt(const SPMAT *A, int i0, int j0, int m0, const VEC *B, 
		     VEC *C)  ;

//   Check column indices of a sparse matrix
extern void check_sparse(SPMAT *C);

//   Transpose of a sparse matrix
extern SPMAT *sp_transp_ea(const SPMAT *A, SPMAT *T);

#endif

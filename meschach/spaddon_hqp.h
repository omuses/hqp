/*
 * spaddon_hqp.h
 *   - some Meschach sparse matrix add-ons and extensions
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

#ifndef SPADDON_HQP_H
#define SPADDON_HQP_H

#include "matrix.h"
#include "sparse.h"

MESCH__BEGIN_DECLS

/*
 * sprow_norm2 -- Euclidean norm of a SPROW.
 */
MESCH_API Real sprow_norm2(const SPROW *r);

/* 
 * sprow_norm1 -- L1 norm of a SPROW.
 */
MESCH_API Real sprow_norm1(const SPROW *r);

/*
 * sp_extract_mat_iv -- Copy a block from sparse matrix src to dense matrix dst
 */
MESCH_API void sp_extract_mat_iv(const SPMAT *src, const IVEC *iv, int j_offs, 
				 MAT *dst);

/*
 * bspm_mlt -- Block from sparse matrix - dense matrix multiplication
 */
MESCH_API MAT *bspm_mlt(const SPMAT *A, int i0, int j0, int m0, 
			const MAT *B, MAT *C);

/*
 * mbsptr_mlt -- Dense matrix - block from transposed sparse matrix 
 *               multiplication
 */
MESCH_API MAT *mbsptr_mlt(const MAT *B, const SPMAT *A, int i0, int j0, 
			  int n0, MAT *C);

/*
 * bspv_mlt -- Block from sparse matrix - vector multiplication
 */
MESCH_API VEC *bspv_mlt(const SPMAT *A, int i0, int j0, int m0, const VEC *B, 
			VEC *C)  ;

/*
 * check_sparse -- Check column indices of a sparse matrix
 */
MESCH_API void check_sparse(SPMAT *C);

MESCH__END_DECLS

#endif

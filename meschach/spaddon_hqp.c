/*
 * spaddon_hqp.c
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

#include <math.h>

#include "spaddon_hqp.h"
#include "sparse2.h"

/*
 * sprow_norm2 -- Euclidean norm of a SPROW
 * E. Arnold 12/12/96
 */
Real sprow_norm2(const SPROW *r)
{
    Real sum = 0.0;
    int i;

    if ( r == (SPROW *) NULL )
	m_error(E_NULL, "sprow_norm2");

    for ( i = 0; i < r->len; i++ )
	sum += r->elt[i].val*r->elt[i].val;

    return sqrt(sum);
}

/* 
 * sprow_norm1 -- L1 norm of a SPROW
 * E. Arnold 2001-08-07
 */
Real sprow_norm1(const SPROW *r)
{
    Real sum = 0.0;
    int i;

    if ( r == (SPROW *) NULL )
	m_error(E_NULL, "sprow_norm1");

    for ( i = 0; i < r->len; i++ )
	sum += fabs(r->elt[i].val);

    return sum;
}

/*
 * sp_extract_mat_iv -- Copy a block from sparse matrix src to dense matrix dst
 *                      using an integer vector iv of row indices.
 *                      dst = src(iv(:),j_offs+(1:size(dst,2))
 * E. Arnold  10/07/96 (adapted from (rf) sp_extract_mat)
 */
void sp_extract_mat_iv(const SPMAT *src, const IVEC *iv, int j_offs, MAT *dst)
{
    SPROW *row;
    int i, j;
    int i_end, j_end, j_idx;

    if ( ( src == SMNULL ) || ( iv == IVNULL ) || ( dst == MNULL ) )
	m_error(E_NULL, "sp_extract_mat_iv");
    if ( iv->dim != dst->m )
	m_error(E_SIZES, "sp_extract_mat_iv");
    m_zero(dst);

    i_end = dst->m;
    j_end = dst->n;
    for (i=0; i<i_end; i++) {
	row = src->row + iv->ive[i];
	j_idx = sprow_idx(row, j_offs);
	if (j_idx < 0) {
	    if (j_idx == -1)
		m_error(E_BOUNDS,"sp_extract_mat_iv");
	    j_idx = -(j_idx + 2);
	}
	while (j_idx < row->len) {
	    j = row->elt[j_idx].col - j_offs;
	    if (j >= j_end)
		break;
	    m_set_val(dst, i, j, row->elt[j_idx].val);
	    j_idx++;
	}
    }  
}

/*
 * bspm_mlt -- Block from sparse matrix - dense matrix multiplication
 * E. Arnold   10/19/96
 */
MAT *bspm_mlt(const SPMAT *A, int i0, int j0, int m0, const MAT *B, MAT *C)
{
    SPROW *row;
    Real val;
    int i, j, j_idx;

    if ( ( A == SMNULL ) || ( B == MNULL ) )
	m_error(E_NULL, "bspm_mlt");
    if ( ( j0+(int)B->m > A->n) || ( i0 + m0 > A->m ) )
	m_error(E_SIZES, "bspm_mlt");
    if ( B == C ) 
	m_error(E_INSITU, "bspm_mlt");

    if ( ( C == MNULL ) || ( (int)C->m != m0 ) || ( C->n != B->n ) )
	C = m_resize(C, m0, B->n);

    C = m_zero(C);
    row = A->row + i0;
    for ( i = 0; i < m0; i++, row++) {
	j_idx = sprow_idx(row, j0);
	if ( j_idx < 0 ) {
	    if ( j_idx == -1 )
		m_error(E_BOUNDS,"bspm_mlt");
	    j_idx = -(j_idx + 2);
	}
	while ( j_idx < row->len ) {
	    j = row->elt[j_idx].col-j0;
	    if ( j >= (int)B->m )
		break;
	    val = row->elt[j_idx].val;
	    /*      for ( k = 0; k < B->n; k++ ) */
	    /*	m_add_val(C, i, k, m_entry(B, j, k)*val); */
	    __mltadd__(C->me[i],B->me[j],val,(int) B->n);
	    j_idx++;
	}
    }
    return C;
}

/*
 * mbsptr_mlt -- Dense matrix - block from transposed sparse matrix 
 *               multiplication
 * E. Arnold   10/24/96
 */
MAT *mbsptr_mlt(const MAT *B, const SPMAT *A, int i0, int j0, int n0, MAT *C)
{
    SPROW *row;
    Real val;
    int i, j, j_idx, k;

    if ( ( A == SMNULL ) || ( B == MNULL ) )
	m_error(E_NULL, "mbsptr_mlt");
    if ( ( i0+(int)B->n > A->n) || ( j0 + n0 > A->m ) )
	m_error(E_SIZES, "mbsptr_mlt");
    if ( B == C ) 
	m_error(E_INSITU, "msptr_mlt");

    if ( ( C == MNULL ) || ( (int)C->n != n0 ) || ( C->m != B->m ) )
	C = m_resize(C, B->m, n0);

    C = m_zero(C);
    row = A->row + j0;
    for ( j = 0; j < n0; j++, row++) {
	j_idx = sprow_idx(row, i0);
	if ( j_idx < 0 ) {
	    if ( j_idx == -1 )
		m_error(E_BOUNDS,"msptr_mlt");
	    j_idx = -(j_idx + 2);
	}
	while ( j_idx < row->len ) {
	    i = row->elt[j_idx].col-i0;
	    if ( i >= (int) B->n )
		break;
	    val = row->elt[j_idx].val;
	    for ( k = 0; k < (int) B->m; k++ ) 
		m_add_val(C, k, j, m_entry(B, k, i)*val);
	    j_idx++;
	}
    }
    return C;
}

/*
 * bspv_mlt -- Block from sparse matrix - vector multiplication
 * E. Arnold   10/24/96
 */
VEC *bspv_mlt(const SPMAT *A, int i0, int j0, int m0, const VEC *B, VEC *C)
{
    SPROW *row;
    Real val;
    int i, j, j_idx;

    if ( ( A == SMNULL ) || ( B == VNULL ) )
	m_error(E_NULL, "bspv_mlt");
    if ( ( j0+(int)B->dim > A->n) || ( i0+m0 > A->m ) )
	m_error(E_SIZES, "bspv_mlt");
    if ( B == C ) 
	m_error(E_INSITU, "bspv_mlt");

    if ( ( C == VNULL ) || ( (int)C->dim != m0 ) )
	C = v_resize(C, m0);

    row = A->row + i0;
    for ( i = 0; i < m0; i++, row++) {
	j_idx = sprow_idx(row, j0);
	if ( j_idx < 0 ) {
	    if ( j_idx == -1 )
		m_error(E_BOUNDS,"bspv_mlt");
	    j_idx = -(j_idx + 2);
	}
	val = 0.0;
	while ( j_idx < row->len ) {
	    j = row->elt[j_idx].col-j0;
	    if ( j >= (int)B->dim )
		break;
	    val += row->elt[j_idx].val*v_entry(B, j);
	    j_idx++;
	}
	v_set_val(C, i, val);
    }
    return C;
}

/*
 * check_sparse -- Check column indices of a sparse matrix
 * E. Arnold   10/25/96
 */
void check_sparse(SPMAT *C)
{ 
    int i, j, j_idx;
    if ( ! C )
	m_error(E_NULL,"check_sparse");
    for ( i = 0; i < C->m; i++ )
	for ( j = 0; j < C->row[i].len; j++ ) {
	    j_idx = C->row[i].elt[j].col;
	    if ( ( j_idx < 0 ) || ( j_idx >= C->n ) )
		printf("check_sparse: i = %d, j = %d, j_idx = %d\n", 
		       i, j, j_idx);
	}
}


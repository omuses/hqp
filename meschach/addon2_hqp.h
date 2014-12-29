/*
 * addon2_hqp.h
 *   - additional Meschach functions for Hqp
 *
 * rf, 8/18/94
 *
 */

/*
    Copyright (C) 1994--2014  Ruediger Franke

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

#ifndef ADDON2_HQP_H
#define ADDON2_HQP_H

#include "matrix.h"
#include "sparse.h"

MESCH__BEGIN_DECLS

/** @name Additional Meschach functions */
//@{
MESCH_API IVEC *iv_set(IVEC *, int);
MESCH_API IVEC *iv_part(IVEC *iv, int offs, int dim, IVEC *header);
MESCH_API IVEC *iv_expand(IVEC *iv, int nel, int granul);
MESCH_API IVEC *iv_copy_elements(const IVEC *, IVEC *);

MESCH_API VEC *v_set(VEC *, Real);
MESCH_API VEC *v_part(VEC *v, int offs, int dim, VEC *header);
MESCH_API VEC *v_expand(VEC *v, int nel, int granul);
MESCH_API VEC *v_copy_elements(const VEC *, VEC *);
MESCH_API VEC *bd_mv_mlt(const BAND *, const VEC *, VEC *);

MESCH_API MAT *m_mltadd(const MAT *, const MAT *, const MAT *, MAT *);
MESCH_API MAT *m_copy_elements(const MAT *, MAT *);

MESCH_API Real sp_norm_inf(SPMAT *);
MESCH_API SPMAT *sp_copy3(const SPMAT *, SPMAT *);
MESCH_API int sp_update_val(SPMAT *, int, int, Real);
MESCH_API void sp_insert_mat(SPMAT *dst, int i_offs, int j_offs, const MAT *src);
MESCH_API void symsp_insert_symmat(SPMAT *dst, int offs, const MAT *src);
MESCH_API void sp_update_mat(SPMAT *dst, int i_offs, int j_offs, const MAT *src);
MESCH_API void sp_extract_mat(const SPMAT *src, int i_offs, int j_offs, MAT *dst);
MESCH_API void symsp_extract_mat(const SPMAT *src, int offs, MAT *dst);
MESCH_API void sp_insert_mrow(SPMAT *dst, int i_offs, int j_offs,
			      const MAT *src, int i);
MESCH_API void sp_update_mrow(SPMAT *dst, int i_offs, int j_offs,
			      const MAT *src, int i);
MESCH_API void sp_extract_mrow(const SPMAT *src, int i_offs, int j_offs,
			       MAT *dst, int i);
MESCH_API SPMAT *sp_ident(SPMAT *);
MESCH_API SPMAT *sp_ones(SPMAT *);
MESCH_API VEC *sp_mv_mltadd(const VEC *v1, const VEC *v2, const SPMAT *A,
			    Real alpha, VEC *out);
MESCH_API VEC *sp_vm_mltadd(const VEC *v1, const VEC *v2, const SPMAT *A,
			    Real alpha, VEC *out);
MESCH_API VEC *sp_mv_symmlt(SPMAT *A, const VEC *v, VEC *out);

MESCH_API Real sprow_inprod(const SPROW *r1, const VEC *inner, const SPROW *r2);
MESCH_API void sprow_zero(SPROW *row);

MESCH_API SPMAT *spLUfactor2(SPMAT *A, PERM *px);
MESCH_API SPMAT *sp_transp(const SPMAT *, SPMAT *);
//@}

/**
 * @name Routines for copying sparse matrices into sparse/band matrices
 *   ("sym" means that only the upper part of a symmetric matrix is filled)
 */
//@{
MESCH_API void sp_into_sp(const SPMAT *src, Real s, SPMAT *dst,
			  const PERM *px, int i_offs, int j_offs);
MESCH_API void spT_into_sp(const SPMAT *src, Real s, SPMAT *dst,
			   const PERM *px, int i_offs, int j_offs);
MESCH_API void sp_into_symsp(const SPMAT *src, Real s, SPMAT *dst,
			     const PERM *px, int i_offs, int j_offs);
MESCH_API void symsp_into_symsp(const SPMAT *src, Real s, SPMAT *dst,
				const PERM *px, int offs);
MESCH_API void spT_into_symsp(const SPMAT *src, Real s, SPMAT *dst,
			      const PERM *px, int i_offs, int j_offs);

MESCH_API void sp_into_bd(const SPMAT *sp, Real s, BAND *bd,
			  const PERM *px, int i_offs, int j_offs);
MESCH_API void spT_into_bd(const SPMAT *sp, Real s, BAND *bd,
			   const PERM *px, int i_offs, int j_offs);
//@}

MESCH__END_DECLS

#endif

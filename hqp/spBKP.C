/*
 * BKP factorization for sparse matrices
 *  - reference: J.R.Bunch, L.Kaufman, and B.N.Parlett:
 *    Decomposition of a Symmetric Matrix, Numer.Math. 27, 95--109 (1976)
 *  - factorize sparse matrix A in situ into P'AP = MDM'
 *  - M is unit lower triangular
 *  - D is block diagonal with blocks 1x1 and 2x2
 *  - delete lower diagonal part, if any (use it for fill-in)
 *  - P may only be used together with M by BKP-solve routine
 *  - always update diag-access
 *  - factor routine performs (i - j) index searches per iteration
 *    with current row i and pivot row j
 *  - solve routine needs no index searches at all
 *  - apply row-header-copy insead of element-copy whenever possible 
 *
 * rf, 9/13/94
 *
 * 7/16/95
 *  - use SPMAT for swap row for Meschach mem_info compatibility
 */

/*
    Copyright (C) 1994--2002  Ruediger Franke

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

#include <stdio.h>
#include <math.h>

#include "Meschach.h"


#define MINROWLEN 	10
#define INSERT_IDX(idx) 	(-(idx + 2))
#define MEM_MOVE(src,dst,len) 	MEM_COPY(src,dst,len)

#if 0
/*
 * memory tests
 */

static int sprow_size(SPROW *r)
{
  if (!r)
    return 0;

  return sizeof(SPROW) + r->maxlen * sizeof(row_elt);
}

static int spmat_size(SPMAT *m)
{
  if (!m)
    return 0;

  int size = sizeof(SPMAT);
  int i, i_end = m->max_m;

  for (i=0; i<i_end; i++)
    size += sprow_size(m->row + i);

  return size;
}
#endif

/* sprow_ins_val -- insert the idx's entry into sparse row r with col j
 * -- derived frome sprow_set_val (Copyright D.E. Steward, Z. Leyk)
 */
static Real sprow_ins_val(SPROW *r, int idx, Real val, int j, int type)
{
  int  new_len;
  
  if (!r)
    m_error(E_NULL,"sprow_ins_val");
  
  /* shift & insert new value */
  if (idx < 0)
    idx = INSERT_IDX(idx);
  if ( r->len >= r->maxlen ) {
    r->len = r->maxlen;
    new_len = max(2*r->maxlen+1,5);
    if (mem_info_is_on()) {
      mem_bytes(type, r->maxlen*sizeof(row_elt),
		new_len*sizeof(row_elt)); 
    }
    
    r->elt = RENEW(r->elt,new_len,row_elt);
    if ( ! r->elt )        /* can't allocate */
      m_error(E_MEM,"sprow_ins_val");
    r->maxlen = new_len;
  }
  if ( idx < r->len )
    MEM_MOVE((char *)(&(r->elt[idx])),(char *)(&(r->elt[idx+1])),
	     (r->len-idx)*sizeof(row_elt));

  r->len++;
  r->elt[idx].col = j;

  return r->elt[idx].val = val;
}

/* spbkp_mltadd -- sets r_out <- r1 + alpha.r2
 * -- derived frome sprow_mltadd (Copyright D.E. Steward, Z. Leyk)
 * -- specialized for BKP: 
 *     + r1 starts with index r1->diag
 *     + start index for r2 is passed as argument (r2_idx0)
 *     + r_out is filled from index 0; r_out->diag is set
 * -- cannot be in situ
 * -- type must be SPMAT or SPROW depending on
 *    whether r_out is a row of a SPMAT structure
 *    or a SPROW variable
 * -- returns r_out
 */
static SPROW *spbkp_mltadd(const SPROW *r1, const SPROW *r2, int r2_idx0,
			   Real s, SPROW *r_out, int type)
{
  int	idx1, idx2, idx_out, len1, len2, len_out;
  row_elt	*elt1, *elt2, *elt_out;
  
  if (!r1 || !r2)
    m_error(E_NULL,"spbkp_mltadd");
  if (r1 == r_out || r2 == r_out)
    m_error(E_INSITU,"spbkp_mltadd");
  if (!r_out)
    /* don't use sprow_get() because of Meschach memory management */
    /* r_out = sprow_get(MINROWLEN); */
    m_error(E_NULL,"spbkp_mltadd");
  
  /* Initialise */
  len1 = r1->len;
  len2 = r2->len;
  len_out = r_out->maxlen;
  idx1    = r1->diag;
  idx2    = r2_idx0;
  idx_out = 0;
  if (idx1 >= 0 || idx2 >= 0)
    r_out->diag = 0;
  else
    r_out->diag = -2;
  idx1    = (idx1 < 0)? INSERT_IDX(idx1): idx1;
  idx2    = (idx2 < 0)? INSERT_IDX(idx2): idx2;
  elt1    = &(r1->elt[idx1]);
  elt2    = &(r2->elt[idx2]);
  elt_out = &(r_out->elt[idx_out]);

  while (idx1 < len1 || idx2 < len2) {
    if (idx_out >= len_out) {
      /* r_out is too small */
      r_out->len = idx_out;
      r_out = sprow_xpd(r_out,0,type);
      len_out = r_out->maxlen;
      elt_out = &(r_out->elt[idx_out]);
    }
    if (idx2 >= len2 || (idx1 < len1 && elt1->col <= elt2->col)) {
      elt_out->col = elt1->col;
      elt_out->val = elt1->val;
      if (idx2 < len2 && elt1->col == elt2->col) {
	elt_out->val += s*elt2->val;
	elt2++;
	idx2++;
      }
      elt1++;
      idx1++;
    }
    else {
      elt_out->col = elt2->col;
      elt_out->val = s*elt2->val;
      elt2++;
      idx2++;
    }
    elt_out++;
    idx_out++;
  }
  r_out->len = idx_out;
  
  return r_out;
}

/* 
 * spbkp_mltadd2 -- sets rowj := rowj + s.row + t.row1
 */
static SPROW *spbkp_mltadd2(SPROW *rowj, const SPROW *row, int idx, Real s,
			    const SPROW *row1, int idx1, Real t,
			    SPROW *swap, int type)
{
  spbkp_mltadd(rowj, row, idx, s, swap, type);
  spbkp_mltadd(swap, row1, idx1, t, rowj, type);

  return rowj;
}

static void interchange(SPMAT *A, int i, int j,
			const IVEC *col_idxs, SPROW *swap)
/*
 * -- interchange row and col j of A and row and col i
 * -- i < j
 * -- A is the reduced matrix of order n-i+1
 * -- build new row i in swap and exchange headers at the end
 */
{
  int	n, k, len, col;
  int	idx, idxi, idxj;
  SPROW *row, *rowi, *rowj, tmp_row;
  row_elt *elt;
  Real 	aii;

  n = A->n;
  rowi = A->row + i;
  rowj = A->row + j;
  idxi = rowi->diag;
  idxj = rowj->diag;

  swap->len = len = 0;

  /*
   * build swap (new row i)
   */

  /* diagonal entry j to swap (i) */
  if (idxj >= 0) {
    /* see if swap needs expanding */
    if (swap->maxlen == 0) {
      swap = sprow_xpd(swap, 0, TYPE_SPMAT);
    }
    swap->elt[0].val = rowj->elt[idxj].val;
    swap->elt[0].col = i;
    swap->diag = 0;
    len = 1;
    idxj ++;
  }
  else {
    swap->diag = -2;
    idxj = INSERT_IDX(idxj);
  }

  /* backing store diagonal entry of row i */
  if (idxi >= 0 && idxi < rowi->len) {
    aii = rowi->elt[idxi].val;
    idxi ++;
  }
  else {
    aii = 0.0;
    idxi = INSERT_IDX(idxi);
  }

  /* exchange column j and row i for i<j */
  k = i+1;
  row = A->row + k;
  if (idxi < rowi->len)
    col = rowi->elt[idxi].col;
  else
    col = n;
  for (; k < j; k++, row++) {
    idx = col_idxs->ive[k];

    /* col j to swap (row i) */
    if (idx >= 0) {
      /* see if swap needs expanding */
      if (swap->maxlen == len) {
	swap->len = len;
	swap = sprow_xpd(swap, 0, TYPE_SPMAT);
      }
      elt = swap->elt + len;
      elt->val = row->elt[idx].val;
      elt->col = k;
      len ++;
    }

    /* row i to col j */
    if (k == col && idxi < rowi->len) {
      sp_set_val(A,k,j,rowi->elt[idxi].val);
      idxi ++;
      if (idxi < rowi->len)
	col = rowi->elt[idxi].col;
      else
	col = n;
    }
    else if (idx >= 0) {
      row->elt[idx].val = 0.0;
    }
  }

  /* element (i,j) from row i to swap */
  if (col == j && idxi < rowi->len) {
    /* see if swap needs expanding */
    if (swap->maxlen == len) {
      swap->len = len;
      swap = sprow_xpd(swap, 0, TYPE_SPMAT);
    }
    elt = swap->elt + len;
    elt->val = rowi->elt[idxi].val;
    elt->col = j;
    len ++;
    idxi ++;
  }

  /* row j (off the diagonal) onto swap (row i) */
  k = rowj->len - idxj;
  if (k > 0) {
    if (swap->maxlen < len + k) {
      swap->len = len;
      swap = sprow_xpd(swap, len + k, TYPE_SPMAT);
    }
    MEM_COPY((char *)(rowj->elt + idxj), (char *)(swap->elt + len),
	     k * sizeof(row_elt));
    len += k;
  }

  swap->len = len;

  /* 
   * build new row j
   */

  /* diagonal entry aii to j */
  if (aii != 0.0) {
    /* see if rowj needs expanding */
    if (rowj->maxlen == 0) {
      rowj = sprow_xpd(rowj, 0, TYPE_SPMAT);
    }
    rowj->elt[0].val = aii;
    rowj->elt[0].col = j;
    rowj->diag = 0;
    len = 1;
  }
  else {
    rowj->diag = -2;
    len = 0;
  }

  /* resting row i onto row j */
      
  k = rowi->len - idxi;
  if (k > 0) {
    if (rowj->maxlen < len + k) {
      rowj->len = len;
      rowj = sprow_xpd(rowj, len + k, TYPE_SPMAT);
    }
    MEM_COPY((char *)(rowi->elt + idxi), (char *)(rowj->elt + len),
	     k * sizeof(row_elt));
    len += k;
  }

  rowj->len = len;

  /*
   * exchange swap and old row i
   */

  MEM_COPY(rowi, &tmp_row, sizeof(SPROW));
  MEM_COPY(swap, rowi, sizeof(SPROW));
  MEM_COPY(&tmp_row, swap, sizeof(SPROW));
}


SPMAT *spBKPfactor(SPMAT *A, PERM *pivot, Real tol)
/*
 * -- factorize A in situ into P'AP = MDM'
 * -- P(i+1) == 0 for (i,i+1) is a 2x2 block
 */
{
  int	i, ip1, ip2, j, j1, k, k_end, n;
  int  	idx, idx1, len;
  Real	aii, aip1, aiip1, lambda, sigma, tmp;
  Real	det, s, t;
  SPMAT *swap_mat;
  SPROW *row, *row1, *swap, tmp_row;
  row_elt *elt, *elt1;
  IVEC	*col_idxs;
  Real alpha;

  if (!A || !pivot)
    m_error(E_NULL, "spBKPfactor");
  if (A->n != (int)pivot->size || A->n != A->m)
    m_error(E_SIZES, "spBKPfactor");
  if (tol < 0.0 || tol > 1.0 )
    m_error(E_RANGE, "spBKPfactor");

  alpha = tol * 0.6403882032022076; /* = tol * (1+sqrt(17))/8 */

  n = A->n;
  px_ident(pivot);
  col_idxs = iv_get(n);
  /* don't use sprow_get because of Meschach memory management */
  /* swap = sprow_get(MINROWLEN); */
  swap_mat = sp_get(1, n, MINROWLEN);
  swap = swap_mat->row;

  if (!A->flag_diag)
    sp_diag_access(A);
  A->flag_col = 0;

  for (i = 0; i < n-1;) {
    row = A->row + i;
    ip1 = i+1;
    ip2 = i+2;

    /*
     * find the maximum element in the first column of the reduced
     * matrix below the diagonal (go through first row as A is symmetric)
     */
    idx = row->diag;
    if (idx >= 0) {
      elt = row->elt + idx;
      aii = fabs(elt->val);
      elt ++;
      idx ++;
    }
    else {
      idx = INSERT_IDX(idx);
      elt = row->elt + idx;
      aii = 0.0;
    }

    lambda = aii;
    j = i;
    k_end = row->len;
    for (; idx < k_end; idx++, elt++) {
      tmp = fabs(elt->val);
      if (tmp > lambda) {
	j = elt->col;
	lambda = tmp;
      }
    }
    if (aii >= alpha * lambda)
      goto onebyone;
    
    /* 
     * - determine the maximum element in the jth column of the 
     *   reduced matrix off the diagonal
     * - cache col indizes for interchange
     */
    sigma = lambda;
    k = ip1;
    row = A->row + k;
    for (; k < j; k++, row++) {
      idx = sprow_idx(row, j);
      col_idxs->ive[k] = idx;
      if (idx >= 0) {
	tmp = fabs(row->elt[idx].val);
	if (tmp > sigma)
	  sigma = tmp;
      }
    }
    row = A->row + j;
    idx = row->diag;
    if (idx >= 0)
      idx ++;
    else
      idx = INSERT_IDX(idx);
    elt = row->elt + idx;
    k_end = row->len;
    for (; idx < k_end; idx++, elt++) {
      tmp = fabs(elt->val);
      if (tmp > sigma)
	sigma = tmp;
    }
    if (sigma * aii >= alpha * lambda * lambda)
      goto onebyone;

    row = A->row + j;
    idx = row->diag;
    if (idx >= 0)
      tmp = fabs(row->elt[idx].val);
    else
      tmp = 0.0;
    if (tmp >= alpha * sigma) {
      interchange(A, i, j, col_idxs, swap);
      pivot->pe[i] = j;
      goto onebyone;
    }

    /*
     * do 2x2 pivot
     */
    row = A->row + i;
    /* get element aii */
    idx = row->diag;
    if (idx >= 0) {
      aii = row->elt[idx].val;
      idx ++;
    }
    else {
      aii = 0.0;
      idx = INSERT_IDX(idx);
    }
    /* perform row/col interchange; get element aiip1 */
    elt = row->elt + idx;
    if (idx < row->len && elt->col == ip1)
      aiip1 = elt->val;
    else
      aiip1 = 0.0;
    if (j > ip1) {
      interchange(A, ip1, j, col_idxs, swap);
      /* exchange A[i][ip1] with A[i][j] */
      /* update A[i][j] first */
      idx1 = sprow_idx(row, j);
      if (idx1 >= 0) {
	tmp = row->elt[idx1].val;
	row->elt[idx1].val = aiip1;
      }
      else {
	tmp = 0.0;
	/* fill in aiip1 */
	if (aiip1 != 0.0)
	  sprow_ins_val(row, idx1, aiip1, j, TYPE_SPMAT);
      }
      /* now update A[i][ip1] */
      aiip1 = tmp;

      elt = row->elt + idx;
      if (idx < row->len && elt->col == ip1) {
	elt->val = aiip1;
      }
      else {
	/* fill in aiip1 */
	if (aiip1 != 0.0)
	  sprow_ins_val(row, idx, aiip1, ip1, TYPE_SPMAT);
      }
    }
    /* set idx to next element after A[i][ip] */
    if (idx < row->len && row->elt[idx].col == ip1) {
      idx ++;
    }
    pivot->pe[i] = j;
    pivot->pe[ip1] = 0;
    /* get A[ip1][ip1], set idx1 behind it */
    row1 = A->row + ip1;
    idx1 = row1->diag;
    if (idx1 >= 0) {
      aip1 = row1->elt[idx1].val;
      idx1 ++;
    }
    else {
      aip1 = 0.0;
      idx1 = INSERT_IDX(idx1);
    }
    det = aii * aip1 - aiip1 * aiip1;
    aii /= det;
    aiip1 /= det;
    aip1 /= det;

    j = (idx < row->len)? row->elt[idx].col: n;
    j1 = (idx1 < row1->len)? row1->elt[idx1].col: n;
    while (j < n || j1 < n) {
      if (j < j1) {
	elt = row->elt + idx;
	s = aip1 * elt->val;
	t = - aiip1 * elt->val;
	spbkp_mltadd2(A->row+j, row, idx, -s,
		      row1, INSERT_IDX(idx1), -t, swap, TYPE_SPMAT);
	elt->val = s;
	idx ++;
	if (t != 0.0) {
	  sprow_ins_val(row1, idx1, t, j, TYPE_SPMAT);
	  idx1 ++;
	}
      }
      else if (j1 < j) {
	elt1 = row1->elt + idx1;
	s = - aiip1 * elt1->val;
	t = aii * elt1->val;
	spbkp_mltadd2(A->row+j1, row, INSERT_IDX(idx), -s,
		      row1, idx1, -t, swap, TYPE_SPMAT);
	if (s != 0.0) {
	  sprow_ins_val(row, idx, s, j1, TYPE_SPMAT);
	  idx ++;
	}
	elt1->val = t;
	idx1 ++;
      }
      else { /* j == j1 */
	elt = row->elt + idx;
	elt1 = row1->elt + idx1;
	s = - aiip1 * elt1->val + aip1 * elt->val;
	t = - aiip1 * elt->val + aii * elt1->val;
	spbkp_mltadd2(A->row+j, row, idx, -s,
		      row1, idx1, -t, swap, TYPE_SPMAT);
	elt->val = s;
	idx ++;
	elt1->val = t;
	idx1 ++;
      }
      j = (idx < row->len)? row->elt[idx].col: n;
      j1 = (idx1 < row1->len)? row1->elt[idx1].col: n;
    }

    i = ip2;
    continue;

  onebyone:
    /*
     * do 1x1 pivot
     */
    row = A->row + i;
    idx = row->diag;
    if (idx >= 0) {
      aii = row->elt[idx].val;
      idx ++;
    }
    else {
      aii = 0.0;
      idx = INSERT_IDX(idx);
    }
    if (aii != 0.0) {
      elt = row->elt + idx;
      len = row->len;
      for (; idx < len; idx++, elt++) {
	j = elt->col;
	s = elt->val / aii;
	if (s != 0.0) {
	  spbkp_mltadd(A->row+j, row, idx, -s, swap, TYPE_SPMAT);
	  /* exchange swap with row #j */
	  MEM_COPY(A->row+j, &tmp_row, sizeof(SPROW));
	  MEM_COPY(swap, A->row+j, sizeof(SPROW));
	  MEM_COPY(&tmp_row, swap, sizeof(SPROW));
	}
	elt->val = s;
      }
    }

    i = ip1;
  }

  iv_free(col_idxs);
  /* sprow_free(swap); */
  sp_free(swap_mat);

  A->flag_diag = 1;

  return A;
}

VEC *spBKPsolve(const SPMAT *A, const PERM *pivot, const VEC *b, VEC *x)
/*
 * - solve A*x = b after A has been factorized by spBKPfactor
 * - raise an E_SING if A is singular
 */	
{
  int i, ii, k, ip1;
  int n, idx, idx1, len;
  Real det, tmp, save;
  Real aiip1, aii, aip1;
  Real *x_ve;
  u_int *p_pe;
  SPROW *row, *row1;
  row_elt *elt;

  if (!A || !pivot || !b) 
    m_error(E_NULL, "spBKPsolve");
  if (!A->flag_diag)
    m_error(E_FORMAT, "spBKPsolve");
  n = A->n;
  if ((int)b->dim != n || (int)pivot->size != n || A->n != A->m)
    m_error(E_SIZES, "spBKPsolve");
  if (!x || (int)x->dim != n)
    x = v_resize(x,n);

  p_pe = pivot->pe;

  /*
   * solve MDy = b for y, where b is stored in x and store y in x
   */
  
  x = v_copy(b,x);
  x_ve = x->ve;
  for (i = 0; i < n-1;) {
    row = A->row + i;
    row1 = row + 1;
    ip1 = i+1;
    save = x_ve[p_pe[i]];
    idx = row->diag;
    if (idx >= 0) {
      aii = row->elt[idx].val;
      idx ++;
    }
    else {
      aii = 0.0;
      idx = INSERT_IDX(idx);
    }
    if (p_pe[ip1] > 0) {
      /*
       * 1x1 pivot
       */
      x_ve[p_pe[i]] = x_ve[i];
      if (aii == 0.0)
	m_error(E_SING, "spBKPsolve");
      x_ve[i] = save / aii;
      elt = row->elt + idx;
      len = row->len;
      for (; idx < len; idx++, elt++) {
	x_ve[elt->col] -= save * elt->val;
      }
      i = ip1;
    }
    else {
      /*
       * 2x2 pivot
       */
      tmp = x_ve[i];
      x_ve[p_pe[i]] = x_ve[ip1];
      if (row->elt[idx].col == ip1) {
	aiip1 = row->elt[idx].val;
	idx ++;
      }
      else
	aiip1 = 0.0;
      idx1 = row1->diag;
      if (idx1 >= 0) {
	aip1 = row1->elt[idx1].val;
	idx1 ++;
      }
      else {
	aip1 = 0.0;
	idx1 = INSERT_IDX(idx1);
      }
      det = aii * aip1 - aiip1 * aiip1;
      if (det == 0.0)
	m_error(E_SING, "spBKPsolve");
      x_ve[i] = (tmp * aip1 - save * aiip1) / det;
      x_ve[ip1] = (save * aii - tmp * aiip1) / det;
      elt = row->elt + idx;
      len = row->len;
      for (; idx < len; idx++, elt++) {
	x_ve[elt->col] -= tmp * elt->val;
      }
      elt = row1->elt + idx1;
      len = row1->len;
      for (; idx1 < len; idx1++, elt++) {
	x_ve[elt->col] -= save * elt->val;
      }
      i += 2;
    }
  }
  if (i == n-1) {
    row = A->row + i;
    idx = row->diag;
    if (idx >= 0)
      aii = row->elt[idx].val;
    else
      aii = 0.0;
    if (aii == 0.0)
      m_error(E_SING, "spBKPsolve");
    x_ve[i] /= aii;
    i = n - 2;
  }
  else
    i = n - 3;

  /*
   * solve M'x = y for x, where y is stored in x
   */
  while (i >= 0) {
    ip1 = i + 1;
    if (p_pe[i] > 0 || i == 0)
      ii = i;       /*  1x1 pivot */
    else
      ii = i-1;     /*  2x2 pivot */
    for (k = ii; k <= i; k++) {
      tmp = x_ve[k];
      row = A->row + k;
      idx = row->diag;
      idx = (idx < 0)? INSERT_IDX(idx): idx;
      elt = row->elt + idx;
      len = row->len;
      /* go to index i+1 */
      while (idx < len && elt->col < ip1) {
	idx ++;
	elt ++;
      }
      for (; idx < len; idx++, elt++)
	tmp -= elt->val * x_ve[elt->col];
      x_ve[k] = tmp;
    }
    if (i != (int)p_pe[ii]) {
      tmp = x_ve[i];
      x_ve[i] = x_ve[p_pe[ii]];
      x_ve[p_pe[ii]] = tmp;
    }
    i = ii - 1;
  }

  return x;
}

/*
 * Meschach.C -- definitions
 *
 * rf, 8/18/94
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

#include <math.h>
#include "Meschach.h"

#define	SPLU_MOD

#define SPROW_IDX(r,k) \
   ((k < r->elt[0].col) ? -2 : \
    ((k > r->elt[r->len-1].col) ? -(r->len+2) : sprow_idx(r,k)))

/* 
 * spLUfactor2 -- derived from spLUfactor, faster for banded matrices
 *   - makro to avoid many binary searches in sprow_idx (+25%)
 *   - no Markowitz strategy as A is globaly permuted (+15%)
 *   - additional continue-check for elimination (+10%)
 */
SPMAT	*spLUfactor2(SPMAT *A, PERM *px)
{
  	int	i, best_i, k, idx, m, n;
	SPROW	*r, *r_piv, tmp_row;
	static	SPROW	*merge = (SPROW *)NULL;
	Real	max_val, tmp;
	static VEC	*col_vals=VNULL;
#ifdef SPLU_MOD
  	int	piv_max_col;
#else
	int	len, best_len;
	double	alpha;
	alpha = 0.1;
#endif
	if ( ! A || ! px )
		m_error(E_NULL,"spLUfctr");
#ifndef SPLU_MOD
	if ( alpha <= 0.0 || alpha > 1.0 )
		m_error(E_RANGE,"alpha in spLUfctr");
#endif
	if ( (int)px->size <= A->m )
		px = px_resize(px,A->m);
	px_ident(px);
	col_vals = v_resize(col_vals,A->m);
	MEM_STAT_REG(col_vals,TYPE_VEC);

	m = A->m;	n = A->n;
#ifndef SPLU_MOD
	if ( ! A->flag_col )
		sp_col_access(A);
	if ( ! A->flag_diag )
		sp_diag_access(A);
#endif
	A->flag_col = A->flag_diag = FALSE;
	if ( ! merge ) {
	   merge = sprow_get(20);
	   MEM_STAT_REG(merge,TYPE_SPROW);
	}

	for ( k = 0; k < n; k++ )
	{
	    /* find pivot row/element for partial pivoting */

	    /* get first row with a non-zero entry in the k-th column */
	    max_val = 0.0;
	    best_i = k;
	    for ( i = k; i < m; i++ )
	    {
		r = &(A->row[i]);
		idx = SPROW_IDX(r,k);
		if ( idx < 0 )
		    tmp = 0.0;
		else {
		    tmp = r->elt[idx].val;
		    if ( fabs(tmp) > max_val ) {
		        max_val = fabs(tmp);
		        best_i = i;
		    }
		}
		col_vals->ve[i] = tmp;
	    }

	    if ( max_val == 0.0 )
		continue;
#ifndef SPLU_MOD
	    best_len = n+1;	/* only if no possibilities */
	    best_i = -1;
	    for ( i = k; i < m; i++ )
	    {
		tmp = fabs(col_vals->ve[i]);
		if ( tmp == 0.0 )
		    continue;
		if ( tmp >= alpha*max_val )
		{
		    r = &(A->row[i]);
		    idx = SPROW_IDX(r,k);
		    len = (r->len) - idx;
		    if ( len < best_len )
		    {
			best_len = len;
			best_i = i;
		    }
		}
	    }
#endif
	    /* swap row #best_i with row #k */
	    MEM_COPY(&(A->row[best_i]),&tmp_row,sizeof(SPROW));
	    MEM_COPY(&(A->row[k]),&(A->row[best_i]),sizeof(SPROW));
	    MEM_COPY(&tmp_row,&(A->row[k]),sizeof(SPROW));
	    /* swap col_vals entries */
	    tmp = col_vals->ve[best_i];
	    col_vals->ve[best_i] = col_vals->ve[k];
	    col_vals->ve[k] = tmp;
	    px_transp(px,k,best_i);

	    r_piv = &(A->row[k]);
#ifdef SPLU_MOD
	    piv_max_col = r_piv->elt[r_piv->len - 1].col;
#endif
	    for ( i = k+1; i < n; i++ )
	    {
		r = &(A->row[i]);
#ifdef SPLU_MOD
		if (r->len == 0 || piv_max_col < r->elt[0].col)
		  continue;
#endif
		/* compute and set multiplier */
		tmp = col_vals->ve[i]/col_vals->ve[k];
		if ( tmp != 0.0 )
		    sp_set_val(A,i,k,tmp);
		else
		    continue;

		/* perform row operations */
		merge->len = 0;
		sprow_mltadd(r,r_piv,-tmp,k+1,merge,TYPE_SPROW);
		idx = SPROW_IDX(r,k+1);
		if ( idx < 0 )
		    idx = -(idx+2);
		/* see if r needs expanding */
		if ( r->maxlen < idx + merge->len )
		    sprow_xpd(r,idx+merge->len,TYPE_SPMAT);
		r->len = idx+merge->len;
		MEM_COPY((char *)(merge->elt),(char *)&(r->elt[idx]),
			merge->len*sizeof(row_elt));
	    }
	}

	return A;
}

//--------------------------------------------------------------------------
Real sprow_inprod(const SPROW *r1, const VEC *inner, const SPROW *r2)
{
  Real sum;
  int j_idx1, j_idx2, j1, j2;
  int len1, len2;
  row_elt *elt1, *elt2;

  j_idx1=0;
  j_idx2=0;
  len1 = r1->len;
  len2 = r2->len;
  elt1 = r1->elt;
  elt2 = r2->elt;
  sum = 0.0;
  while (j_idx1 < len1 && j_idx2 < len2) {
    j1 = elt1->col;
    j2 = elt2->col;
    if (j1 < j2) {
      j_idx1 ++;
      elt1 ++;
    }
    else if (j1 > j2) {
      j_idx2 ++;
      elt2 ++;
    }
    else {
      sum += elt1->val * inner->ve[j1] * elt2->val;
      j_idx1 ++;
      elt1 ++;
      j_idx2 ++;
      elt2 ++;
    }
  }

  return sum;
}

//-------------------------------------------------------------------------
void sprow_zero(SPROW *row)
{
  int i, len = row->len;
  row_elt *elt = row->elt;

  for (i = 0; i < len; i++, elt++)
    elt->val = 0.0;
}

//-------------------------------------------------------------------------
VEC *v_part(VEC *v, int offs, int dim, VEC *head)
{
  if (!v)
    m_error(E_NULL,"v_part");
  if ((int)v->dim < offs + dim)
    m_error(E_SIZES,"v_part");

  head->dim = dim;
  head->max_dim = 0;
  head->ve = v->ve + offs;

  return head;
}

//-------------------------------------------------------------------------
VEC *v_set(VEC *v, Real val)
{
  int i, i_end;
  Real *ve;

  if (v != VNULL) {
    ve = v->ve;
    i_end = v->dim;
    for (i=0; i<i_end; i++)
      *ve++ = val;
  }

  return v;
}

//-------------------------------------------------------------------------
VEC *v_expand(VEC *v, int nel, int granul)
{
  VEC *vbak;
  VEC vv;

  granul = max(granul, nel);
    
  if (!v) {
    v = v_get(granul);
    v->dim = nel;
  } else {
    if (v->dim + nel <= v->max_dim)
      v->dim += nel;
    else {
      vbak = v_copy(v, VNULL);
      v_resize(v, v->dim + granul);
      v_copy(vbak, v_part(v, 0, vbak->dim, &vv));
      v->dim = vbak->dim + nel;
      v_free(vbak);
    }
  }
    
  return v;
}

//-------------------------------------------------------------------------
IVEC *iv_part(IVEC *iv, int offs, int dim, IVEC *head)
{
  if (!iv)
    m_error(E_NULL,"iv_part");
  if ((int)iv->dim < offs + dim)
    m_error(E_SIZES,"iv_part");

  head->dim = dim;
  head->max_dim = 0;
  head->ive = iv->ive + offs;

  return head;
}

//-------------------------------------------------------------------------
IVEC *iv_set(IVEC *iv, int val)
{
  int i, i_end;
  int *ive;

  if (iv != IVNULL) {
    ive = iv->ive;
    i_end = iv->dim;
    for (i=0; i<i_end; i++)
      *ive++ = val;
  }

  return iv;
}

//-------------------------------------------------------------------------
IVEC *iv_expand(IVEC *iv, int nel, int granul)
{
  IVEC *ivbak;
  IVEC ivv;

  granul = max(granul, nel);
    
  if (!iv) {
    iv = iv_get(granul);
    iv->dim = nel;
  } else {
    if (iv->dim + nel <= iv->max_dim)
      iv->dim += nel;
    else {
      ivbak = iv_copy(iv, IVNULL);
      iv_resize(iv, iv->dim + granul);
      iv_copy(ivbak, iv_part(iv, 0, ivbak->dim, &ivv));
      iv->dim = ivbak->dim + nel;
      iv_free(ivbak);
    }
  }
    
  return iv;
}

//-------------------------------------------------------------------------
VEC *bd_mv_mlt(const BAND *A, const VEC *x, VEC *out)
{
  int i, j, j_end, k;
  int start_idx, end_idx;
  int n, m, lb, ub;
  Real **A_me;
  Real *x_ve;
  Real sum;

  if (!A || !x)
    m_error(E_NULL,"bd_mv_mlt");
  if (x->dim != A->mat->n)
    m_error(E_SIZES,"bd_mv_mlt");
  if (!out || out->dim != A->mat->n)
    out = v_resize(out, A->mat->n);
  if (out == x)
    m_error(E_INSITU,"bd_mv_mlt");

  n = A->mat->n;
  m = A->mat->m;
  lb = A->lb;
  ub = A->ub;
  A_me = A->mat->me;
  start_idx = lb;
  end_idx = m + n-1 - ub;
  for (i=0; i<n; i++, start_idx--, end_idx--) {
    j = max(0, start_idx);
    k = max(0, -start_idx);
    j_end = min(m, end_idx);
    x_ve = x->ve + k;
    sum = 0.0;	     
    for (; j<j_end; j++, k++)
      sum += A_me[j][k] * *x_ve++;
    out->ve[i] = sum;
  }

  return out;
}

//-------------------------------------------------------------------------
// m_mltadd	-- out = in + a*b 
//
MAT *m_mltadd(const MAT *IN, const MAT *A, const MAT *B, MAT *OUT)
{
  u_int	i, k, m, n, p;
  Real	**A_v, **B_v;

  if (A == (MAT *)NULL || B == (MAT *)NULL || IN == (MAT *)NULL)
    m_error(E_NULL,"m_mltadd");
  if (A->n != B->m || IN->m != A->m || IN->n != B->n)
    m_error(E_SIZES,"m_mltadd");
  if (A == OUT || B == OUT)
    m_error(E_INSITU,"m_mltadd");

  if (IN != OUT)
    OUT = m_copy(IN, OUT);

  m = A->m;
  n = A->n;
  p = B->n;
  A_v = A->me;
  B_v = B->me;

  for (i = 0; i < m; i++)
    for (k = 0; k < n; k++) {
      if (A_v[i][k] != 0.0)
	__mltadd__(OUT->me[i], B_v[k], A_v[i][k], (int)p);
    }

  return OUT;
}

//-------------------------------------------------------------------------
SPMAT *sp_copy3(const SPMAT *src, SPMAT *dst)
{
  SPROW	*rs, *rd;
  int	i, n, m;

  if (!src)
    m_error(E_NULL, "sp_copy3");
  n = src->n;
  m = src->m;
  if (!dst || dst->n != n || dst->m != m)
    dst = sp_resize(dst, m, n);

  rd = dst->row;
  rs = src->row;
  for (i=0; i<m; i++, rd++, rs++) {
    if (rd->maxlen < rs->len)
      sprow_resize(rd, rs->len, TYPE_SPMAT);
    rd->len = rs->len;
    rd->diag = rs->diag;
    MEM_COPY(rs->elt, rd->elt, rs->len * sizeof(row_elt));
  }

  dst->flag_col = src->flag_col;
  dst->flag_diag = src->flag_diag;
  
  return dst;
}

//-------------------------------------------------------------------------
// sp_update_val:
//  set an entry of a sparse matrix that should already exist
//  return  0  if existing entry updated
//          1  if new entry allocated
//
int sp_update_val(SPMAT *A, int i, int j, Real val)
{
  SPROW	*row;
  int	idx;
   
  if ( A == SMNULL )
    m_error(E_NULL,"sp_update_val");
  if ( i < 0 || i >= A->m || j < 0 || j >= A->n )
    m_error(E_SIZES,"sp_update_val");
   
  row = A->row+i;
  idx = sprow_idx(row, j);
  if ( idx >= 0 ) {
    row->elt[idx].val = val;
    return 0;
  }
  else {
    // allocate new entry
    sp_set_val(A, i, j, val);
    return -1;
  }
}

//---------------------------------------------------------------------------
// sp_insert_mrow:
//  insert all non-zero entries of a dense matrix--row into a sparse matrix
//     
void sp_insert_mrow(SPMAT *dst, int i_offs, int j_offs, const MAT *src, int i)
{
  int j, jend = src->n;
  double val;

  if (!dst || !src)
    m_error(E_NULL, "sp_insert_mrow");
  if (i < 0 || i >= (int)src->m)
    m_error(E_BOUNDS, "sp_insert_mrow");

  for(j=0; j<jend; j++) {
    val = src->me[i][j];
    if (val != 0.0)
      sp_set_val(dst, i+i_offs, j+j_offs, val);
  }
}

//---------------------------------------------------------------------------
// sp_update_mrow:
//  update entries of a sparse matrix with elements of a dense matrix--row
//
void sp_update_mrow(SPMAT *dst, int i_offs, int j_offs, const MAT *src, int i)
{
  SPROW	*row;
  int	j, j_idx, j_end;

  if (!dst || !src)
    m_error(E_NULL, "sp_update_mrow");
  if (i < 0 || i >= (int)src->m)
    m_error(E_BOUNDS, "sp_update_mrow");

  j_end = src->n;
  row = dst->row + i_offs + i;
  j_idx = sprow_idx(row, j_offs);
  if (j_idx < 0) {
    if (j_idx == -1)
      m_error(E_BOUNDS,"sp_update_mrow");
    j_idx = -(j_idx + 2);
  }
  while (j_idx < row->len) {
    j = row->elt[j_idx].col - j_offs;
    if (j >= j_end)
      break;
    row->elt[j_idx].val = src->me[i][j];
    j_idx++;
  }
}

//----------------------------------------------------------------------------
void sp_extract_mrow(const SPMAT *src, int i_offs, int j_offs, MAT *dst, int i)
{
  SPROW *row;
  int j, j_end, j_idx;

  if (!src || !dst)
    m_error(E_NULL, "sp_extract_mrow");
  if (i < 0 || i >= (int)dst->m)
    m_error(E_BOUNDS, "sp_extract_mrow");

  j_end = dst->n;
  for (j=0; j<j_end; j++)
    dst->me[i][j] = 0;

  row = src->row + i_offs + i;
  j_idx = sprow_idx(row, j_offs);
  if (j_idx < 0) {
    if (j_idx == -1)
      m_error(E_BOUNDS,"sp_extract_mrow");
    j_idx = -(j_idx + 2);
  }
  while (j_idx < row->len) {
    j = row->elt[j_idx].col - j_offs;
    if (j >= j_end)
      break;
    dst->me[i][j] = row->elt[j_idx].val;
    j_idx++;
  }
}

//--------------------------------------------------------------------------
// sp_insert_mat:
//  insert all non-zero entries of a dense matrix into a sparse matrix
//     
void sp_insert_mat(SPMAT *dst, int i_offs, int j_offs, const MAT *src)
{
  int i, iend = src->m;
  int j, jend = src->n;
  double val;

  for (i=0; i<iend; i++)
    for(j=0; j<jend; j++) {
      val = src->me[i][j];
      if (val != 0.0)
	sp_set_val(dst, i+i_offs, j+j_offs, val);
    }
}

//--------------------------------------------------------------------------
// symsp_insert_symmat:
//  insert upper diagonal non-zero entries of a symmetric dense matrix 
//  into a symmetric sparse matrix
//     
void symsp_insert_symmat(SPMAT *dst, int offs, const MAT *src)
{
  int i, iend = src->m;
  int j, jend = src->n;
  double val;

  for (i = 0; i < iend; i++)
    for(j = i; j < jend; j++) {
      val = src->me[i][j];
      if (val != 0.0)
	sp_set_val(dst, i + offs, j + offs, val);
    }
}

//--------------------------------------------------------------------------
// sp_update_mat:
//  update entries of a sparse matrix with elements of a dense matrix
//
void sp_update_mat(SPMAT *dst, int i_offs, int j_offs, const MAT *src)
{
  SPROW	*row;
  int	i, i_end;
  int	j, j_idx, j_end;

  i_end = src->m;
  j_end = src->n;
  for(i=0; i<i_end; i++) {
    row = dst->row + i_offs + i;
    j_idx = sprow_idx(row, j_offs);
    if (j_idx < 0) {
      if (j_idx == -1)
	m_error(E_BOUNDS,"sp_update_mat");
      j_idx = -(j_idx + 2);
    }
    while (j_idx < row->len) {
      j = row->elt[j_idx].col - j_offs;
      if (j >= j_end)
	break;
      row->elt[j_idx].val = src->me[i][j];
      j_idx++;
    }
  }
}

//-------------------------------------------------------------------------
void sp_extract_mat(const SPMAT *src, int i_offs, int j_offs, MAT *dst)
{
  SPROW *row;
  int i, j;
  int i_end, j_end, j_idx;

  m_zero(dst);

  i_end = dst->m;
  j_end = dst->n;
  for (i=0; i<i_end; i++) {
    row = src->row + i_offs + i;
    j_idx = sprow_idx(row, j_offs);
    if (j_idx < 0) {
      if (j_idx == -1)
	m_error(E_BOUNDS,"sp_extract_mat");
      j_idx = -(j_idx + 2);
    }
    while (j_idx < row->len) {
      j = row->elt[j_idx].col - j_offs;
      if (j >= j_end)
	break;
      dst->me[i][j] = row->elt[j_idx].val;
      j_idx++;
    }
  }  
}

//-------------------------------------------------------------------------
void symsp_extract_mat(const SPMAT *src, int offs, MAT *dst)
{
  SPROW *row;
  Real val;
  int i, j;
  int i_end, j_end, j_idx;

  m_zero(dst);

  i_end = dst->m;
  j_end = dst->n;
  for (i = 0; i < i_end; i++) {
    row = src->row + offs + i;
    j_idx = sprow_idx(row, offs + i);
    if (j_idx < 0) {
      if (j_idx == -1)
	m_error(E_BOUNDS,"sp_extract_mat");
      j_idx = -(j_idx + 2);
    }
    while (j_idx < row->len) {
      j = row->elt[j_idx].col - offs;
      if (j >= j_end)
	break;
      val = row->elt[j_idx].val;
      dst->me[i][j] = val;
      if (i != j)
	dst->me[j][i] = val;
      j_idx++;
    }
  }  
}

//-------------------------------------------------------------------------
SPMAT *sp_ident(SPMAT *A)
{
  int i;
  int i_end = min(A->n, A->m);

  sp_zero(A);
  for (i=0; i<i_end; i++)
    sp_set_val(A, i, i, 1.0);

  return A;
}

/* sp_norm_inf -- \max_i \sum_j |a_{ij}| */
Real sp_norm_inf(SPMAT *A)
{
  int	i, idx, len;
  row_elt	*elt;
  Real rsum, norm = 0.0;
  
  if ( ! A )
    m_error(E_NULL,"sp_norm_inf");
  
  for (i = 0; i < A->m; i++) {
    elt = A->row[i].elt;
    len = A->row[i].len;
    rsum = 0.0;
    for (idx = 0; idx < len; idx++, elt++)
      rsum += fabs((*elt).val);
    norm = max(norm, rsum);
  }	
  
  return norm;
}

/* sp_ones -- set all the (represented) elements of a sparse matrix to 1.0 */
SPMAT	*sp_ones( SPMAT *A)
{
  int	i, idx, len;
  row_elt	*elt;
  
  if ( ! A )
    m_error(E_NULL,"sp_ones");
  
  for ( i = 0; i < A->m; i++ )
    {
      elt = A->row[i].elt;
      len = A->row[i].len;
      for ( idx = 0; idx < len; idx++ )
	(*elt++).val = 1.0;
    }
  
  return A;
}

/*
 * sp_transp:
 *   -- build the transpose of sparse matrix A in T
 *   -- return T
 *   -- routine may not work in-situ
 */
SPMAT *sp_transp(const SPMAT *A, SPMAT *T)
{
  int n, m;
  int i, j_idx;
  int slen, dlen;
  row_elt *selt, *delt;
  SPROW *srow, *drow;

  if (!A)
    m_error(E_NULL, "sp_transp");
  if (A == T)
    m_error(E_INSITU, "sp_transp");

  n = A->n;
  m = A->m;

  if (!T)
    T = sp_copy(A);
  if (T->m != n || T->n != m)
    T = sp_resize(T, n, m);

  for (i=0; i<n; i++)
    T->row[i].len = 0;
  T->flag_col = 0;
  T->flag_diag = 0;

  for (i=0; i<m; i++) {
    srow = A->row + i;
    slen = srow->len;
    selt = srow->elt;
    for (j_idx = 0; j_idx < slen; j_idx++, selt++) {
      drow = T->row + selt->col;
      dlen = drow->len;
#if 1
      if (dlen == drow->maxlen) {
	drow = sprow_xpd(drow, (3*dlen)/2, TYPE_SPMAT);
	drow->len = dlen;
      }
#else
      if (dlen == drow->maxlen) {
	//	drow = sprow_xpd(drow, (3*dlen)/2, TYPE_SPMAT);
	drow = sprow_xpd(drow, dlen+10, TYPE_SPMAT);
	drow->len = dlen;
      }
#endif
      delt = drow->elt + dlen;
      delt->val = selt->val;
      delt->col = i;
      drow->len ++;
    }      
  }

  return T;
}

/* sp_mv_mltadd -- sparse matrix/dense vector multiply and add
 *  --  returns out == v1 + alpha*A*v2
 *  --  if out==NULL on entry then the result vector is created
 */
VEC	*sp_mv_mltadd(const VEC *v1, const VEC *v2, const SPMAT *A,
		      Real alpha, VEC *out)
{
  int	i, j_idx, m, n, max_idx;
  Real	sum, *v2_ve;
  SPROW	*r;
  row_elt	*elts;
  
  if ( !A || !v1 || !v2 )
    m_error(E_NULL,"sp_mv_mltadd");
  if ( (int)v1->dim != A->m || (int)v2->dim != A->n )
    m_error(E_SIZES,"sp_mv_mltadd");
  if ( !out || (int)out->dim != A->m )
    out = v_resize(out, A->m);

  if (out != v1)
    out = v_copy(v1, out);

  m = A->m;
  n = A->n;
  v2_ve = v2->ve;
  
  for ( i = 0; i < m; i++ )
    {
      sum = 0.0;
      r = &(A->row[i]);
      max_idx = r->len;
      elts    = r->elt;
      for ( j_idx = 0; j_idx < max_idx; j_idx++, elts++ )
	sum += elts->val*v2_ve[elts->col];
      out->ve[i] += alpha * sum;
    }

  return out;
}

/* sp_mv_symmlt -- symmetric sparse matrix/dense vector multiply
 *  --  only upper part of A is used (uA is unit upper diagonal)
 *  --  returns out == (uA' + diag(A) + uA) * v
 *  --  if out==NULL on entry then the result vector is created
 */
VEC	*sp_mv_symmlt(SPMAT *A, const VEC *v, VEC *out)
{
  int	i, j_idx, m, n, max_idx;
  Real	sum, tmp, *v_ve, *out_ve;
  SPROW	*r;
  row_elt *elts;
  
  if ( !A || !v )
    m_error(E_NULL, "sp_mv_symmlt");
  if ( (int)v->dim != A->m )
    m_error(E_SIZES, "sp_mv_symmlt");
  if ( !out || (int)out->dim != A->m )
    out = v_resize(out, A->m);
  if (out == v)
    m_error(E_INSITU, "sp_mv_symmlt");

  if (!A->flag_diag)
    sp_diag_access(A);

  m = A->m;
  n = A->n;
  v_ve = v->ve;
  out_ve = out->ve;

  v_zero(out);
  for ( i = 0; i < m; i++ )
    {
      tmp = v_ve[i];
      r = &(A->row[i]);
      max_idx = r->len;
      j_idx = r->diag;
      if (j_idx < 0) {
	sum = 0.0;
	j_idx = -(j_idx + 2);
      }
      else {
	sum = r->elt[j_idx].val * v_ve[i];
	j_idx ++;
      }
      elts = &(r->elt[j_idx]);
      for (; j_idx < max_idx; j_idx++, elts++ ) {
	sum += elts->val * v_ve[elts->col];
	out_ve[elts->col] += elts->val * tmp;
      }
      out_ve[i] += sum;
    }

  return out;
}

/* sp_vm_mltadd -- vector-matrix multiply and add
 *	-- may not be in situ (out != v2)
 *	-- returns out' == v1' + alpha*v2'*A
 */
VEC *sp_vm_mltadd(const VEC *v1, const VEC *v2, const SPMAT *A,
		  Real alpha, VEC *out)
{
  int	i_idx, j, m, len;
  Real	tmp, *out_ve;
  SPROW	*row;
  row_elt *elt;

  if ( ! v1 || ! v2 || ! A )
    m_error(E_NULL,"sp_vm_mltadd");
  if ( v2 == out )
    m_error(E_INSITU,"sp_vm_mltadd");
  if ( (int)v1->dim != A->n || A->m != (int)v2->dim )
    m_error(E_SIZES,"sp_vm_mltadd");
  if ( !out || (int)out->dim != A->n )
    out = v_resize(out, A->n);

  if (out != v1)
    out = v_copy(v1, out);

  out_ve = out->ve;
  m = A->m;
  for (j=0; j<m; j++) {
    tmp = v2->ve[j]*alpha;
    if ( tmp != 0.0 ) {
      row = &(A->row[j]);
      len = row->len;
      elt = row->elt;
      for (i_idx=0; i_idx<len; i_idx++, elt++)
	out_ve[elt->col] += elt->val * tmp;
    }
  }
  return out;
}

/*
 * routines for copying sparse matrices into band matrices
 */

//-------------------------------------------------------------------------
void sp_into_bd(const SPMAT *sp, Real s, BAND *bd, const PERM *px,
		 int i_offs, int j_offs)
{
  int i, i_bd, i_end;
  int j_bd, j_idx, j_end;
  const row_elt *elt;
  int lb = bd->lb;
  Real **bde = bd->mat->me;

  i_end = sp->m;
  for (i=0; i<i_end; i++) {
    i_bd = px->pe[i + i_offs];
    elt = sp->row[i].elt;
    j_end = sp->row[i].len;
    for (j_idx=0; j_idx<j_end; j_idx++, elt++) {
      j_bd = px->pe[elt->col + j_offs];
      bde[lb + j_bd - i_bd][j_bd] = s * elt->val;
    }
  }
}

//-------------------------------------------------------------------------
void spT_into_bd(const SPMAT *sp, Real s, BAND *bd, const PERM *px,
		 int i_offs, int j_offs)
{
  int j, j_bd, j_end;
  int i_bd, i_idx, i_end;
  const row_elt *elt;
  int lb = bd->lb;
  Real **bde = bd->mat->me;

  j_end = sp->m;
  for (j=0; j<j_end; j++) {
    j_bd = px->pe[j + j_offs];
    elt = sp->row[j].elt;
    i_end = sp->row[j].len;
    for (i_idx=0; i_idx<i_end; i_idx++, elt++) {
      i_bd = px->pe[elt->col + i_offs];
      bde[lb + j_bd - i_bd][j_bd] = s * elt->val;
    }
  }
}

//-------------------------------------------------------------------------
void sp_into_sp(const SPMAT *src, Real s, SPMAT *dst, const PERM *px,
		int i_offs, int j_offs)
{
  int i, i_dst, i_end;
  int j_dst, j_idx, j_end;
  const row_elt *elt;

  i_end = src->m;
  for (i=0; i<i_end; i++) {
    i_dst = px->pe[i + i_offs];
    elt = src->row[i].elt;
    j_end = src->row[i].len;
    for (j_idx=0; j_idx<j_end; j_idx++, elt++) {
      j_dst = px->pe[elt->col + j_offs];
      sp_set_val(dst, i_dst, j_dst, s * elt->val);
    }
  }
}

//-------------------------------------------------------------------------
void spT_into_sp(const SPMAT *src, Real s, SPMAT *dst, const PERM *px,
		 int i_offs, int j_offs)
{
  int j, j_dst, j_end;
  int i_dst, i_idx, i_end;
  const row_elt *elt;

  j_end = src->m;
  for (j=0; j<j_end; j++) {
    j_dst = px->pe[j + j_offs];
    elt = src->row[j].elt;
    i_end = src->row[j].len;
    for (i_idx=0; i_idx<i_end; i_idx++, elt++) {
      i_dst = px->pe[elt->col + i_offs];
      sp_set_val(dst, i_dst, j_dst, s * elt->val);
    }
  }
}

//-------------------------------------------------------------------------
void sp_into_symsp(const SPMAT *src, Real s, SPMAT *dst, const PERM *px,
		   int i_offs, int j_offs)
{
  int i, i_dst, i_end;
  int j_dst, j_idx, j_end;
  const row_elt *elt;

  i_end = src->m;
  for (i=0; i<i_end; i++) {
    i_dst = px->pe[i + i_offs];
    elt = src->row[i].elt;
    j_end = src->row[i].len;
    for (j_idx=0; j_idx<j_end; j_idx++, elt++) {
      j_dst = px->pe[elt->col + j_offs];
      if (i_dst <= j_dst)
	sp_set_val(dst, i_dst, j_dst, s * elt->val);
    }
  }
}

//-------------------------------------------------------------------------
void symsp_into_symsp(const SPMAT *src, Real s, SPMAT *dst, const PERM *px,
		      int offs)
{
  int i, j, i_dst, i_end;
  int j_dst, j_idx, j_end;
  const row_elt *elt;

  i_end = src->m;
  for (i = 0; i < i_end; i++) {
    i_dst = px->pe[i + offs];
    elt = src->row[i].elt;
    j_end = src->row[i].len;
    for (j_idx=0; j_idx<j_end; j_idx++, elt++) {
      j = elt->col;
      /* take only diagonal and upper diagonal entries */
      if (j >= i) {
	j_dst = px->pe[j + offs];
	if (j_dst >= i_dst)
	  sp_set_val(dst, i_dst, j_dst, s * elt->val);
	else
	  sp_set_val(dst, j_dst, i_dst, s * elt->val);
      }
    }
  }
}

//-------------------------------------------------------------------------
void spT_into_symsp(const SPMAT *src, Real s, SPMAT *dst, const PERM *px,
		    int i_offs, int j_offs)
{
  int j, j_dst, j_end;
  int i_dst, i_idx, i_end;
  const row_elt *elt;

  j_end = src->m;
  for (j=0; j<j_end; j++) {
    j_dst = px->pe[j + j_offs];
    elt = src->row[j].elt;
    i_end = src->row[j].len;
    for (i_idx=0; i_idx<i_end; i_idx++, elt++) {
      i_dst = px->pe[elt->col + i_offs];
      if (i_dst <= j_dst)
	sp_set_val(dst, i_dst, j_dst, s * elt->val);
    }
  }
}


//=========================================================================

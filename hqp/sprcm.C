/*
 * Reverse Cuthill McKee ordering of sparse matrices and KKT systems.
 *
 * rf, 11/15/95
 */

/*
    Copyright (C) 1994--2000  Ruediger Franke

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

#include <stdlib.h>

#include "Meschach.h"

/*
 * define a struct to describe a node
 */
typedef struct Node {
  int node;	/* position of the node in the original matrix */
  int deg;	/* degree of the node */
  int deg2;	/* sum of degrees of the node's neighbours */
} Node;

/*
 * compare two nodes
 */
static int cmpNode(const Node *n1, const Node *n2)
{
  int ret;
  ret = n1->deg - n2->deg;
  if (ret == 0) {
    ret = n1->deg2 - n2->deg2;
  }
  return ret;
}

/*
 * sp_rcm_scan:
 *   -- analyze sparsity pattern of a Kuhn-Tucker-system
 *      Q A'C'
 *      A
 *      C
 *   -- A or C may be NULL
 *   -- only the upper diagonal part is inspected
 *   -- allocate and return concatted list of neighbours of all nodes
 */
IVEC *sp_rcm_scan(const SPMAT *Q, const SPMAT *A, const SPMAT *C,
		  IVEC *degree, IVEC *neigh_start, IVEC *neighs)
{
  int i, i_end, j, j_idx, j_end, n, offs;
  int *d_ive, *s_ive, *n_ive;
  SPROW *row;
  row_elt *elt;

  if (!Q || !degree || !neigh_start)
    error(E_NULL, "sp_rcm_scan");
  if (Q->n != Q->m)
    error(E_SIZES, "sp_rcm_scan");

  n = Q->n;
  if (A)
    n += A->m;
  if (C)
    n += C->m;

  if ((int)degree->dim != n || (int)neigh_start->dim != n + 1)
    error(E_SIZES, "sp_rcm_scan");

  /*
   * first pass:
   *   -- go through Q,A,C and count neighbours of each node (degree)
   *   -- initialize neigh_start to point into neighs
   *   -- allocate neighs
   */

  iv_zero(degree);

  d_ive = degree->ive;
  s_ive = neigh_start->ive;

  row = Q->row;
  i_end = Q->m;
  for (i = 0; i < i_end; i++, row++) {
    elt = row->elt;
    j_end = row->len;
    for (j_idx = 0; j_idx < j_end; j_idx++, elt++) {
      j = elt->col;
      if (j > i) {
	d_ive[i] ++;
	d_ive[j] ++;
      }
    }
  }

  if (A) {
    offs = Q->n;
    row = A->row;
    i_end = offs + A->m;
    for (i = offs; i < i_end; i++, row++) {
      elt = row->elt;
      j_end = row->len;
      for (j_idx = 0; j_idx < j_end; j_idx++, elt++) {
	j = elt->col;
	d_ive[i] ++;
	d_ive[j] ++;
      }
    }
  }

  if (C) {
    offs = Q->n;
    if (A)
      offs += A->m;
    row = C->row;
    i_end = offs + C->m;
    for (i = offs; i < i_end; i++, row++) {
      elt = row->elt;
      j_end = row->len;
      for (j_idx = 0; j_idx < j_end; j_idx++, elt++) {
	j = elt->col;
	d_ive[i] ++;
	d_ive[j] ++;
      }
    }
  }

  offs = 0;
  for (i = 0; i < n; i++) {
    s_ive[i] = offs;
    offs += d_ive[i];
    d_ive[i] = 0;
  }
  s_ive[n] = offs;

  neighs = iv_resize(neighs, offs);
  n_ive = neighs->ive;

  /*
   * second pass:
   *   -- initialize neighs as a concatted list of neighbours of all nodes
   */

  row = Q->row;
  i_end = Q->m;
  for (i = 0; i < i_end; i++, row++) {
    elt = row->elt;
    j_end = row->len;
    for (j_idx = 0; j_idx < j_end; j_idx++, elt++) {
      j = elt->col;
      if (j > i) {
	n_ive[s_ive[i] + d_ive[i]] = j;
	n_ive[s_ive[j] + d_ive[j]] = i;
	d_ive[i] ++;
	d_ive[j] ++;
      }
    }
  }

  if (A) {
    offs = Q->n;
    row = A->row;
    i_end = offs + A->m;
    for (i = offs; i < i_end; i++, row++) {
      elt = row->elt;
      j_end = row->len;
      for (j_idx = 0; j_idx < j_end; j_idx++, elt++) {
	j = elt->col;
	n_ive[s_ive[i] + d_ive[i]] = j;
	n_ive[s_ive[j] + d_ive[j]] = i;
	d_ive[i] ++;
	d_ive[j] ++;
      }
    }
  }

  if (C) {
    offs = Q->n;
    if (A)
      offs += A->m;
    row = C->row;
    i_end = offs + C->m;
    for (i = offs; i < i_end; i++, row++) {
      elt = row->elt;
      j_end = row->len;
      for (j_idx = 0; j_idx < j_end; j_idx++, elt++) {
	j = elt->col;
	n_ive[s_ive[i] + d_ive[i]] = j;
	n_ive[s_ive[j] + d_ive[j]] = i;
	d_ive[i] ++;
	d_ive[j] ++;
      }
    }
  }

  return neighs;
}

/*
 * Find Reverse-Cuthill-McKee order for a sparse matrix
 * given by a concatted list of neighbours of all nodes.
 *  degree:       number of neighbours for n nodes
 *  neigh_start:  start index into neighs for each node i, 0 <= i <= n
 *                (neigh_start->ive[n] is the first free index)
 *  neighs:       concatted list of neighbours of n nodes
 *  order:        determined permutation
 *
 * Reference:
 *  J. Weissinger. Sp"arlich besetzte Gleichungssysteme,
 *  Wissenschaftsverlag, Mannheim, 1990.
 */
PERM *sp_rcm_order(const IVEC *degree, const IVEC *neigh_start,
		   const IVEC *neighs, PERM *order)
{
  int n, root;
  int count, node;
  char *marks;
  char *glob_marks;
  Node *levels;
  int *n_ive, *s_ive, *d_ive;
  int l_begin, l_end;
  int nlevels, nlevels_old;
  int i, j, j_end, k, k_end;
  int deg;
  int cluster_start;

  if (!degree || !neigh_start || !neighs)
    error(E_NULL, "sp_rcm_order");
  n = degree->dim;
  if (n + 1 != (int)neigh_start->dim
      || neigh_start->ive[n] != (int)neighs->dim)
    error(E_SIZES, "sp_rcm_order");

  if (!order || (int)order->size != n)
    order = px_resize(order, n);

  marks = (char *)malloc(n);
  glob_marks = (char *)malloc(n);
  levels = (Node *)calloc(n, sizeof(Node));
  if (!marks || !glob_marks || !levels)
    error(E_MEM, "sp_rcm_order");

  n_ive = neighs->ive;
  s_ive = neigh_start->ive;
  d_ive = degree->ive;
  
  for (i = 0; i < n; i++)
    glob_marks[i] = 1;
  root = 0;
  count = 0;

  while (root < n) {

    nlevels = 0;
    cluster_start = count;

    do { /* ... while (nlevels > nlevels_old) */

      /* find level structure for current root node */

      count = cluster_start;

      for (i = 0; i < n; i++)
	marks[i] = glob_marks[i];

      nlevels_old = nlevels;
      nlevels = 0;
      l_begin = count;
      l_end = count;
      levels[count++].node = root;	/* add root node */
      marks[root] = 0;			/* unmark root node */
      /* decrement degree of neighbours */
      k_end = s_ive[root + 1];
      for (k = s_ive[root]; k < k_end; k++) {
	d_ive[n_ive[k]]--;
      }

      do { /* while (count > l_end) */

	/* add nodes for the next level */
	
	l_begin = l_end;
	l_end = count;
	nlevels ++;
	for (i = l_begin; i < l_end; i++) {
	  node = levels[i].node;
	  j_end = s_ive[node + 1];
	  for (j = s_ive[node]; j < j_end; j++) {
	    node = n_ive[j];
	    if (marks[node]) {
	      levels[count++].node = node;	/* add node */
	      marks[node] = 0;			/* unmark added node */
	      /* decrement degree of neighbours */
	      k_end = s_ive[node + 1];
	      for (k = s_ive[node]; k < k_end; k++) {
		d_ive[n_ive[k]]--;
	      }
	    }
	  }

	  /*
	   * complete nodes of current level with degrees and sort them
	   */

	  if (count - l_end > 1) {
	    for (k = l_end; k < count; k++) {
	      node = levels[k].node;
	      levels[k].deg = d_ive[node];
	      deg = 0;
	      j_end = s_ive[node + 1];
	      for (j = s_ive[node]; j < j_end; j++) {
		if (marks[n_ive[j]]) {
		  deg += d_ive[n_ive[j]];
		}
	      }
	      levels[k].deg2 = deg;
	    }

	    qsort( (void *)(levels + l_end),
		  count - l_end, sizeof(Node),
		  (int (*)(const void *, const void *))cmpNode);
	  }
	}
      } while (count > l_end);

      if (nlevels_old > nlevels)
	warning(WARN_UNKNOWN, "sp_rcm_order: Levels decreased!");

      /*
       * choose the node with the smallest degree in the last level
       * as root for the next iteration
       */
      root = levels[l_begin].node;
      deg = d_ive[root];
      for (i = l_begin + 1; i < l_end; i++) {
	if (d_ive[levels[i].node] < deg) {
	  root = levels[i].node;
	  deg = d_ive[root];
	}
      }

      /* restore degree */
      for (i = 0; i < n; i++)
	d_ive[i] = s_ive[i + 1] - s_ive[i];

    } while (nlevels > nlevels_old);

    root = n;
    for (i = 0; i < n; i++) {
      glob_marks[i] = marks[i];
      if (marks[i])
	root = i;
    }
  }

  /*
   * fill order
   *  -- order->pe[i] says, where to put row/col i into the reordered matrix
   *  -- order is reverted
   */

  for (i = 0; i < n; i++)
    order->pe[levels[i].node] = n - i - 1;

  free(levels);
  free(glob_marks);
  free(marks);

  return order;
}

/*
 * sp_rcm_sbw:
 *   -- return semi band width according order
 *      (the number of upper diagonals)
 */
int sp_rcm_sbw(const IVEC *neigh_start, const IVEC *neighs, const PERM *order)
{
  int i, i_idx, n;
  int j_idx, j_end;
  int max_neigh;
  int ub;
  int node, dist;
  int *n_ive, *s_ive;

  n = order->size;
  n_ive = neighs->ive;
  s_ive = neigh_start->ive;

  ub = 0;
  for (i_idx = 0; i_idx < n; i_idx++) {
    i = order->pe[i_idx];
    max_neigh = i;
    j_end = s_ive[i_idx + 1];
    for (j_idx = s_ive[i_idx]; j_idx < j_end; j_idx++) {
      node = order->pe[n_ive[j_idx]];
      if (node > max_neigh)
	max_neigh = node;
    }
    dist = max_neigh - i;
    if (dist > ub)
      ub = dist;
  }

  return ub;
}

/*
 * sp_symrcm:
 *   -- analyze sparsity pattern of symmetric matrix Q
 *      and return Reverse-Cuthil-McKee order
 *   -- only the upper diagonal part of Q is inspected
 */
PERM *sp_symrcm(const SPMAT *Q, PERM *order)
{
  int n;
  IVEC *degree, *neigh_start, *neighs;

  if (!Q)
    error(E_NULL, "sp_symrcm");
  n = Q->n;
  if (n != Q->m)
    error(E_SIZES, "sp_symrcm");

  degree = iv_get(n);
  neigh_start = iv_get(n + 1);
  neighs = sp_rcm_scan(Q, SMNULL, SMNULL, degree, neigh_start, IVNULL);

  order = sp_rcm_order(degree, neigh_start, neighs, order);

  iv_free(degree);
  iv_free(neigh_start);
  iv_free(neighs);

  return order;
}

/*
 * sp_kktrcm:
 *   -- analyze sparsity pattern of symmetric Kuhn-Tucker-system
 *      and return Reverse-Cuthil-McKee order
 *   -- only the upper diagonal part is inspected
 */
PERM *sp_kktrcm(const SPMAT *Q, const SPMAT *A, const SPMAT *C,
		PERM *order)
{
  int n;
  IVEC *degree, *neigh_start, *neighs;

  if (!Q)
    error(E_NULL, "sp_kktrcm");
  if (Q->n != Q->m)
    error(E_SIZES, "sp_kktrcm");

  n = Q->n;
  if (A)
    n += A->m;
  if (C)
    n += C->m;

  degree = iv_get(n);
  neigh_start = iv_get(n + 1);
  neighs = sp_rcm_scan(Q, A, C, degree, neigh_start, IVNULL);

  order = sp_rcm_order(degree, neigh_start, neighs, order);

  iv_free(degree);
  iv_free(neigh_start);
  iv_free(neighs);

  return order;
}

/*
 * sp_copy_rcm:
 *   -- copy src into dst
 *   -- the structure of dst is destroyed (in contrast to sp_copy2)
 */
static SPMAT *sp_copy_rcm(const SPMAT *src, SPMAT *dst)
{
  SPROW	*rs, *rd;
  int	i, n, m;

  if (!src)
    error(E_NULL, "sp_copy_rcm");
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

/*
 * pxinv_spcols:
 *   -- permute columns of a sparse matrix
 *      dst = src * px'
 */

static int cmp_row_elt(const row_elt *elt1, const row_elt *elt2)
{
  int ret;
  ret = elt1->col - elt2->col;
  if (ret == 0)
    /* invalid multiple occurence of same column index */
    error(E_INTERN, "pxinv_spcols");
  return ret;
}

SPMAT *pxinv_spcols(const PERM *px, const SPMAT *src, SPMAT *dst)
{
  int i, j_idx, len, m;
  SPROW *row;
  row_elt *elt;

  if (!src || !px)
    error(E_NULL, "pxinv_spcols");
  m = src->m;
  if (m != (int)px->size)
    error(E_SIZES, "pxinv_spcols");

  /*
   * copy src into dst
   */

  if (dst != src) {
    if (!dst)
      dst = sp_copy(src);
    else
      dst = sp_copy_rcm(src, dst);
  }

  /*
   * perform permutation in situ
   * (insert new column indezes and sort row entries)
   */

  dst->flag_diag = 0;

  row = dst->row;
  for (i = 0; i < m; i++, row++) {
    len = row->len;
    elt = row->elt;
    for (j_idx = 0; j_idx < len; j_idx++, elt++)
      elt->col = px->pe[elt->col];
    if (len > 1)
      qsort(row->elt, len, sizeof(row_elt),
	    (int (*)(const void *, const void *))cmp_row_elt);
  }    

  return dst;
}

/*
 * pxinv_sprows:
 *   -- permute rows of a sparse matrix
 *      dst = px' * src
 */
SPMAT *pxinv_sprows(const PERM *px, const SPMAT *src, SPMAT *dst)
{
  int start, size, i, i_new, m;
  SPROW tmp, tmp_new;

  if (!src || !px)
    error(E_NULL, "pxinv_sprows");
  m = src->m;
  size = px->size;
  if (m != size)
    error(E_SIZES, "pxinv_sprows");

  /*
   * copy src into dst
   */

  if (dst != src) {
    if (!dst)
      dst = sp_copy(src);
    else
      dst = sp_copy_rcm(src, dst);
  }

  if (size == 0)
    return dst;

  /*
   * perform permutation in situ
   */

  dst->flag_col = dst->flag_diag = 0;

  start = 0;
  do {
    tmp_new = dst->row[start];
    i_new = start;
    do {
      i = i_new;
      i_new = px->pe[i];
      if (i_new >= size)
	/* invalid multiple occurence of same row index */
	error(E_INTERN, "pxinv_sprows");
      if (i != i_new) {
	tmp = tmp_new;
	tmp_new = dst->row[i_new];
	dst->row[i_new] = tmp;
      }
      px->pe[i] += size;
    } while (i_new != start);
    for (start = 0; (int)px->pe[start] >= size && start < size; start++);
  } while (start < size);

  for (i = 0; i < size; i++)
    px->pe[i] -= size;

  return dst;
}

/*
 * pxinv_vec:
 *   -- permute a vector
 *      dst = px' * src
 */
VEC *pxinv_vec(const PERM *px, const VEC *src, VEC *dst)
{
  int start, size, i, i_new, dim;
  Real tmp, tmp_new;

  if (!px || !src)
    error(E_NULL, "pxinv_vec");
  dim = src->dim;
  size = px->size;
  if (dim != size)
    error(E_SIZES, "pxinv_vec");

  /*
   * copy src into dst
   */

  if (dst != src) {
    dst = v_copy(src, dst);
    dst->dim = src->dim;
  }

  /*
   * perform permutation in situ
   */

  if (size == 0)
    return dst;

  start = 0;
  do {
    tmp_new = dst->ve[start];
    i_new = start;
    do {
      i = i_new;
      i_new = px->pe[i];
      if (i_new >= size)
	/* invalid multiple occurence of same row index */
	error(E_INTERN, "pxinv_vec");
      if (i != i_new) {
	tmp = tmp_new;
	tmp_new = dst->ve[i_new];
	dst->ve[i_new] = tmp;
      }
      px->pe[i] += size;
    } while (i_new != start);
    for (start = 0; (int)px->pe[start] >= size && start < size; start++);
  } while (start < size);

  for (i = 0; i < size; i++)
    px->pe[i] -= size;

  return dst;
}

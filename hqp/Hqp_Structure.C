/*
 * Hqp_Structure.C
 *
 * rf, 5/19/94
 */

/*
    Copyright (C) 1994--1998  Ruediger Franke

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

#include <assert.h>
#include <stdlib.h>

#include "Hqp_Structure.h"
#include "Hqp_Program.h"


//------------------------------------------------------------------------
// some static stuff for sorting

static int *_criteria;
static int *_criteria2;
static int _compare(const int *p1, const int *p2)
{
  int ret;
  ret = _criteria[*p1] - _criteria[*p2];
  if (ret == 0) {
    ret = _criteria2[*p1] - _criteria2[*p2];
  }
  return ret;
}

//------------------------------------------------------------------------
Hqp_Structure::Hqp_Structure()
{
  _n = 0;
  _nentries = 0;
  _neigh_size = 0;
  _order = PNULL;
}

//------------------------------------------------------------------------
Hqp_Structure::~Hqp_Structure()
{
  if (_n > 0) {
    delete [] _neigh_start;
    delete [] _degree;
  }
  if (_neigh_size > 0) {
    delete [] _neighbours;
  }
  px_free(_order);
}

//------------------------------------------------------------------------
void Hqp_Structure::realloc(int n, int neigh_size)
{
  assert (n > 0 && neigh_size > 0);

  if (n != _n) {
    if (_n > 0) {
      delete [] _neigh_start;
      delete [] _degree;
    }
    _n = n;
    _neigh_start = new int [n+1];
    _degree = new int [n];
    _order = px_resize(_order, n);
  }

  if (neigh_size != _neigh_size) {
    if (_neigh_size > 0) {
      delete [] _neighbours;
    }
    _neigh_size = neigh_size;
    _neighbours = new int [neigh_size];
  }
}

//------------------------------------------------------------------------
void Hqp_Structure::neigh_grow()
{
  int cur_size = _neigh_size;
  int *cur_neigh = _neighbours;

  _neigh_size += _neigh_size / 2;
  _neighbours = new int [_neigh_size];
  for (int i=0; i<cur_size; i++)
    _neighbours[i] = cur_neigh[i];
  delete [] cur_neigh;
}

//------------------------------------------------------------------------
void Hqp_Structure::neigh_sort()
{
  int i;
  int nneigh;

  //
  // sort neighbours according their degrees
  //

  _criteria = _degree;
  _criteria2 = _degree;
  for (i=0; i<_n; i++) {
    nneigh = _neigh_start[i+1] - _neigh_start[i];
    if (nneigh > 1)
      qsort( (void *)(_neighbours + _neigh_start[i]),
	    nneigh, sizeof(int),
	    (int (*)(const void *, const void *))_compare);
  }
}

//------------------------------------------------------------------------
int Hqp_Structure::add_row(const SPMAT *mat, int i, int node, int offs)
{
  int j_idx, jend;
  const row_elt *elt;
  int neigh_idx = _neigh_start[node+1];

  elt = mat->row[i].elt;
  jend = mat->row[i].len;
  for (j_idx=0; j_idx<jend; j_idx++, elt++) {
    if (elt->col + offs != node) {
      if (neigh_idx == _neigh_size)
	neigh_grow();
      _neighbours[neigh_idx++] = elt->col + offs;
      _degree[node] ++;
    }
  }
  _neigh_start[node+1] = neigh_idx;

  return jend;
}
#ifdef SPARSE_COL_ACCESS
//------------------------------------------------------------------------
int Hqp_Structure::add_col(const SPMAT *mat, int j, int node, int offs)
{
  int i, j_idx;
  const row_elt *elt;
  int neigh_idx = _neigh_start[node+1];
  int nadd = 0;

  i = mat->start_row[j];
  j_idx = mat->start_idx[j];
  while (i >= 0) {
    elt = &(mat->row[i].elt[j_idx]);
    if (i + offs != node) {
      if (neigh_idx == _neigh_size)
	neigh_grow();
      _neighbours[neigh_idx++] = i + offs;
      _degree[node] ++;
    }
    i = elt->nxt_row;
    j_idx = elt->nxt_idx;
    nadd++;
  }
  _neigh_start[node+1] = neigh_idx;

  return nadd;
}
#endif
//------------------------------------------------------------------------
void Hqp_Structure::init_spmat(const SPMAT *M)
{
  int i, n;

  // reallocate data

  n = M->m;
  realloc(n, 10*n);

  // scan M's sparsity structure (M should be symmetric)

  _nentries = 0;
  _neigh_start[0] = 0;
  for (i=0; i<n; i++) {
    _neigh_start[i+1] = _neigh_start[i];
    _degree[i] = 0;

    _nentries += add_row(M, i, i);
  }

  // sort neighbours

  neigh_sort();
}
#ifdef SPARSE_COL_ACCESS
//------------------------------------------------------------------------
// init_QAC -- suppose symmetric structure:
//  
//   Q A' C'
//   A
//   C   W/Z
//
void Hqp_Structure::init_QAC(const SPMAT *Q, const SPMAT *A,
			     const SPMAT *C)
{
  int node, i, iend;
  int n;

  assert(Q && A && C);
  assert(Q->m == Q->n && Q->m == A->n && Q->m == C->n);

  if (!A->flag_col)
    sp_col_access(A);
  if (!C->flag_col)
    sp_col_access(C);

  // reallocate data

  n = Q->m + A->m + C->m;
  realloc(n, 10*n);

  //
  // scan sparsity structure
  //

  _nentries = C->m;	// values of W/Z that are not looked for

  // submatrix  Q A' C'

  node = 0;
  _neigh_start[0] = 0;
  iend = Q->m;
  for (i=0; i<iend; i++, node++) {
    _neigh_start[node+1] = _neigh_start[node];
    _degree[node] = 0;

    _nentries += add_row(Q, i, node);
    _nentries += add_col(A, i, node, Q->n);
    _nentries += add_col(C, i, node, Q->n + A->m);
  }

  // matrix A by rows

  iend = A->m;
  for (i=0; i<iend; i++, node++) {
    _neigh_start[node+1] = _neigh_start[node];
    _degree[node] = 0;

    _nentries += add_row(A, i, node);
  }  

  // matrix C by rows

  iend = C->m;
  for (i=0; i<iend; i++, node++) {
    _neigh_start[node+1] = _neigh_start[node];
    _degree[node] = 0;

    _nentries += add_row(C, i, node);
  }  

  //
  // sort neighbours
  //

  neigh_sort();
}

//------------------------------------------------------------------------
// init_QA -- suppose symmetric structure:
//  
//   Q A'
//   A
//
void Hqp_Structure::init_QA(const SPMAT *Q, const SPMAT *A)
{
  int node, i, iend;
  int n;

  assert(Q && A);
  assert(Q->m == Q->n && Q->m == A->n);

  if (!A->flag_col)
    sp_col_access(A);

  // reallocate data

  n = Q->m + A->m;
  realloc(n, 10*n);

  //
  // scan sparsity structure
  //

  _nentries = 0;

  // submatrix  Q A'

  node = 0;
  _neigh_start[0] = 0;
  iend = Q->m;
  for (i=0; i<iend; i++, node++) {
    _neigh_start[node+1] = _neigh_start[node];
    _degree[node] = 0;

    _nentries += add_row(Q, i, node);
    _nentries += add_col(A, i, node, Q->n);
  }

  // matrix A by rows

  iend = A->m;
  for (i=0; i<iend; i++, node++) {
    _neigh_start[node+1] = _neigh_start[node];
    _degree[node] = 0;

    _nentries += add_row(A, i, node);
  }  

  //
  // sort neighbours
  //

  neigh_sort();
}
#endif
//------------------------------------------------------------------------
void Hqp_Structure::order_rcm()
{
  int root = 0;
  int count, node;
  char *marks = new char [_n];
  char *glob_marks = new char [_n];
  int *deg = new int [_n];
  int *deg2 = new int [_n];
  u_int *levels = _order->pe;
  int l_begin, l_end;
  int nlevels, nlevels_old;
  int i, j, jend, k, kend;
  int min_deg;
  int cluster_start;

  _criteria = deg;
  _criteria2 = deg2;

  for (i=0; i<_n; i++)
    glob_marks[i] = 1;
  root = 0;
  count = 0;

  while (root < _n) {

    nlevels = 0;
    cluster_start = count;

    do { /* ... while (nlevels > nlevels_old) */

      // find level structure for current root node

      count = cluster_start;

      for (i=0; i<_n; i++) {
	marks[i] = glob_marks[i];
	deg[i] = _degree[i];
      }
      nlevels_old = nlevels;
      nlevels = 0;
      l_begin = count;
      l_end = count;
      levels[count++] = root;		// add root node
      marks[root] = 0;			// unmark root node
      kend = _neigh_start[root+1];	// decrement degree of neighbours
      for (k=_neigh_start[root]; k<kend; k++) {
	deg[_neighbours[k]] --;
      }

      do { /* while (count > l_end) */

	// add nodes for next level
	
	l_begin = l_end;
	l_end = count;
	nlevels ++;
	for (i=l_begin; i<l_end; i++) {
	  node = levels[i];
	  jend = _neigh_start[node+1];
	  for (j=_neigh_start[node]; j<jend; j++) {
	    node = _neighbours[j];
	    if (marks[node]) {
	      levels[count++] = node;	// add node
	      marks[node] = 0;		// unmark added node
	      kend = _neigh_start[node+1];// decrement degree of neighbours
	      for (k=_neigh_start[node]; k<kend; k++) {
		deg[_neighbours[k]] --;
	      }
	    }
	  }

	  // determine degrees for next level
	  // sort added nodes

	  if (count - l_end > 1) {
	    for (k=l_end; k<count; k++) {
	      node = levels[k];
	      deg2[node] = 0;
	      jend = _neigh_start[node+1];
	      for (j=_neigh_start[node]; j<jend; j++) {
		if (marks[_neighbours[j]]) {
		  deg2[node] += deg[_neighbours[j]];
		}
	      }
	    }

	    qsort( (void *)(levels + l_end),
		  count - l_end, sizeof(int),
		  (int (*)(const void *, const void *))_compare);
	  }
	}
      } while (count > l_end);

      if (nlevels_old > nlevels)
	fprintf(stderr, "Hqp_Structure::order_rcm: Levels decreased!\n");

      // choice of the node with smallest degree in the last level
      // as root for next iteration
      root = levels[l_begin];
      min_deg = _degree[root];
      for (i=l_begin+1; i<l_end; i++) {
	if (_degree[levels[i]] < min_deg) {
	  root = levels[i];
	  min_deg = _degree[root];
	}
      }
    } while (nlevels > nlevels_old);

    root = _n;
    for (i=0; i<_n; i++) {
      glob_marks[i] = marks[i];
      if (marks[i])
	root = i;
    }
  }
  px_inv(_order, _order);

  // _order[i] says, where to put row/col i into reordered matrix

  // revert _order
  for (i=0; i<_n; i++)
    _order->pe[i] = _n - _order->pe[i] - 1;

  delete [] marks;
  delete [] glob_marks;
  delete [] deg;
  delete [] deg2;
}

//------------------------------------------------------------------------
PERM *Hqp_Structure::px_get(PERM *px)
{
  px = px_copy(_order, px);
  return px;
}

//------------------------------------------------------------------------
const IVEC *Hqp_Structure::neighbours(int node)
{
  assert(node < _n);

  _neigh_out.ive = _neighbours + _neigh_start[node];
  _neigh_out.dim = _neigh_start[node+1] - _neigh_start[node];
  _neigh_out.max_dim = 0;

  return &_neigh_out;
}

//------------------------------------------------------------------------
const Hqp_BandSize *Hqp_Structure::bd_size()
{
  int i, i_idx;
  int j_idx, jend;
  int min_neigh, max_neigh;
  int lb, ub;
  int node, dist;

  // determine band size according current ordering

  lb = ub = 0;
  for (i_idx=0; i_idx<_n; i_idx++) {
    i = _order->pe[i_idx];
    min_neigh = max_neigh = i;
    jend = _neigh_start[i_idx+1];
    for (j_idx=_neigh_start[i_idx]; j_idx<jend; j_idx++) {
      node = _order->pe[_neighbours[j_idx]];
      if (node < min_neigh)
	min_neigh = node;
      else if (node > max_neigh)
	max_neigh = node;
    }

    dist = i - min_neigh;
    if (dist > lb)
      lb = dist;

    dist = max_neigh - i;
    if (dist > ub)
      ub = dist;
  }

  // fill and return _bandSize

  _bandSize.lb = lb;
  _bandSize.ub = ub;
  _bandSize.n  = _n;

  return &_bandSize;
}

/*
 * Hqp_IpRedSpBKP.C -- class definition
 *
 * rf, 9/13/94
 */

/*
    Copyright (C) 1994--2001  Ruediger Franke

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
#include <math.h>

extern "C" {
#include <sparse2.h>
}

#include <If_Int.h>
#include <If_Real.h>

#include "sprcm.h"
#include "Hqp_Program.h"
#include "Hqp_IpRedSpBKP.h"

IF_CLASS_DEFINE("RedSpBKP", Hqp_IpRedSpBKP, Hqp_IpMatrix);

int sp_size(const SPMAT *A)
{
  if (!A)
    return 0;
  int i, m = A->m;
  int size = sizeof(SPMAT);

  for (i=0; i<m; i++) {
    size += sizeof(SPROW);
    size += A->row[i].maxlen * sizeof(row_elt);
  }

  return size;
}
  

//--------------------------------------------------------------------------
Hqp_IpRedSpBKP::Hqp_IpRedSpBKP()
{
  _n = _me = _m = 0;
  _sbw = -1;
  _tol = 1.0;
  _CT = SMNULL;
  _J = SMNULL;
  _J_raw = SMNULL;
  _QP2J = PNULL;
  _J2QP = PNULL;
  _pivot = PNULL;
  _blocks = PNULL;
  _zw = VNULL;
  _scale = VNULL;
  _r12 = VNULL;
  _xy = VNULL;
  _CTC_degree = IVNULL;
  _CTC_neigh_start = IVNULL;
  _CTC_neighs = IVNULL;

  _ifList.append(new If_Int("mat_sbw", &_sbw));
  _ifList.append(new If_Real("mat_tol", &_tol));
}

//--------------------------------------------------------------------------
Hqp_IpRedSpBKP::~Hqp_IpRedSpBKP()
{
  sp_free(_CT);
  sp_free(_J);
  sp_free(_J_raw);
  px_free(_QP2J);
  px_free(_J2QP);
  px_free(_pivot);
  px_free(_blocks);
  v_free(_zw);
  v_free(_scale);
  v_free(_r12);
  v_free(_xy);
  iv_free(_CTC_degree);
  iv_free(_CTC_neigh_start);
  iv_free(_CTC_neighs);
}

//--------------------------------------------------------------------------
SPMAT *Hqp_IpRedSpBKP::sub_CTC(const PERM *px, SPMAT *Q)
// return Q - _CT * diag(_zw) * _CT'
// read:  _CT, _zw
// write: _scale
{
  int i, j, j_idx, j_end;
  int qi, qj, qj_idx;
  SPROW *crow, *qrow;
  Real sum, val;
  IVEC neigh_header;
  IVEC *neigh = &neigh_header;

  assert(Q->n == Q->m);
  assert((int)px->size == Q->m && Q->m >= _n);

  if (!Q->flag_diag)
    sp_diag_access(Q);

  neigh_header.max_dim = 0;

  crow = _CT->row;
  for (i=0; i<_n; i++, crow++) {

    qrow = Q->row + px->pe[i];
    if (crow->len <= 0) {
      val = qrow->elt[qrow->diag].val;
      _scale->ve[i] = min(1.0, sqrt(-1.0 / val));
    }
    else {

      // calculate diagonal entry
      sum = sprow_inprod(crow, _zw, crow);
      j_idx = qrow->diag;
      val = qrow->elt[j_idx].val -= sum;
      _scale->ve[i] = min(1.0, sqrt(-1.0 / val));

      // calculate resting entries
      neigh->ive = _CTC_neighs->ive + _CTC_neigh_start->ive[i];
      neigh->dim = _CTC_neigh_start->ive[i + 1] - _CTC_neigh_start->ive[i];
      j_end = neigh->dim;
      for (j_idx = 0; j_idx < j_end; j_idx++) {
	j = neigh->ive[j_idx];
	if (j < i) {
	  sum = sprow_inprod(crow, _zw, _CT->row + j);
	  qi = px->pe[i];
	  qj = px->pe[j];

	  // substract sum from Qij or Qji (entry from upper part only)

	  if (qi < qj) {
	    qrow = Q->row + qi;
	    qj_idx = sprow_idx(qrow, qj);
	    if (qj_idx < 0) {
	      // the structure must already have been allocated in init()
	      error(E_INTERN, "Hqp_IpRedSpBKP");
	    }
	    else {
	      qrow->elt[qj_idx].val -= sum;
	    }
	  }
	  else {
	    qrow = Q->row + qj;
	    qj_idx = sprow_idx(qrow, qi);
	    if (qj_idx < 0) {
	      // the structure must already have been allocated in init()
	      error(E_INTERN, "Hqp_IpRedSpBKP");
	    }
	    else {
	      qrow->elt[qj_idx].val -= sum;
	    }
	  }
	}
      }
    }
  }

  return Q;  
}

//--------------------------------------------------------------------------
void Hqp_IpRedSpBKP::init(const Hqp_Program *qp)
{
  IVEC *degree, *neigh_start, *neighs;
  SPMAT *QCTC;
  SPROW *r1, *r2;
  int i, j;
  int len, dim;
  Real sum;

  _n = qp->c->dim;
  _me = qp->b->dim;
  _m = qp->d->dim;
  dim = _n + _me;

  // reallocations

  _pivot = px_resize(_pivot, dim);
  _blocks = px_resize(_blocks, dim);
  _zw = v_resize(_zw, _m);
  _scale = v_resize(_scale, _n);
  _r12 = v_resize(_r12, dim);
  _xy = v_resize(_xy, dim);

  // store C' for further computations
  // analyze structure of C'*C

  _CT = sp_transp(qp->C, _CT);
  sp_ones(_CT);
  v_ones(_zw);
  QCTC = sp_get(_n, _n, 10);
  r1 = _CT->row;
  for (i=0; i<_n; i++, r1++) {
    r2 = r1;
    for (j=i; j<_n; j++, r2++) {
      sum = sprow_inprod(r1, _zw, r2);
      if (sum != 0.0) {
	sp_set_val(QCTC, i, j, sum);
	if (i != j)
	  sp_set_val(QCTC, j, i, sum);
      }
    }
  }
  _CTC_degree = iv_resize(_CTC_degree, _n);
  _CTC_neigh_start = iv_resize(_CTC_neigh_start, _n + 1);
  _CTC_neighs = sp_rcm_scan(QCTC, SMNULL, SMNULL,
			    _CTC_degree, _CTC_neigh_start, _CTC_neighs);

  // initialize structure of reduced qp

  QCTC = sp_add(qp->Q, QCTC, QCTC);

  // determine RCM ordering

  degree = iv_get(dim);
  neigh_start = iv_get(dim + 1);
  neighs = sp_rcm_scan(QCTC, qp->A, SMNULL, degree, neigh_start, IVNULL);

  _QP2J = sp_rcm_order(degree, neigh_start, neighs, _QP2J);
  _sbw = sp_rcm_sbw(neigh_start, neighs, _QP2J);
  _J2QP = px_inv(_QP2J, _J2QP);

  iv_free(degree);
  iv_free(neigh_start);
  iv_free(neighs);

  len = 1 + (int)(log((double)dim) / log(2.0));
  sp_free(_J);
  sp_free(_J_raw);
  _J_raw = sp_get(dim, dim, len);
  _J = SMNULL;

  // fill up data (to allocate _J_raw)
  sp_into_symsp(QCTC, -1.0, _J_raw, _QP2J, 0, 0);
  spT_into_symsp(qp->A, 1.0, _J_raw, _QP2J, 0, _n);
  sp_into_symsp(qp->A, 1.0, _J_raw, _QP2J, _n, 0);

  sp_free(QCTC);

  // prepare iterations

  update(qp);
}

//--------------------------------------------------------------------------
void Hqp_IpRedSpBKP::update(const Hqp_Program *qp)
{
  // prepare _J_raw
  sp_zero(_J_raw);
  symsp_into_symsp(qp->Q, -1.0, _J_raw, _QP2J, 0);
  spT_into_symsp(qp->A, 1.0, _J_raw, _QP2J, 0, _n);
  sp_into_symsp(qp->A, 1.0, _J_raw, _QP2J, _n, 0);

  // update _CT
  _CT = sp_transp(qp->C, _CT);
}

//--------------------------------------------------------------------------
void Hqp_IpRedSpBKP::factor(const Hqp_Program *qp,
			    const VEC *z, const VEC *w)
{
  assert((int)z->dim == _m && (int)w->dim == _m);

  int	i, i_end;
  int	j, j_end, k;
  Real	scale;
  SPROW *row;
  row_elt *elt;

  // copy _J_raw to _J
  _J = sp_copy3(_J_raw, _J);

  // augment _J
  v_slash(w, z, _zw);
  sub_CTC(_QP2J, _J);

  // diagonal scaling
  j_end = _n + _me;
  for (j=0; j<j_end; j++) {
    k = _QP2J->pe[j];
    row = &(_J->row[k]);
    elt = row->elt;
    scale = j<_n? _scale->ve[j]: 1.0;
    i_end = row->len;
    for (i=0; i<i_end; i++, elt++) {
      elt->val *= scale;
      k = _J2QP->pe[elt->col];
      if (k < _n)
	elt->val *= _scale->ve[k];
    }
  }

  // delete zeros
  sp_compact(_J, 0.0);

  // factorization of _J
  spBKPfactor(_J, _pivot, _tol);
}

//--------------------------------------------------------------------------
void Hqp_IpRedSpBKP::step(const Hqp_Program *qp, const VEC *z, const VEC *w,
			  const VEC *r1, const VEC *r2, const VEC *r3,
			  const VEC *r4, VEC *dx, VEC *dy, VEC *dz, VEC *dw)
{
  VEC v;

  assert((int)r1->dim == _n && (int)dx->dim == _n);
  assert((int)r2->dim == _me && (int)dy->dim == _me);
  assert((int)r3->dim == _m && (int)dz->dim == _m);
  assert((int)r4->dim == _m && (int)dw->dim == _m);

  // augment, copy, scale and permutate [r1;r2;r3] into _r12
  // calculate, permutate and scale x
  
  // augment r1
  // temporary store (W^{-1}r_4 + ZW^{-1}r_3) in dz
  v_part(_r12, 0, _n, &v);
  v_slash(w, r4, dw);
  v_star(_zw, r3, dz);
  v_add(dw, dz, dz);
  sp_mv_mlt(_CT, dz, &v);
  v_sub(r1, &v, &v);
  v_star(&v, _scale, &v);

  v_copy(r2, v_part(_r12, _n, _me, &v));

  px_vec(_J2QP, _r12, _r12);

  spBKPsolve(_J, _pivot, _r12, _xy);

  px_vec(_QP2J, _xy, _xy);

  v_star(v_part(_xy, 0, _n, &v), _scale, dx);
  v_copy(v_part(_xy, _n, _me, &v), dy);

  sp_vm_mlt(_CT, dx, dw);
  v_star(_zw, dw, dw);
  v_sub(dz, dw, dz);

  // calculate dw

  sv_mlt(-1.0, r3, dw);
  sp_mv_mltadd(dw, dx, qp->C, 1.0, dw);
  // usage of _CT is numerically worse!
  //sp_vm_mltadd(dw, dx, _CT, 1.0, dw);
}


//==========================================================================

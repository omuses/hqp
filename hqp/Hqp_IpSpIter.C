/*
 * Hqp_IpSpIter.C -- class definition
 *
 * rf, 9/13/94
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
#include <math.h>

extern "C" {
#include <sparse2.h>
#include <iter.h>
}

#include "Hqp_Program.h"
#include "Hqp_Structure.h"
#include "Hqp_IpSpIter.h"


//--------------------------------------------------------------------------
Hqp_IpSpIter::Hqp_IpSpIter()
{
  _n = _me = _m = 0;
  _CT = SMNULL;
  _J = SMNULL;
  _J_raw = SMNULL;
  _J_fct = SMNULL;
  _QP2J = PNULL;
  _J2QP = PNULL;
  _pivot = PNULL;
  _blocks = PNULL;
  _zw = VNULL;
  _scale = VNULL;
  _r12 = VNULL;
  _xy = VNULL;
  _test = VNULL;
  _CTC_structure = new Hqp_Structure;
  _iter = NULL;
}

//--------------------------------------------------------------------------
Hqp_IpSpIter::~Hqp_IpSpIter()
{
  sp_free(_CT);
  sp_free(_J);
  sp_free(_J_raw);
  sp_free(_J_fct);
  px_free(_QP2J);
  px_free(_J2QP);
  px_free(_pivot);
  px_free(_blocks);
  v_free(_zw);
  v_free(_scale);
  v_free(_r12);
  v_free(_xy);
  v_free(_test);
  delete _CTC_structure;
  iter_free(_iter);
}

//--------------------------------------------------------------------------
SPMAT *Hqp_IpSpIter::sub_CTC(const PERM *px, SPMAT *Q)
// return Q - _CT * diags(_zw) * _CT'
// read:  _CT, _zw
// write: _scale
{
  int i, j, j_idx, j_end;
  int qi, qj, qj_idx;
  SPROW *crow, *qrow;
  Real sum, val;
  const IVEC *neigh;

  assert(Q->n == Q->m);
  assert(px->size == Q->m && Q->m >= _n);

  if (!Q->flag_diag)
    sp_diag_access(Q);

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
      neigh = _CTC_structure->neighbours(i);
      j_end = neigh->dim;
      for (j_idx=0; j_idx<j_end; j_idx++) {
	j = neigh->ive[j_idx];
	if (j < i) {
	  sum = sprow_inprod(crow, _zw, _CT->row + j);
	  qi = px->pe[i];
	  qj = px->pe[j];

	  // substract sum from Qij and Qji

	  qrow = Q->row + qi;
	  qj_idx = sprow_idx(qrow, qj);
	  if (qj_idx < 0) {
	    error(E_INTERN, "Hqp_IpSpIter");
	  }
	  qrow->elt[qj_idx].val -= sum;

	  qrow = Q->row + qj;
	  qj_idx = sprow_idx(qrow, qi);
	  if (qj_idx < 0) {
	    error(E_INTERN, "Hqp_IpSpIter");
	  }
	  qrow->elt[qj_idx].val -= sum;
	}
      }
    }
  }

  return Q;  
}

//--------------------------------------------------------------------------
void Hqp_IpSpIter::init(const Hqp_Program *qp)
{
  Hqp_Structure *sp = new Hqp_Structure;
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
  _test = v_resize(_test, dim);
  
  iter_free(_iter);
  _iter = iter_get(dim, dim);
  _iter->info = NULL;
  v_rand(_xy);

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
  _CTC_structure->init_spmat(QCTC);

  // initialize structure of reduced qp

  QCTC = sp_add(qp->Q, QCTC, QCTC);
  sp->init_QA(QCTC, qp->A);
  sp->order_rcm();
  _QP2J = sp->px_get(_QP2J);
  _J2QP = px_inv(_QP2J, _J2QP);
  len = 2 * sp->nentries() / dim;
  sp_free(_J);
  sp_free(_J_raw);
  sp_free(_J_fct);
  _J_raw = sp_get(dim, dim, len);
  _J_fct = SMNULL;
  _J = SMNULL;

  // fill up data (to allocate _J_raw)
  sp_into_sp(QCTC, -1.0, _J_raw, _QP2J, 0, 0);
  spT_into_sp(qp->A, 1.0, _J_raw, _QP2J, 0, _n);
  sp_into_sp(qp->A, 1.0, _J_raw, _QP2J, _n, 0);

  delete sp;
  sp_free(QCTC);

  // prepare iterations

  update(qp);
}

//--------------------------------------------------------------------------
void Hqp_IpSpIter::update(const Hqp_Program *qp)
{
  // prepare _J_raw
  sp_zero(_J_raw);
  sp_into_sp(qp->Q, -1.0, _J_raw, _QP2J, 0, 0);
  spT_into_sp(qp->A, 1.0, _J_raw, _QP2J, 0, _n);
  sp_into_sp(qp->A, 1.0, _J_raw, _QP2J, _n, 0);

  // update _CT
  _CT = sp_transp(qp->C, _CT);

}

//--------------------------------------------------------------------------
void Hqp_IpSpIter::factor(const VEC *z, const VEC *w)
{
  assert(z->dim == _m && w->dim == _m);

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

  // incomplete Cholesky factorization
  _J_fct = sp_copy3(_J, _J_fct);
  spILUfactor(_J_fct, 0.1);
  px_ident(_pivot);
}

//--------------------------------------------------------------------------
Real Hqp_IpSpIter::solve(const VEC *r1, const VEC *r2, const VEC *r3,
			 VEC *dx, VEC *dy, VEC *dz)
{
  static VEC v;
  Real residuum = 0.0;

  assert(r1->dim == _n && dx->dim == _n);
  assert(r2->dim == _me && dy->dim == _me);
  assert(r3->dim == _m && dz->dim == _m);

  // augment, copy, scale and permutate [r1;r2;r3] into _r12
  // calculate, permutate and scale x
  
  // augment r1
  v_part(_r12, 0, _n, &v);
  v_star(_zw, r3, dz);
  sp_mv_mlt(_CT, dz, &v);
  v_sub(r1, &v, &v);
  v_star(&v, _scale, &v);

  v_copy(r2, v_part(_r12, _n, _me, &v));

  px_vec(_J2QP, _r12, _r12);

  iter_Ax(_iter, &sp_mv_mlt, _J);
  iter_ATx(_iter, &sp_vm_mlt, _J);
  iter_Bx(_iter, &iter_solve, this);
  _iter->b = v_copy(_r12, _iter->b);
  _iter->eps = 1e-9;
  _iter->limit = 10000;
  _iter->k = 20;
  v_rand(_xy);
  iter_mgcr(_iter);
  _xy = v_copy(_iter->x, _xy);
//  _xy = iter_splsqr(_J, _r12, 1e-9, _xy, 1000, &steps);
//if (_iter->steps >= _iter->limit)
  printf("%d\n", _iter->steps);

//  spBKPsolve(_J_fct, _pivot, _r12, _xy);
/*
  sp_mv_mlt(_J, _xy, _test);
  v_sub(_r12, _test, _test);
  residuum = v_norm2(_test);
*/
  px_vec(_QP2J, _xy, _xy);

  v_star(v_part(_xy, 0, _n, &v), _scale, dx);
  v_copy(v_part(_xy, _n, _me, &v), dy);

  sp_vm_mlt(_CT, dx, dz);
  v_sub(r3, dz, dz);
  v_star(_zw, dz, dz);

  return residuum;
}

VEC *Hqp_IpSpIter::iter_solve(void *cld, VEC *x, VEC *out)
{
  Hqp_IpSpIter *obj = (Hqp_IpSpIter *)cld;

  spLUsolve(obj->_J_fct, PNULL, x, out);
  
  return out;
}


//==========================================================================

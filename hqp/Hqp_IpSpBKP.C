/*
 * Hqp_IpSpBKP.C -- class definition
 *
 * rf, 10/12/95
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
}
#include "sprcm.h"

#include <If_Int.h>
#include <If_Float.h>

#include "Hqp_Program.h"
#include "Hqp_IpSpBKP.h"

IF_CLASS_DEFINE("SpBKP", Hqp_IpSpBKP, Hqp_IpMatrix);

//--------------------------------------------------------------------------
Hqp_IpSpBKP::Hqp_IpSpBKP()
{
  _n = _me = _m = 0;
  _sbw = -1;
  _tol = 1.0;
  _J = SMNULL;
  _J_raw = SMNULL;
  _QP2J = PNULL;
  _J2QP = PNULL;
  _pivot = PNULL;
  _scale = VNULL;
  _r123 = VNULL;
  _xyz = VNULL;

  _ifList.append(new If_Int("mat_sbw", &_sbw));
  _ifList.append(new If_Float("mat_tol", &_tol));
}

//--------------------------------------------------------------------------
Hqp_IpSpBKP::~Hqp_IpSpBKP()
{
  sp_free(_J);
  sp_free(_J_raw);
  px_free(_QP2J);
  px_free(_J2QP);
  px_free(_pivot);
  v_free(_scale);
  v_free(_r123);
  v_free(_xyz);
}

//--------------------------------------------------------------------------
void Hqp_IpSpBKP::init(const Hqp_Program *qp)
{
  int dim;
  IVEC *degree, *neigh_start, *neighs;

  _n = qp->c->dim;
  _me = qp->b->dim;
  _m = qp->d->dim;
  dim = _n + _me + _m;

  // determine RCM ordering

  degree = iv_get(dim);
  neigh_start = iv_get(dim + 1);
  neighs = sp_rcm_scan(qp->Q, qp->A, qp->C, degree, neigh_start, IVNULL);

  _QP2J = sp_rcm_order(degree, neigh_start, neighs, _QP2J);
  _sbw = sp_rcm_sbw(neigh_start, neighs, _QP2J);
  _J2QP = px_inv(_QP2J, _J2QP);

  iv_free(degree);
  iv_free(neigh_start);
  iv_free(neighs);

  // allocate memory

  sp_free(_J);
  sp_free(_J_raw);
  _J_raw = sp_get(dim, dim, 1 + (int)(log((double)dim) / log(2.0)));
  _J = SMNULL;
  _pivot = px_resize(_pivot, dim);
  _scale = v_resize(_scale, _m);
  _r123 = v_resize(_r123, dim);
  _xyz = v_resize(_xyz, dim);

  // fill up data

  update(qp);
}

//--------------------------------------------------------------------------
void Hqp_IpSpBKP::update(const Hqp_Program *qp)
{
  int i, i_dst;
  int dim = _n+_me+_m;

  sp_zero(_J_raw);
  symsp_into_symsp(qp->Q, -1.0, _J_raw, _QP2J, 0);
  spT_into_symsp(qp->A, 1.0, _J_raw, _QP2J, 0, _n);
  spT_into_symsp(qp->C, 1.0, _J_raw, _QP2J, 0, _n+_me);
  sp_into_symsp(qp->A, 1.0, _J_raw, _QP2J, _n, 0);
  sp_into_symsp(qp->C, 1.0, _J_raw, _QP2J, _n+_me, 0);

  for (i=_n+_me; i<dim; i++) {
    i_dst = _QP2J->pe[i];
    sp_set_val(_J_raw, i_dst, i_dst, 1.0);
  }

  // delete zeros
  sp_compact(_J_raw, 0.0);
}

//--------------------------------------------------------------------------
void Hqp_IpSpBKP::factor(const Hqp_Program *, const VEC *z, const VEC *w)
{
  assert((int)z->dim == _m && (int)w->dim == _m);

  int	i, i_end;
  int	j, j_end, k;
  int 	n_me = _n + _me;
  Real	wz;
  SPROW *row;
  row_elt *elt;

  // copy _J_raw to _J
  _J = sp_copy3(_J_raw, _J);

  // insert slacks
  j_end = _m;
  for (j = 0; j < j_end; j++) {
    wz = w->ve[j] / z->ve[j];
    k = _QP2J->pe[n_me + j];
    sp_set_val(_J, k, k, wz);	// insert diagonal entry
    _scale->ve[j] = min(1.0, sqrt(1.0 / wz));
  }

  // diagonal scaling
  j_end = _n + _me + _m;
  for (j=0; j<j_end; j++) {
    k = _QP2J->pe[j];
    row = &(_J->row[k]);
    elt = row->elt;
    wz = (j >= n_me)? _scale->ve[j - n_me]: 1.0;
    i_end = row->len;
    for (i=0; i<i_end; i++, elt++) {
      elt->val *= wz;
      k = _J2QP->pe[elt->col];
      if (k >= n_me)
	elt->val *= _scale->ve[k - n_me];
    }
  }

  // factorization of _J
  spBKPfactor(_J, _pivot, _tol);
}

//--------------------------------------------------------------------------
void Hqp_IpSpBKP::step(const Hqp_Program *qp, const VEC *z, const VEC *w,
		       const VEC *r1, const VEC *r2, const VEC *r3,
		       const VEC *r4, VEC *dx, VEC *dy, VEC *dz, VEC *dw)
{
  VEC v;

  assert((int)r1->dim == _n && (int)dx->dim == _n);
  assert((int)r2->dim == _me && (int)dy->dim == _me);
  assert((int)r3->dim == _m && (int)dz->dim == _m);
  assert((int)r4->dim == _m && (int)dw->dim == _m);

  // copy, scale and permutate [r1;r2;r3] into _r123
  // calculate, permutate and scale x
  
  v_copy(r1, v_part(_r123, 0, _n, &v));
  v_copy(r2, v_part(_r123, _n, _me, &v));
  v_part(_r123, _n+_me, _m, &v);
  v_slash(z, r4, &v);
  v_add(&v, r3,	&v);
  v_star(&v, _scale, &v);

  px_vec(_J2QP, _r123, _r123);

  spBKPsolve(_J, _pivot, _r123, _xyz);

  px_vec(_QP2J, _xyz, _xyz);

  v_copy(v_part(_xyz, 0, _n, &v), dx);
  v_copy(v_part(_xyz, _n, _me, &v), dy);
  v_star(v_part(_xyz, _n+_me, _m, &v), _scale, dz);

  // calculate dw

  sv_mlt(-1.0, r3, dw);
  sp_mv_mltadd(dw, dx, qp->C, 1.0, dw);
}


//==========================================================================

/*
 * Hqp_IpFullSpLU.C -- class definition
 *
 * rf, 6/2/94
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

#include "Hqp_Program.h"
#include "Hqp_Structure.h"
#include "Hqp_IpFullSpLU.h"


//--------------------------------------------------------------------------
Hqp_IpFullSpLU::Hqp_IpFullSpLU()
{
  _n = _me = _m = 0;
  _J = SMNULL;
  _J_raw = SMNULL;
  _J_fct = SMNULL;
  _QP2J = PNULL;
  _J2QP = PNULL;
  _pivot = PNULL;
  _scale = VNULL;
  _r123 = VNULL;
  _xyz = VNULL;
  _test = VNULL;
}

//--------------------------------------------------------------------------
Hqp_IpFullSpLU::~Hqp_IpFullSpLU()
{
  sp_free(_J);
  sp_free(_J_raw);
  sp_free(_J_fct);
  px_free(_QP2J);
  px_free(_J2QP);
  px_free(_pivot);
  v_free(_scale);
  v_free(_r123);
  v_free(_xyz);
  v_free(_test);
}

//--------------------------------------------------------------------------
void Hqp_IpFullSpLU::init(const Hqp_Program *qp)
{
  Hqp_Structure *sp = new Hqp_Structure;
  int len;
  int dim;

  _n = qp->c->dim;
  _me = qp->b->dim;
  _m = qp->d->dim;
  dim = _n + _me + _m;

  // allocate memory

  sp->init_QAC(qp->Q, qp->A, qp->C);
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
  _pivot = px_resize(_pivot, dim);
  _scale = v_resize(_scale, _m);
  _r123 = v_resize(_r123, dim);
  _xyz = v_resize(_xyz, dim);
  _test = v_resize(_test, dim);

  delete sp;

  // fill up data

  update(qp);
}

//--------------------------------------------------------------------------
void Hqp_IpFullSpLU::update(const Hqp_Program *qp)
{
  int i, i_dst;
  int dim = _n+_me+_m;

  sp_zero(_J_raw);
  sp_into_sp(qp->Q, -1.0, _J_raw, _QP2J, 0, 0);
  spT_into_sp(qp->A, 1.0, _J_raw, _QP2J, 0, _n);
  spT_into_sp(qp->C, 1.0, _J_raw, _QP2J, 0, _n+_me);
  sp_into_sp(qp->A, 1.0, _J_raw, _QP2J, _n, 0);
  sp_into_sp(qp->C, 1.0, _J_raw, _QP2J, _n+_me, 0);

  for (i=_n+_me; i<dim; i++) {
    i_dst = _QP2J->pe[i];
    sp_set_val(_J_raw, i_dst, i_dst, 0.0);
  }
}

//--------------------------------------------------------------------------
void Hqp_IpFullSpLU::factor(const VEC *z, const VEC *w)
{
  assert(z->dim == _m && w->dim == _m);

  int	i, i_end;
  int	j, j_end, k, j_idx;
  int 	dim = _n+_me+_m;
  Real	wz;
  SPROW *row;
  row_elt *elt;

  // copy _J_raw to _J
  _J = sp_copy3(_J_raw, _J);

  // insert slacks, diagonal scaling

  j_end = _m;
  for (j=0; j<j_end; j++) {
    wz = w->ve[j] / z->ve[j];
    k = _QP2J->pe[j+dim-_m];
    sp_set_val(_J, k, k, wz);	// insert diagonal entry
    wz = min(1.0, sqrt(1.0 / wz));
    _scale->ve[j] = wz;
    row = &(_J->row[k]);	// scale k'th row of _J
    elt = row->elt;
    i_end = row->len;
    for (i=0; i<i_end; i++, elt++) {
      elt->val *= wz;
    }
  }

  if (!_J->flag_col)
    sp_col_access(_J);
  for (j=0; j<j_end; j++) {	// scale columns of _J
    wz = _scale->ve[j];
    k = _QP2J->pe[j+dim-_m];
    i = _J->start_row[k];
    j_idx = _J->start_idx[k];
    while (i >= 0) {
      elt = &(_J->row[i].elt[j_idx]);
      elt->val *= wz;
      i = elt->nxt_row;
      j_idx = elt->nxt_idx;
    }
  }  

  // delete zeros
  sp_compact(_J, 0.0);

  // copy _J to _J_fct for factorization
  _J_fct = sp_copy3(_J, _J_fct);

  // factorization of _J_fct
//  spLUfactor2(_J_fct, _pivot);
  spLUfactor(_J_fct, _pivot, 0.1);
}

//--------------------------------------------------------------------------
Real Hqp_IpFullSpLU::solve(const VEC *r1, const VEC *r2, const VEC *r3,
			   VEC *dx, VEC *dy, VEC *dz)
{
  static VEC v;
  Real residuum;

  assert(r1->dim == _n && dx->dim == _n);
  assert(r2->dim == _me && dy->dim == _me);
  assert(r3->dim == _m && dz->dim == _m);

  // copy, scale and permutate [r1;r2;r3] into _r123
  // calculate, permutate and scale x
  
  v_copy(r1, v_part(_r123, 0, _n, &v));
  v_copy(r2, v_part(_r123, _n, _me, &v));
  v_star(r3, _scale, v_part(_r123, _n+_me, _m, &v));

  px_vec(_J2QP, _r123, _r123);

  spLUsolve(_J_fct, _pivot, _r123, _xyz);

  sp_mv_mlt(_J, _xyz, _test);
  v_sub(_r123, _test, _test);
  residuum = v_norm2(_test);

  px_vec(_QP2J, _xyz, _xyz);

  v_copy(v_part(_xyz, 0, _n, &v), dx);
  v_copy(v_part(_xyz, _n, _me, &v), dy);
  v_star(v_part(_xyz, _n+_me, _m, &v), _scale, dz);

  return residuum;
}


//==========================================================================

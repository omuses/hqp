/*
 * Hqp_IpFullBdLU.C -- class definition
 *
 * rf, 6/1/94
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
#include <matrix2.h>
}

#include "Hqp_Program.h"
#include "Hqp_Structure.h"
#include "Hqp_IpFullBdLU.h"

//--------------------------------------------------------------------------
Hqp_IpFullBdLU::Hqp_IpFullBdLU()
{
  _n = _me = _m = 0;
  _J = BDNULL;
  _J_raw = BDNULL;
  _J_fct = BDNULL;
  _QP2J = PNULL;
  _J2QP = PNULL;
  _pivot = PNULL;
  _scale = VNULL;
  _r123 = VNULL;
  _xyz = VNULL;
  _test = VNULL;
}

//--------------------------------------------------------------------------
Hqp_IpFullBdLU::~Hqp_IpFullBdLU()
{
  bd_free(_J);
  bd_free(_J_raw);
  bd_free(_J_fct);
  px_free(_QP2J);
  px_free(_J2QP);
  px_free(_pivot);
  v_free(_scale);
  v_free(_r123);
  v_free(_xyz);
  v_free(_test);
}

//--------------------------------------------------------------------------
void Hqp_IpFullBdLU::init(const Hqp_Program *qp)
{
  Hqp_Structure *sp = new Hqp_Structure;
  const Hqp_BandSize *bds;
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
  bds = sp->bd_size();
  _J = bd_resize(_J, bds->lb, bds->ub, dim);
  _J_raw = bd_resize(_J_raw, bds->lb, bds->ub, dim);
  _J_fct = bd_resize(_J_fct, bds->lb, bds->ub + bds->lb, dim);
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
void Hqp_IpFullBdLU::update(const Hqp_Program *qp)
{
  m_zero(_J_raw->mat);
  sp_into_bd(qp->Q, -1.0, _J_raw, _QP2J, 0, 0);
  spT_into_bd(qp->A, 1.0, _J_raw, _QP2J, 0, _n);
  spT_into_bd(qp->C, 1.0, _J_raw, _QP2J, 0, _n+_me);
  sp_into_bd(qp->A, 1.0, _J_raw, _QP2J, _n, 0);
  sp_into_bd(qp->C, 1.0, _J_raw, _QP2J, _n+_me, 0);
}

//--------------------------------------------------------------------------
void Hqp_IpFullBdLU::factor(const VEC *z, const VEC *w)
{
  assert(z->dim == _m && w->dim == _m);

  int	i, i_end;
  int	j, j_end, k, l;
  int	lb, ub, bw;
  Real	wz;
  Real	**Jmat;
  int	dim = _n+_me+_m;

  // copy _J_raw into _J

  bd_copy(_J_raw, _J);

  // insert slacks, diagonal scaling

  lb = _J->lb;
  ub = _J->ub;
  bw = _J->mat->m;
  Jmat = _J->mat->me;
  j_end = _m;
  for (j=0; j<j_end; j++) {
    wz = w->ve[j] / z->ve[j];
    k = _QP2J->pe[j+dim-_m];
    Jmat[lb][k] = wz;		// insert diagonal entry
    wz = min(1.0, sqrt(1.0 / wz));
    _scale->ve[j] = wz;
    i_end = _J_raw->mat->m;
    for (i=0; i<i_end; i++) {
      Jmat[i][k] *= wz;		// scale k'th col of _J
    }
    i = max(0, lb-k);
    l = max(0, k-lb);
    i_end = min(bw, bw + dim-k-1 - ub);
    for (; i<i_end; i++, l++) {
      Jmat[i][l] *= wz;		// scale k'th row of _J
    }
  }  

  // copy _J to _J_fct for factorization
  // (_J is hold for calculating residuum)

  bd_resize(_J_fct, _J->lb, _J->ub, dim);	// no memory reallocation!
  bd_copy(_J, _J_fct);

  // factorization of _J_fct

  bdLUfactor(_J_fct, _pivot);
}

//--------------------------------------------------------------------------
Real Hqp_IpFullBdLU::solve(const VEC *r1, const VEC *r2, const VEC *r3,
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

  bdLUsolve(_J_fct, _pivot, _r123, _xyz);

  bd_mv_mlt(_J, _xyz, _test);
  v_sub(_r123, _test, _test);
  residuum = v_norm2(_test);

  px_vec(_QP2J, _xyz, _xyz);

  v_copy(v_part(_xyz, 0, _n, &v), dx);
  v_copy(v_part(_xyz, _n, _me, &v), dy);
  v_star(v_part(_xyz, _n+_me, _m, &v), _scale, dz);

  return residuum;
}


//==========================================================================

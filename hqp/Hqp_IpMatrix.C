/*
 * Hqp_IpMatrix.C -- class definition
 *
 * rf, 9/14/96
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

#include "Hqp_IpMatrix.h"
#include "Hqp_Program.h"
#include <If_Float.h>

IF_BASE_DEFINE(Hqp_IpMatrix);

//--------------------------------------------------------------------------
Hqp_IpMatrix::Hqp_IpMatrix()
{
  _r1 = v_resize(v_get(1), 0);
  _r2 = v_resize(v_get(1), 0);
  _r3 = v_resize(v_get(1), 0);
  _r4 = v_resize(v_get(1), 0);

  _dx = v_resize(v_get(1), 0);
  _dy = v_resize(v_get(1), 0);
  _dz = v_resize(v_get(1), 0);
  _dw = v_resize(v_get(1), 0);

  _eps = 1e-10;

  _ifList.append(new If_Float("mat_eps", &_eps));
}

//--------------------------------------------------------------------------
Hqp_IpMatrix::~Hqp_IpMatrix()
{
  v_free(_r1);
  v_free(_r2);
  v_free(_r3);
  v_free(_r4);

  v_free(_dx);
  v_free(_dy);
  v_free(_dz);
  v_free(_dw);
}

//--------------------------------------------------------------------------
Real Hqp_IpMatrix::solve(const Hqp_Program *qp, const VEC *z, const VEC *w,
			 const VEC *r1, const VEC *r2, const VEC *r3,
			 const VEC *r4, VEC *dx, VEC *dy, VEC *dz, VEC *dw)
{
  int i;
  Real res, res_last;
  Real alpha;

  if (_dx->dim != dx->dim) {
    v_resize(_dx, dx->dim);
    v_resize(_r1, r1->dim);
  }
  if (_dy->dim != dy->dim) {
    v_resize(_dy, dy->dim);
    v_resize(_r2, r2->dim);
  }
  if (_dz->dim != dz->dim) {
    v_resize(_dz, dz->dim);
    v_resize(_dw, dw->dim);
    v_resize(_r3, r3->dim);
    v_resize(_r4, r4->dim);
  }

  // calculate the solution
  step(qp, z, w, r1, r2, r3, r4, dx, dy, dz, dw);
  res = residuum(qp, z, w, r1, r2, r3, r4, dx, dy, dz, dw);
  //fprintf(stderr, " %g", res);

  // iterative refinement
  for (i = 0; i < 5 && res > _eps; i++) {

    res_last = res;

    // recalculate the solution
    step(qp, z, w, _r1, _r2, _r3, _r4, _dx, _dy, _dz, _dw);

    alpha = 1.0;
    do {
      v_mltadd(dx, _dx, alpha, dx);
      v_mltadd(dy, _dy, alpha, dy);
      v_mltadd(dz, _dz, alpha, dz);
      v_mltadd(dw, _dw, alpha, dw);

      res = residuum(qp, z, w, r1, r2, r3, r4, dx, dy, dz, dw);
      //fprintf(stderr, " %g", res);

      if (res > res_last) {
	v_mltadd(dx, _dx, -alpha, dx);
	v_mltadd(dy, _dy, -alpha, dy);
	v_mltadd(dz, _dz, -alpha, dz);
	v_mltadd(dw, _dw, -alpha, dw);
	alpha -= 0.3;
      }
    }
    while (res > res_last && alpha > 0.0);

    if (alpha <= 0.0) {
      break;
    }
  }
  //fprintf(stderr, "\n");

  return res;
}

//--------------------------------------------------------------------------
Real Hqp_IpMatrix::residuum(const Hqp_Program *qp,
			    const VEC *z, const VEC *w,
			    const VEC *r1, const VEC *r2, const VEC *r3,
			    const VEC *r4,
			    VEC *dx, VEC *dy, VEC *dz, VEC *dw)
{
  Real res, res_part;

  if (_dx->dim != dx->dim) {
    v_resize(_dx, dx->dim);
    v_resize(_r1, r1->dim);
  }
  if (_dy->dim != dy->dim) {
    v_resize(_dy, dy->dim);
    v_resize(_r2, r2->dim);
  }
  if (_dz->dim != dz->dim) {
    v_resize(_dz, dz->dim);
    v_resize(_dw, dw->dim);
    v_resize(_r3, r3->dim);
    v_resize(_r4, r4->dim);
  }

  sp_mv_symmlt(qp->Q, dx, _r1);
  sp_vm_mltadd(_r1, dy, qp->A, -1.0, _r1);
  sp_vm_mltadd(_r1, dz, qp->C, -1.0, _r1);
  v_add(r1, _r1, _r1);

  sp_mv_mlt(qp->A, dx, _r2);
  v_sub(r2, _r2, _r2);

  v_add(v_star(z, dw, _r3), v_star(w, dz, _r4), _r4);
  v_sub(r4, _r4, _r4);

  sp_mv_mlt(qp->C, dx, _r3);
  v_sub(_r3, dw, _r3);
  v_sub(r3, _r3, _r3);

  res = v_norm_inf(_r1);
  res_part = v_norm_inf(_r2);
  res = max(res, res_part);
  res_part = v_norm_inf(_r3);
  res = max(res, res_part);
  res_part = v_norm_inf(_r4);
  res = max(res, res_part);

  return res;
}

//--------------------------------------------------------------------------
void Hqp_IpMatrix::step(const Hqp_Program *, const VEC *, const VEC *,
			const VEC *r1, const VEC *r2, const VEC *r3,
			const VEC *r4, VEC *dx, VEC *dy, VEC *dz, VEC *dw)
{
  v_copy(r1, dx);
  v_copy(r2, dy);
  v_copy(r3, dz);
  v_copy(r4, dw);
}


//==========================================================================

/*
 * Hqp_SqpPowell.C -- class definition
 *
 * rf, 6/8/94
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

#include <math.h>

#include <If_Int.h>
#include <If_Bool.h>
#include <If_Method.h>

#include "Hqp_SqpPowell.h"
#include "Hqp_Solver.h"
#include "Hqp_SqpProgram.h"
#include "Hqp_Program.h"

typedef If_Method<Hqp_SqpPowell> If_Cmd;

IF_CLASS_DEFINE("Powell", Hqp_SqpPowell, Hqp_SqpSolver);

//--------------------------------------------------------------------------
Hqp_SqpPowell::Hqp_SqpPowell()
{
  _relaxed = false;
  _watchdog_iter = -1;
  _watchdog_start = 10;
  _watchdog_credit = 0;
  _watchdog_logging = true;
  _damped_multipliers = false;
  _re = v_resize(v_get(1), 0);
  _r = v_resize(v_get(1), 0);
  _sy_y = v_resize(v_get(1), 0);
  _sz_z = v_resize(v_get(1), 0);
  _x0 = v_resize(v_get(1), 0);
  _y0 = v_resize(v_get(1), 0);
  _z0 = v_resize(v_get(1), 0);
  _xk = v_resize(v_get(1), 0);
  _xl = v_resize(v_get(1), 0);
  _qp_xl = v_resize(v_get(1), 0);
  _yl = v_resize(v_get(1), 0);
  _zl = v_resize(v_get(1), 0);

  _ifList.append(new If_Int("sqp_watchdog_start", &_watchdog_start));
  _ifList.append(new If_Int("sqp_watchdog_credit", &_watchdog_credit));
  _ifList.append(new If_Bool("sqp_watchdog_logging", &_watchdog_logging));
  _ifList.append(new If_Bool("sqp_damped_multipliers",
			     &_damped_multipliers));
  _ifList.append(new If_Cmd("sqp_init", &Hqp_SqpPowell::init, this));
}

//--------------------------------------------------------------------------
Hqp_SqpPowell::~Hqp_SqpPowell()
{
  v_free(_re);
  v_free(_r);
  v_free(_sy_y);
  v_free(_sz_z);
  v_free(_x0);
  v_free(_y0);
  v_free(_z0);
  v_free(_xk);
  v_free(_xl);
  v_free(_qp_xl);
  v_free(_yl);
  v_free(_zl);
}

//--------------------------------------------------------------------------
int Hqp_SqpPowell::init(IF_CMD_ARGS)
{
  int ret = Hqp_SqpSolver::init();

  Hqp_Program *qp = _prg->qp();

  _x0 = v_resize(_x0, qp->c->dim);
  _y0 = v_resize(_y0, qp->b->dim);
  _z0 = v_resize(_z0, qp->d->dim);
  _xk = v_resize(_xk, qp->c->dim);
  _xl = v_resize(_xl, qp->c->dim);
  _qp_xl = v_resize(_qp_xl, qp->c->dim);
  _yl = v_resize(_yl, qp->b->dim);
  _zl = v_resize(_zl, qp->d->dim);
  _re = v_resize(_re, qp->b->dim);
  _r = v_resize(_r, qp->d->dim);
  _sy_y = v_resize(_sy_y, qp->b->dim);
  _sz_z = v_resize(_sz_z, qp->d->dim);

  v_zero(_re);
  v_zero(_r);

  _relaxed = false;
  _watchdog_iter = -1;

  return ret;
}

#if 1
// Powell's algorithm
//--------------------------------------------------------------------------
VEC *Hqp_SqpPowell::update_r(const VEC *z, VEC *r)
{
  int i, m;
  const Real *ve1;
  Real *ve2;
  Real val1, val2;

  m = z->dim;

  if (_iter == 0) {
    if (!r || (int)r->dim != m)
      r = v_resize(r, m);
    ve1 = z->ve;
    ve2 = r->ve;
    for (i=0; i<m; i++, ve1++, ve2++) {
      *ve2 = fabs(*ve1);
    }
  }
  else {
    ve1 = z->ve;
    ve2 = r->ve;
    for (i=0; i<m; i++, ve1++, ve2++) {
      val1 = fabs(*ve1);
      val2 = *ve2;
      if (val1 > val2)
	*ve2 = val1;
      else
	*ve2 = 0.5 * (val1 + val2);
    }
  }

  return r;
}

#else
// Han's algorithm
//--------------------------------------------------------------------------
VEC *Hqp_SqpPowell::update_r(const VEC *z, VEC *r)
{
  int i, m;
  const Real *ve1;
  Real *ve2;
  Real val1, val2;

  m = z->dim;

  if (_iter == 0) {
    if (!r || (int)r->dim != m)
      r = v_resize(r, m);
    ve1 = z->ve;
    ve2 = r->ve;
    for (i=0; i<m; i++, ve1++, ve2++) {
      *ve2 = 2.0 * fabs(*ve1);
    }
  }
  else {
    ve1 = z->ve;
    ve2 = r->ve;
    for (i=0; i<m; i++, ve1++, ve2++) {
      val1 = fabs(*ve1);
      val2 = *ve2;
      if (val2 < 1.5 * val1)
	*ve2 = 2.0 * val1;
    }
  }

  return r;
}
#endif

//--------------------------------------------------------------------------
Real Hqp_SqpPowell::phi()
{
  int i, i_end;
  const Real *ve1, *ve2;
  Real ret;

  ret = _prg->f();

  ve1 = _prg->qp()->b->ve;
  ve2 = _re->ve;
  i_end = _re->dim;
  for (i=0; i<i_end; i++, ve1++, ve2++)
    ret += *ve2 * fabs(*ve1);

  ve1 = _prg->qp()->d->ve;
  ve2 = _r->ve;
  i_end = _r->dim;
  for (i=0; i<i_end; i++, ve1++, ve2++)
    ret -= *ve2 * min(0.0, *ve1);

  return ret;
}

//--------------------------------------------------------------------------
Real Hqp_SqpPowell::phi1()
{
  Hqp_Program *qp = _prg->qp();
  int i, i_end;
  const Real *ve1, *ve2;
  Real ret;
  VEC  *g1 = VNULL;

  ret = _prg->f();
  ret += in_prod(qp->c, qp->x);

  g1 = sp_mv_mlt(qp->A, qp->x, g1);
  v_add(g1, qp->b, g1);
  ve1 = g1->ve;
  ve2 = _re->ve;
  i_end = _re->dim;
  for (i=0; i<i_end; i++, ve1++, ve2++)
    ret += *ve2 * fabs(*ve1);

  g1 = v_resize(g1, qp->C->m);		// suspect in sp_mv_mlt of Mesch1.2b
  g1 = sp_mv_mlt(qp->C, qp->x, g1);
  v_add(g1, qp->d, g1);
  ve1 = g1->ve;
  ve2 = _r->ve;
  i_end = _r->dim;
  for (i=0; i<i_end; i++, ve1++, ve2++)
    ret -= *ve2 * min(0.0, *ve1);

  v_free(g1);

  return ret;
}

//--------------------------------------------------------------------------
void Hqp_SqpPowell::update_vals()
{
  Hqp_Program *qp = _prg->qp();
  Real dphi0, phi0, phik;
  Real n_alpha;

  // update penalty coeffizients

  if (_damped_multipliers) {
    v_copy(_y, _y0);
    v_copy(_z, _z0);
    v_sub(_solver->y(), _y, _sy_y);
    v_sub(_solver->z(), _z, _sz_z);
  }
  v_copy(_solver->y(), _y);
  v_copy(_solver->z(), _z);
  update_r(_solver->y(), _re);
  update_r(_solver->z(), _r);

  // calculate penalty function

  _x0 = v_copy(_prg->x(), _x0);
  phi0 = phik = phi();
  dphi0 = phi1() - phi0;

  if (dphi0 > 0.0) {
    //printf("No descending direction (%g).\n", dphi0);
    _alpha = _min_alpha;
  }
  else {
    _alpha = 1.0;
  }

  // watchdog
  if (_iter == 0)
    _phil = phi0; // initialize _phil in first iteration
  if (_watchdog_iter < 0) {
    _phil_test = _phil;	// phi0 from last iteration
    _phil = phi0;
  }
  if (_watchdog_credit > 0 && _iter >= _watchdog_start) {
    // Step 3
    if (phi0 <= _phil_test) {
      _relaxed = true;
      if (_watchdog_logging) {
        printf("r");
	fflush(stdout);
      }
      // Step 4 (same test as for Step 3)
      _watchdog_iter = _iter;
      _xl = v_copy(_x0, _xl);
      _qp_xl = v_copy(qp->x, _qp_xl);
      _yl = v_copy(_y, _yl);
      _zl = v_copy(_z, _zl);
      _phil = phi0;
      if (dphi0 < 0.0) // treat positive dphi0, e.g. near an optimizer
        _phil_test += 0.1 * _min_alpha * dphi0;
    }
    else {
      _relaxed = false;
      if (_watchdog_logging) {
        printf("s");
	fflush(stdout);
      }
    }
    // Step 5
    if (_watchdog_iter >= 0
	&& _iter >= _watchdog_iter + _watchdog_credit) {

      // backing store values from _watchdog_iter
      _x0 = v_copy(_xl, _x0);
      _prg->x(_x0);
      _y = v_copy(_yl, _y);
      _z = v_copy(_zl, _z);
      _prg->update(_y, _z);
      qp->x = v_copy(_qp_xl, qp->x);
      hela_restart();

      // update penalties
      if (_damped_multipliers) {
	v_copy(_y, _y0);
	v_copy(_z, _z0);
	v_zero(_sy_y);
	v_zero(_sz_z);
      }
      update_r(_y, _re);
      update_r(_z, _r);
      phi0 = phik = phi();
      dphi0 = phi1() - phi0;

      _phil = phi0;
      _relaxed = false;
      _watchdog_iter = -1;

      if (_watchdog_logging) {
	printf("b");
	fflush(stdout);
      }
    }
  }

  // obtain step size

  while (1) {
    _d = sv_mlt(_alpha, qp->x, _d);
    v_add(_x0, _d, _xk);
    if (_damped_multipliers && _alpha < 1.0) {
      v_mltadd(_y0, _sy_y, _alpha, _y);
      v_mltadd(_z0, _sz_z, _alpha, _z);
    }
    _prg->x(_xk);
    if (_alpha <= _min_alpha)
      break;
    if (_relaxed && _watchdog_credit > 0)
      // accept step length 1
      break;
    _prg->update_fbd();
    if (!is_finite(_prg->f())) {
      _alpha *= 0.1;
      continue;
    }
    phik = phi();
    if (phik <= (phi0 + 0.1 * _alpha * dphi0) || fabs(dphi0) <= _eps)
      break;
    n_alpha = 0.5*dphi0*_alpha*_alpha / (dphi0 * _alpha - (phik - phi0));
    if (fabs(_alpha - n_alpha) < _min_alpha)
      break;
    _alpha *= 0.1;
    _alpha = max(_alpha, n_alpha);
    _alpha = max(_alpha, _min_alpha);
  }

  _dphi = dphi0;
  _phi = phi0;
}

//--------------------------------------------------------------------------
Real Hqp_SqpPowell::val_L()
{
  Hqp_Program *qp = _prg->qp();
  Real ret;

  ret = _prg->f();
  ret -= in_prod(_y, qp->b);
  ret -= in_prod(_z, qp->d);

  return ret;
}

//==========================================================================

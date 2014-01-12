/*
 * Hqp_IpsFranke.C -- class definition
 *
 * rf, 5/28/94
 * 
 * References:
 *     Franke, R.:    Anwendung von Interior-Point-Methoden zur Loesung
 *                    zeitdiskreter Optimalsteuerungsprobleme.
 *                    Diploma thesis, TU Ilmenau, Germany, 1994.
 *     Wright, S.J.:  Interior point methods for optimal control of
 *                    discrete time systems JOTA 77(1):161-187, 1993.
 *
 * rf, 11/5/95: inserted _max_warm_iters
 * rf, 2/12/97: treat singular matrix errors
 * rf, 4/20/97: treat numerical overflow for _qp->x and _gap
 * rf, 8/13/98
 *   - rename Hqp_IpSolver to Hqp_IpsFranke
 *   - make Hqp_IpsFranke an exchangeable interface class
 */

/*
    Copyright (C) 1994--2014  Ruediger Franke

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

#include <If_Int.h>
#include <If_Real.h>
#include <If_Method.h>
#include <If_Module.h>

#include "Hqp_Program.h"
#include "Hqp_IpRedSpBKP.h"
#include "Hqp_IpsFranke.h"

IF_CLASS_DEFINE("Franke", Hqp_IpsFranke, Hqp_Solver);

typedef If_Method<Hqp_IpsFranke> If_Cmd;

//--------------------------------------------------------------------------
Hqp_IpsFranke::Hqp_IpsFranke()
{
  _n = _me = _m = 0;

  _w = VNULL;
  _r1 = VNULL;
  _r2 = VNULL;
  _r3 = VNULL;
  _r4 = VNULL;
  _a1 = VNULL;
  _a2 = VNULL;
  _a3 = VNULL;
  _dx = VNULL;
  _dy = VNULL;
  _dz = VNULL;
  _dw = VNULL;

  _matrix = new Hqp_IpRedSpBKP;

  _alpha = 1.0;
  _mu0 = 0.0;
  _beta = 0.995;
  _Ltilde = 0.0;
  _fail_iters = 0;
  _max_warm_iters = 15;

  _ifList.append(new If_Real("qp_gap", &_gap));
  _ifList.append(new If_Real("qp_alpha", &_alpha));
  _ifList.append(new If_Real("qp_beta", &_beta));
  _ifList.append(new If_Real("qp_rhomin", &_rhomin));
  _ifList.append(new If_Real("qp_mu0", &_mu0));
  _ifList.append(new If_Real("qp_Ltilde", &_Ltilde));
  _ifList.append(new If_Int("qp_fail_iters", &_fail_iters));
  _ifList.append(new If_Int("qp_max_warm_iters", &_max_warm_iters));
  _ifList.append(new IF_MODULE("qp_mat_solver", &_matrix, Hqp_IpMatrix));

  _ifList.append(new If_Cmd("qp_init", &Hqp_IpsFranke::init, this));
  _ifList.append(new If_Cmd("qp_update", &Hqp_IpsFranke::update, this));
  _ifList.append(new If_Cmd("qp_cold_start", &Hqp_IpsFranke::cold_start, this));
  _ifList.append(new If_Cmd("qp_hot_start", &Hqp_IpsFranke::hot_start, this));
  _ifList.append(new If_Cmd("qp_step", &Hqp_IpsFranke::step, this));
  _ifList.append(new If_Cmd("qp_solve", &Hqp_IpsFranke::solve, this));
}

//--------------------------------------------------------------------------
Hqp_IpsFranke::~Hqp_IpsFranke()
{
  v_free(_w);
  v_free(_r1);
  v_free(_r2);
  v_free(_r3);
  v_free(_r4);
  v_free(_a1);
  v_free(_a2);
  v_free(_a3);
  v_free(_dx);
  v_free(_dy);
  v_free(_dz);
  v_free(_dw);
  delete _matrix;
}

//--------------------------------------------------------------------------
int Hqp_IpsFranke::init(IF_CMD_ARGS)
{
  assert(_qp != NULL);

  // allocate matrices and vectors

  _n = _qp->Q->n;
  _me = _qp->A->m;
  _m = _qp->C->m;

  _y = v_resize(_y, _me);
  _z = v_resize(_z, _m);
  _w = v_resize(_w, _m);
  _r1 = v_resize(_r1, _n);
  _r2 = v_resize(_r2, _me);
  _r3 = v_resize(_r3, _m);
  _r4 = v_resize(_r4, _m);
  _a1 = v_resize(_a1, _n);
  _a2 = v_resize(_a2, _me);
  _a3 = v_resize(_a3, _m);
  _dx = v_resize(_dx, _n);
  _dy = v_resize(_dy, _me);
  _dz = v_resize(_dz, _m);
  _dw = v_resize(_dw, _m);

  // fill up internal data

  _matrix->init(_qp);

  return IF_OK;
}

//--------------------------------------------------------------------------
int Hqp_IpsFranke::update(IF_CMD_ARGS)
{
  _matrix->update(_qp);

  return IF_OK;
}

//--------------------------------------------------------------------------
int Hqp_IpsFranke::cold_start(IF_CMD_ARGS)
{
//  Real Ltilde;
  Real min_d;
  int  i;

  if (_m > 0) {
    _rhomin = 1000.0 * _m;
    min_d = v_min(_qp->d, NULL);

    if (_mu0 > 0) {
      // choose Ltilde according _mu0
      Real mean_d_h;
      mean_d_h = 0.5 * v_sum(_qp->d) / (Real) _m;
      _Ltilde = - mean_d_h + sqrt(mean_d_h*mean_d_h + (Real)_m*_rhomin*_mu0);
      _Ltilde = max(_Ltilde, -min_d);
    }
    else {
      // choose Ltilde according Wright
      // (rf, 10/29/95: let Ltilde be greater than 1e2*_m)
      Real norm_d;
      norm_d = v_norm_inf(_qp->d);
//      norm_d = max(1.0, norm_d);
      _Ltilde = max(norm_d, -min_d);	// actually useless with v_norm_inf()
      _Ltilde = max(_Ltilde, 1e2*_m);
    }
    _zeta = 1.0;

    v_zero(_qp->x);
    v_zero(_y);
    v_set(_z, _Ltilde / ((Real)_m * _m));
    v_set(_w, _Ltilde);
    v_add(_w, _qp->d, _w);
    for (i=0; i<_m; i++)
      _w->ve[i] += 1e-10;

    v_copy(_qp->c, _a1);
    sp_vm_mltadd(_a1, _z, _qp->C, -1.0, _a1);
    sv_mlt(1.0/_zeta, _a1, _a1);

    sv_mlt(-1.0/_zeta, _qp->b, _a2);

    v_set(_a3, _Ltilde);
    sv_mlt(1.0/_zeta, _a3, _a3);

    _gap = in_prod(_z, _w);
  }

  // initialize a program without inequaltiy constraints
  else {
    v_zero(_qp->x);
    v_zero(_y);
    v_copy(_qp->c, _a1);
    sv_mlt(-1.0, _qp->b, _a2);
    _zeta = 1.0;
    _gap = 0.0;
  }

  _iter = 0;
  _alpha = 1.0;

  _hot_started = 0;
  _result = Hqp_Infeasible;

  return IF_OK;
}

//--------------------------------------------------------------------------
// hot_start:
//  - just correct slack vectors a1, a2 and a3
//
int Hqp_IpsFranke::hot_start(IF_CMD_ARGS)
{
  int i;

  if (_m > 0) {
    _zeta = 1.0;

    for (i=0; i<_m; i++)
      _w->ve[i] += 1e-10;

    sp_mv_symmlt(_qp->Q, _qp->x, _a1);
    v_add(_a1, _qp->c, _a1);
    sp_vm_mltadd(_a1, _y, _qp->A, -1.0, _a1);
    sp_vm_mltadd(_a1, _z, _qp->C, -1.0, _a1);

    sp_mv_mlt(_qp->A, _qp->x, _a2);
    v_add(_a2, _qp->b, _a2);
    sv_mlt(-1.0, _a2, _a2);

    sp_mv_mlt(_qp->C, _qp->x, _a3);
    v_add(_a3, _qp->d, _a3);
    v_sub(_a3, _w, _a3);
    sv_mlt(-1.0, _a3, _a3);
    
    _gap = in_prod(_z, _w) + 1.0;
  }

  // initialize a program without inequaltiy constraints
  else {
    v_zero(_qp->x);
    v_zero(_y);
    v_copy(_qp->c, _a1);
    sv_mlt(-1.0, _qp->b, _a2);
    _zeta = 1.0;
    _gap = 0.0;
  }

  _iter = 0;
  _alpha = 1.0;

  _hot_started = 1;
  _result = Hqp_Infeasible;

  return IF_OK;
}

//--------------------------------------------------------------------------
int Hqp_IpsFranke::step(int, const char *[], const char **)
{
  int	i, i_end;
  Real 	mu;
  Real	val1, val2;
  Real	*ve1, *ve2;
  Real 	residuum;

  // determine mu according current duality _gap and _zeta
  if (_iter == 0)
    _alphabar = 1.0;

  val1 = _gap;
  if (1.0 / val1 < _rhomin || _alpha < 1.0) {
    mu = _alphabar * val1 / _rhomin;	   	// potential reduction
    mu += (1.0 - _alphabar) * val1 / (Real)_m;	// centering
  } else {
    mu = val1 * val1;				// quadratic converence
  }

  // prepare right hand side vectors

  sv_mlt(-_zeta, _a1, _r1);
  sv_mlt(-_zeta, _a2, _r2);
  sv_mlt(-_zeta, _a3, _r3);
  ve1 = _r4->ve;
  i_end = _m;
  for (i=0; i<i_end; i++)
    *ve1++ = _z->ve[i] * _w->ve[i] - mu;

  // solve the system

  m_catch(E_SING,
	  // try
	  _matrix->factor(_qp, _z, _w);
	  residuum = _matrix->solve(_qp, _z, _w,
				    _r1, _r2, _r3, _r4, _dx, _dy, _dz, _dw),
	  // catch(E_SING)
	  _result = Hqp_Degenerate;
	  return IF_OK);

  // step size determination (find maximal feasible step)

  val1 = 2.0;
  i_end = _m;
  ve1 = _dz->ve;
  ve2 = _z->ve;
  for (i = 0; i < i_end; i++) {
    if (*ve1 >= 0.0 && *ve2 < val1 * *ve1)
      val1 = *ve2 / *ve1;
    ve1++;
    ve2++;
  }
  ve1 = _dw->ve;
  ve2 = _w->ve;
  for (i = 0; i < i_end; i++) {
    if (*ve1 >= 0.0 && *ve2 < val1 * *ve1)
      val1 = *ve2 / *ve1;
    ve1++;
    ve2++;
  }
  
  val1 *= _beta;
  _alpha = min(1.0, val1);

  _alphabar = 0.5 * _alphabar + 0.5 * _alpha;

  if (_alphabar == 1.0) {
    _rhomin *= 2.0;
  } else if (_alphabar < 0.5 && _rhomin > 100.0 * _m) {
    _rhomin /= 2.0;
  }

  // perform step

  v_mltadd(_qp->x, _dx, -_alpha, _dx);
  v_mltadd(_y, _dy, -_alpha, _y);
  v_mltadd(_z, _dz, -_alpha, _z);
  v_mltadd(_w, _dw, -_alpha, _w);
  _zeta *= (1.0 - _alpha);

  _gap = in_prod(_z, _w);

  if (!is_finite(_gap) || !is_finite(v_norm_inf(_dx))) {
    _result = Hqp_Degenerate;
    return IF_OK;
  }

  v_copy(_dx, _qp->x); 
  _iter ++;

  // return result
  // (strange comparisons to filter out NaN)

  if (!(_zeta < _eps)) {
    if (_alpha < _eps)
      _result = Hqp_Suboptimal;
    else
      _result = Hqp_Infeasible;
  }
  //else if (!(_gap < _eps))
  else if (!(_gap < _eps) || !(residuum < _eps))
    _result = Hqp_Feasible;
  else {
    //fprintf(stderr, "Res %g\n", residuum);
    _result = Hqp_Optimal;
  }

  return IF_OK;
}

//--------------------------------------------------------------------------
int Hqp_IpsFranke::solve(int, const char *[], const char ** /*ret*/)
{
  Real 	gap1 = 0.0;

  _fail_iters = 0;

  do {
    do {
      step();
      if (_hot_started) {
	if (_iter == 1)
	  gap1 = _gap;
	else
	  if (_gap > gap1) {
	    //printf("Restarted cold after %d iters (%g)\n", _iter, _gap);
	    _fail_iters += _iter;
	    cold_start();
	  }
      }

      if (_iter + _fail_iters >= _max_iters) break;
      if (_hot_started && _iter >= _max_warm_iters) break;

    } while (_result != Hqp_Optimal
	     && _result != Hqp_Suboptimal && _result != Hqp_Degenerate);

    if (_hot_started && _result != Hqp_Optimal) {
      //fprintf(stderr, "Bad hot-start, lost %d iters\n", _iter);
      _fail_iters += _iter;
      cold_start();
    } else
      break;

  } while (1);

  _iter += _fail_iters;

  /*
  if (_result == Hqp_Insolvable) {
    if (ret)
      *ret = hqp_result_strings[_result];
    return IF_ERROR;
  }
  else
  */
  return IF_OK;
}

//==========================================================================

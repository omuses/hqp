/*
 * Hqp_SqpSolver.C -- class definition
 *
 * rf, 6/6/94
 *
 * rf, 8/13/98
 *   - make Hqp_Solver an exchangeable interface class
 *
 * rf, 12/23/99
 *   - extended member access methods
 *
 */

/*
    Copyright (C) 1994--2002  Ruediger Franke

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
#include <If_Real.h>
#include <If_RealVec.h>
#include <If_Method.h>
#include <If_Module.h>

#include "Hqp_SqpSolver.h"
#include "Hqp_Program.h"
#include "Hqp_SqpProgram.h"
#include "Hqp_IpsFranke.h"
#include "Hqp_Client.h"

#include "Hqp_HL_BFGS.h"

extern Hqp_SqpProgram *theSqpProgram;

typedef If_Method<Hqp_SqpSolver> If_Cmd;

IF_BASE_DEFINE(Hqp_SqpSolver);

//--------------------------------------------------------------------------
Hqp_SqpSolver::Hqp_SqpSolver()
{
  _prg = NULL;
  _x0 = v_resize(v_get(1), 0);
  _xk = v_resize(v_get(1), 0);
  _d = v_resize(v_get(1), 0);
  _y = v_resize(v_get(1), 0);
  _z = v_resize(v_get(1), 0);
  _dL = v_resize(v_get(1), 0);
  _grd_L = v_resize(v_get(1), 0);
  _solver = new Hqp_IpsFranke;
  _hela = new Hqp_HL_BFGS;
  _iter = 0;
  _max_iters = 500;
  _inf_iters = 0;
  _max_inf_iters = 10;
  _alpha = 1.0;
  _min_alpha = 1e-10;
  _phi = 0.0;
  _dphi = 0.0;
  _sQs = 0.0;
  _xQx = 0.0;
  _eps = 1e-6;
  _solver->eps(1e-10);
  _norm_dx = 0.0;
  _norm_x = 0.0;
  _norm_inf = 0.0;
  _norm_grd_L = 0.0;
  _norm_df = 0.0;
  _f_bak = 0.0;
  _logging = false;

  _hot_started = false;
  _qp_Q_hot = SMNULL;

  if (theSqpProgram)
    set_prg(theSqpProgram);

  _ifList.append(new If_Int("sqp_iter", &_iter));
  _ifList.append(new If_Int("sqp_max_iters", &_max_iters));
  _ifList.append(new If_Int("sqp_inf_iters", &_inf_iters));
  _ifList.append(new If_Int("sqp_max_inf_iters", &_max_inf_iters));
  _ifList.append(new If_Bool("sqp_logging", &_logging));
  _ifList.append(new If_Real("sqp_eps", &_eps));
  _ifList.append(new If_Real("sqp_alpha", &_alpha));
  _ifList.append(new If_Real("sqp_min_alpha", &_min_alpha));
  _ifList.append(new If_Real("sqp_norm_s", &_norm_dx));
  _ifList.append(new If_Real("sqp_norm_x", &_norm_x));
  _ifList.append(new If_Real("sqp_norm_inf", &_norm_inf));
  _ifList.append(new If_Real("sqp_norm_grd_L", &_norm_grd_L));
  _ifList.append(new If_Real("sqp_norm_df", &_norm_df));
  _ifList.append(new If_Real("sqp_dphi", &_dphi));
  _ifList.append(new If_Real("sqp_phi", &_phi));
  _ifList.append(new If_Real("sqp_sQs", &_sQs));
  _ifList.append(new If_Real("sqp_xQx", &_xQx));
  _ifList.append(new If_RealVec("sqp_y", &_y));
  _ifList.append(new If_RealVec("sqp_z", &_z));
  _ifList.append(new If_RealVec("sqp_grd_L", &_grd_L));

  _ifList.append(new IF_MODULE("sqp_qp_solver", &_solver, Hqp_Solver));
  _ifList.append(new IF_MODULE("sqp_hela", &_hela, Hqp_HL));

  _ifList.append(new If_Cmd("sqp_init", &Hqp_SqpSolver::init, this));
  _ifList.append(new If_Cmd("sqp_qp_update",
			    &Hqp_SqpSolver::qp_update, this));
  _ifList.append(new If_Cmd("sqp_qp_solve", &Hqp_SqpSolver::qp_solve, this));
  _ifList.append(new If_Cmd("sqp_step", &Hqp_SqpSolver::step, this));
  _ifList.append(new If_Cmd("sqp_solve", &Hqp_SqpSolver::solve, this));
  _ifList.append(new If_Cmd("sqp_hela_restart",
			    &Hqp_SqpSolver::hela_restart, this));
  _ifList.append(new If_Cmd("sqp_qp_reinit_bd",
			    &Hqp_SqpSolver::qp_reinit_bd, this));
}

//--------------------------------------------------------------------------
Hqp_SqpSolver::~Hqp_SqpSolver()
{
  delete _hela;
  delete _solver;
  if (_qp_Q_hot)
    sp_free(_qp_Q_hot);
  v_free(_x0);
  v_free(_xk);
  v_free(_d);
  v_free(_y);
  v_free(_z);
  v_free(_dL);
  v_free(_grd_L);
}

//--------------------------------------------------------------------------
void Hqp_SqpSolver::set_prg(Hqp_SqpProgram *prg)
{
  _prg = prg;
  _solver->qp(prg->qp());
}

//--------------------------------------------------------------------------
Real Hqp_SqpSolver::norm_inf(const Hqp_Program *qp)
{
  Real ret, tmp;

  if (qp->b->dim > 0)
    ret = v_norm_inf(qp->b);
  else
    ret = 0.0;

  if (qp->d->dim > 0) {
    tmp = -v_min(qp->d, NULL);
    ret = max(ret, tmp);
  }

  return ret;
}

//--------------------------------------------------------------------------
int Hqp_SqpSolver::init(IF_CMD_ARGS)
{
  if (!_prg) {
    m_error(E_NULL, "Hqp_SqpSolver::init");
  }

  Hqp_Program *qp = _prg->qp();

  if (!qp || !(VEC *)qp->c || !(VEC *)qp->b || !(VEC *)qp->d) {
    m_error(E_NULL, "Hqp_SqpSolver::init");
  }

  _hela->setup(_prg);

  _solver->qp(qp);
  _solver->init();

  _x0 = v_resize(_x0, qp->c->dim);
  _xk = v_resize(_xk, qp->c->dim);
  _d = v_resize(_d, qp->c->dim);
  _y = v_zero(v_resize(_y, qp->b->dim));
  _z = v_zero(v_resize(_z, qp->d->dim));
  _dL = v_zero(v_resize(_dL, qp->c->dim));
  _grd_L = v_zero(v_resize(_grd_L, qp->c->dim));

  _iter = 0;
  _inf_iters = 0;
  _alpha = 1.0;

  _hot_started = false;

  return IF_OK;
}

//--------------------------------------------------------------------------
int Hqp_SqpSolver::qp_update(int, char *[], char **)
{
  if (!_prg) {
    m_error(E_NULL, "Hqp_SqpSolver::qp_update");
  }

  Hqp_Program *qp = _prg->qp();

  if (!qp || !(VEC *)qp->c || !(VEC *)qp->b || !(VEC *)qp->d) {
    m_error(E_NULL, "Hqp_SqpSolver::qp_update");
  }

  if (_iter == 0) {
    _prg->update(_y, _z);

    // init Hessian
    _hela->init(_y, _z, _prg);

    // calculate x'Qx (use _xk just for temporal storage)
    _xk = sp_mv_symmlt(qp->Q, _prg->x(), _xk);
    _xQx = in_prod(_xk, _prg->x());
    _sQs = _xQx;

    _norm_inf = norm_inf(qp);
    _norm_df = 0.0;
    _norm_grd_L = v_norm_inf(qp->c);
    _norm_x = v_norm_inf(_prg->x());
    _dphi = 0.0;
    _phi = _prg->f();
  }
  else {
    /*
    // Hessian update
    if (_status != Hqp_Optimal || _alpha <= _min_alpha) {
      // reinit Hessian
      _prg->update(_y, _z);
      _hela->init(_prg);
      _grd_L = grd_L(qp, _grd_L);
    }
    else {
    */
      _dL = grd_L(qp, _dL);
      // update problem derivatives
      if (!_hot_started)
        _prg->update(_y, _z);
      // update gradient of Lagrangian
      _grd_L = grd_L(qp, _grd_L);
      _dL = v_sub(_grd_L, _dL, _dL);
      // update Hessian of Lagrangian
      _hela->update(_d, _dL, _alpha, _prg);

      // calculate x'Qx (use _xk just for temporal storage)
      _xk = sp_mv_symmlt(qp->Q, _prg->x(), _xk);
      _xQx = in_prod(_xk, _prg->x());
    /*
    }
    */
    _norm_inf = norm_inf(qp);
    _norm_df = fabs(_f_bak - _prg->f());
    _norm_grd_L = v_norm_inf(_grd_L);
  }

  return IF_OK;
}

//--------------------------------------------------------------------------
int Hqp_SqpSolver::qp_solve(int, char *[], char **)
{
  if (!_prg) {
    m_error(E_NULL, "Hqp_SqpSolver::qp_solve");
  }

  Hqp_Program *qp = _prg->qp();

  if (!qp || !(VEC *)qp->c || !(VEC *)qp->b || !(VEC *)qp->d) {
    m_error(E_NULL, "Hqp_SqpSolver::qp_update");
  }

  _f_bak = _prg->f();

  if (_iter == 0) {
    _solver->update();
    _solver->cold_start();
  }
  else {
    _solver->update();
    if (_status == Hqp_Optimal && _alpha > _min_alpha)
      _solver->hot_start();
    else
      _solver->cold_start();
  }

  _solver->solve();

  // calculate s'Qs and ||s|| (use _xk just for temporal storage)
  _xk = sp_mv_symmlt(qp->Q, qp->x, _xk);
  _sQs = in_prod(_xk, qp->x);
  _norm_dx = v_norm_inf(qp->x);

  return IF_OK;
}

//--------------------------------------------------------------------------
int Hqp_SqpSolver::hela_restart(int, char *[], char **)
{
  if (!_prg) {
    m_error(E_NULL, "Hqp_SqpSolver::hela_restart");
  }

  sp_zero(_prg->qp()->Q);
  _hela->init(_y, _z, _prg);

  //sp_ident(_prg->qp()->Q);

  //v_zero(_y);
  //v_zero(_z);

  return IF_OK;
}

//--------------------------------------------------------------------------
int Hqp_SqpSolver::qp_reinit_bd(int, char *[], char **)
{
  if (!_prg) {
    m_error(E_NULL, "Hqp_SqpSolver::qp_reinit_bd");
  }

  // reinit bounds in qp->b and qp->d
  _prg->reinit_bd();
  _norm_inf = norm_inf(_prg->qp());

  if (!_hot_started) {
    // keep final Hessian of last cold solution
    _qp_Q_hot = sp_copy3(_prg->qp()->Q, _qp_Q_hot);
    _hot_started = true;
  }
  else {
    // restore Hessian of last cold solution
    sp_copy3(_qp_Q_hot, _prg->qp()->Q);
  }

  return IF_OK;
}

//--------------------------------------------------------------------------
void Hqp_SqpSolver::feasible_vals()
{
  Hqp_Program *qp = _prg->qp();
  Real old_norm_inf = max(_norm_inf, _eps);

  v_zero(_y);
  v_zero(_z);

  _x0 = v_copy(_prg->x(), _x0);
  _alpha = 1.0;

  while (1) {
    _d = sv_mlt(_alpha, qp->x, _d);
    v_add(_x0, _d, _xk);
    _prg->x(_xk);
    _prg->update_fbd();
    _norm_inf = norm_inf(qp);

    if (is_finite(_prg->f()) && _norm_inf < 1e2 * old_norm_inf)
      break;

    _alpha *= 0.5;

    if (_alpha <= _min_alpha)
      break;
  }
}

//--------------------------------------------------------------------------
int Hqp_SqpSolver::step(int, char *[], char **)
{
  if (!_prg) {
    m_error(E_NULL, "Hqp_SqpSolver::step");
  }

  Hqp_Program *qp = _prg->qp();

  if (!qp || !(VEC *)qp->c || !(VEC *)qp->b || !(VEC *)qp->d) {
    m_error(E_NULL, "Hqp_SqpSolver::step");
  }

  _status = _solver->result();

  //if (_status == Hqp_Suboptimal || _status == Hqp_Degenerate)
  if (_status == Hqp_Suboptimal)
    feasible_vals();
  else {
    update_vals();
    if (_alpha <= _min_alpha)
      feasible_vals();
  }

  // calculate _norm_x
  _norm_x = v_norm_inf(_prg->x());
  _norm_inf = norm_inf(qp);

  _iter++;

  if (_status != Hqp_Optimal && _status != Hqp_Feasible)
    _inf_iters++;
  else
    _inf_iters = 0;

  return IF_OK;
}

//--------------------------------------------------------------------------
int Hqp_SqpSolver::solve(int, char *[], char **ret)
{
  init();

  do {
    qp_update();
    if (_iter > 0 && _norm_inf < _eps) {
      if (_iter >= _max_iters) break;
      if (_norm_df < _eps) break;
      if (_norm_dx < _eps) break;
      if (_norm_grd_L < _eps) break;
      //if (_status == Hqp_Insolvable) return IF_ERROR;
    }
    qp_solve();
    step();
    if (_alpha <= _min_alpha) break;
  } while (1);

  if (_status != Hqp_Optimal) {
    if (ret)
      *ret = hqp_result_strings[_status];
    return IF_ERROR;
  }
  else
    return IF_OK;
}

//--------------------------------------------------------------------------
VEC *Hqp_SqpSolver::grd_L(const Hqp_Program *qp, VEC *out)
{
  int n = qp->Q->m;

  if (!out || (int)out->dim != n) {
    out = v_resize(out, n);
  }

  // calculate gradient of Lagrangian
  v_zero(out);
  sp_vm_mltadd(out, _y, qp->A, 1.0, out);
  sp_vm_mltadd(out, _z, qp->C, 1.0, out);
  v_sub(qp->c, out, out);

  return out;
}

//==========================================================================

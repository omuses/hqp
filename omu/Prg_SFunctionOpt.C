/*
 * Prg_SFunctionOpt.C -- class definition
 *
 */

/*
    Copyright (C) 1997--2002  Ruediger Franke

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

#include "Prg_SFunctionOpt.h"

#include <stdlib.h>

#include <Hxi_MEX_SFunction.h>

#include <If_Int.h>
#include <If_Bool.h>
#include <If_Real.h>
#include <If_RealVec.h>
#include <If_RealMat.h>
#include <If_Int.h>
#include <If_IntVec.h>

// redefine assert to throw an error instead of aborting
#undef assert
#define assert(expr) if (!(expr)) error(E_INTERN, "assert(" #expr ")");

IF_CLASS_DEFINE("SFunctionOpt", Prg_SFunctionOpt, Omu_Program);

// define own Inf so that iftcl works
static const double INF = 8e88;

//--------------------------------------------------------------------------
Prg_SFunctionOpt::Prg_SFunctionOpt()
{
  _K = 1;
  _KK = 1;
  _sps = 1;

  _mdl_u_active = iv_get(_mdl_nu);
  _mdl_u_nominal = v_get(_mdl_nu);
  _mdl_u_min = v_get(_mdl_nu);
  _mdl_u_max = v_get(_mdl_nu);
  _mdl_u_ref = v_get(_mdl_nu);
  _mdl_u_weight = v_get(_mdl_nu);
  _mdl_der_u_min = v_get(_mdl_nu);
  _mdl_der_u_max = v_get(_mdl_nu);
  _mdl_der_u_weight = v_get(_mdl_nu);
  iv_zero(_mdl_u_active);
  v_set(_mdl_u_nominal, 1.0);
  v_set(_mdl_u_min, -INF);
  v_set(_mdl_u_max, INF);
  v_set(_mdl_u_ref, 0.0);
  v_set(_mdl_u_weight, 0.0);
  v_set(_mdl_der_u_min, -INF);
  v_set(_mdl_der_u_max, INF);
  v_set(_mdl_der_u_weight, 0.0);

  _mdl_y_active = iv_get(_mdl_ny);
  _mdl_y_nominal = v_get(_mdl_ny);
  _mdl_y_bias = v_get(_mdl_ny);
  _mdl_y_min = v_get(_mdl_ny);
  _mdl_y_max = v_get(_mdl_ny);
  _mdl_y_ref = v_get(_mdl_ny);
  _mdl_y_weight = v_get(_mdl_ny);
  _mdl_y_min_soft = v_get(_mdl_ny);
  _mdl_y_min_weight1 = v_get(_mdl_ny);
  _mdl_y_min_weight2 = v_get(_mdl_ny);
  _mdl_y_max_soft = v_get(_mdl_ny);
  _mdl_y_max_weight1 = v_get(_mdl_ny);
  _mdl_y_max_weight2 = v_get(_mdl_ny);
  iv_zero(_mdl_y_active);
  v_set(_mdl_y_nominal, 1.0);
  v_set(_mdl_y_bias, 0.0);
  v_set(_mdl_y_min, -INF);
  v_set(_mdl_y_max, INF);
  v_set(_mdl_y_ref, 0.0);
  v_set(_mdl_y_weight, 0.0);
  v_set(_mdl_y_min_soft, -INF);
  v_set(_mdl_y_min_weight1, 0.0);
  v_set(_mdl_y_min_weight2, 0.0);
  v_set(_mdl_y_max_soft, INF);
  v_set(_mdl_y_max_weight1, 0.0);
  v_set(_mdl_y_max_weight2, 0.0);

  // numbers of optimization variables
  _nu = 0;
  _nx = _nu + _mdl_nx;
  _nc = 0;
  _ns = 0;

  _mdl_us = m_get(_KK+1, _mdl_nu);
  _mdl_ys = m_get(_KK+1, _mdl_ny);

  _mdl_x_min = v_get(_mdl_nx);
  _mdl_x_nominal = v_get(_mdl_nx);
  _mdl_x_max = v_get(_mdl_nx);
  v_set(_mdl_x_nominal, 1.0);
  v_set(_mdl_x_min, -INF);
  v_set(_mdl_x_max, INF);

  _ifList.append(new If_Int("prg_sps", &_sps));

  _ifList.append(new If_IntVec("mdl_u_active", &_mdl_u_active));
  _ifList.append(new If_RealVec("mdl_u_nominal", &_mdl_u_nominal));
  _ifList.append(new If_RealVec("mdl_u_min", &_mdl_u_min));
  _ifList.append(new If_RealVec("mdl_u_max", &_mdl_u_max));
  _ifList.append(new If_RealVec("mdl_u_ref", &_mdl_u_ref));
  _ifList.append(new If_RealVec("mdl_u_weight", &_mdl_u_weight));
  _ifList.append(new If_RealVec("mdl_der_u_min", &_mdl_der_u_min));
  _ifList.append(new If_RealVec("mdl_der_u_max", &_mdl_der_u_max));
  _ifList.append(new If_RealVec("mdl_der_u_weight", &_mdl_der_u_weight));

  _ifList.append(new If_RealVec("mdl_y_nominal", &_mdl_y_nominal));
  _ifList.append(new If_RealVec("mdl_y_bias", &_mdl_y_bias));
  _ifList.append(new If_RealVec("mdl_y_min", &_mdl_y_min));
  _ifList.append(new If_RealVec("mdl_y_max", &_mdl_y_max));
  _ifList.append(new If_RealVec("mdl_y_ref", &_mdl_y_ref));
  _ifList.append(new If_RealVec("mdl_y_weight", &_mdl_y_weight));
  _ifList.append(new If_RealVec("mdl_y_min_soft", &_mdl_y_min_soft));
  _ifList.append(new If_RealVec("mdl_y_min_weight1", &_mdl_y_min_weight1));
  _ifList.append(new If_RealVec("mdl_y_min_weight2", &_mdl_y_min_weight2));
  _ifList.append(new If_RealVec("mdl_y_max_soft", &_mdl_y_max_soft));
  _ifList.append(new If_RealVec("mdl_y_max_weight1", &_mdl_y_max_weight1));
  _ifList.append(new If_RealVec("mdl_y_max_weight2", &_mdl_y_max_weight2));

  _ifList.append(new If_RealVec("mdl_x_nominal", &_mdl_x_nominal));
  _ifList.append(new If_RealVec("mdl_x_min", &_mdl_x_min));
  _ifList.append(new If_RealVec("mdl_x_max", &_mdl_x_max));

  _ifList.append(new If_RealMat("mdl_us", &_mdl_us));
  _ifList.append(new If_RealMat("mdl_ys", &_mdl_ys));
}

//--------------------------------------------------------------------------
Prg_SFunctionOpt::~Prg_SFunctionOpt()
{
  iv_free(_mdl_y_active);
  iv_free(_mdl_u_active);
  m_free(_mdl_ys);
  m_free(_mdl_us);
  v_free(_mdl_x_min);
  v_free(_mdl_x_max);
  v_free(_mdl_x_nominal);
  v_free(_mdl_y_max_weight2);
  v_free(_mdl_y_max_weight1);
  v_free(_mdl_y_max_soft);
  v_free(_mdl_y_min_weight2);
  v_free(_mdl_y_min_weight1);
  v_free(_mdl_y_min_soft);
  v_free(_mdl_y_ref);
  v_free(_mdl_y_weight);
  v_free(_mdl_y_min);
  v_free(_mdl_y_max);
  v_free(_mdl_y_bias);
  v_free(_mdl_y_nominal);
  v_free(_mdl_der_u_min);
  v_free(_mdl_der_u_max);
  v_free(_mdl_der_u_weight);
  v_free(_mdl_u_ref);
  v_free(_mdl_u_weight);
  v_free(_mdl_u_min);
  v_free(_mdl_u_max);
  v_free(_mdl_u_nominal);
}

//--------------------------------------------------------------------------
void Prg_SFunctionOpt::setup_stages(IVECP ks, VECP ts)
{
  int i;

  // setup optimization problem
  _K = _KK/_sps;
  stages_alloc(ks, ts, _K, _sps);

  // setup S-function
  setup_sfun();

  // check for optional S-function methods that are required
  assert(ssGetmdlDerivatives(_S) != NULL);

  // adapt sizes of model vectors
  iv_resize(_mdl_y_active, _mdl_ny);
  iv_zero(_mdl_y_active);

  v_resize(_mdl_x_nominal, _mdl_nx);
  v_resize(_mdl_x_min, _mdl_nx);
  v_resize(_mdl_x_max, _mdl_nx);
  v_set(_mdl_x_nominal, 1.0);
  v_set(_mdl_x_min, -INF);
  v_set(_mdl_x_max, INF);

  iv_resize(_mdl_u_active, _mdl_nu);
  v_resize(_mdl_u_nominal, _mdl_nu);
  v_resize(_mdl_u_min, _mdl_nu);
  v_resize(_mdl_u_max, _mdl_nu);
  v_resize(_mdl_u_ref, _mdl_nu);
  v_resize(_mdl_u_weight, _mdl_nu);
  v_resize(_mdl_der_u_min, _mdl_nu);
  v_resize(_mdl_der_u_max, _mdl_nu);
  v_resize(_mdl_der_u_weight, _mdl_nu);
  iv_zero(_mdl_u_active);
  v_set(_mdl_u_nominal, 1.0);
  v_set(_mdl_u_min, -INF);
  v_set(_mdl_u_max, INF);
  v_set(_mdl_u_ref, 0.0);
  v_set(_mdl_u_weight, 0.0);
  v_set(_mdl_der_u_min, -INF);
  v_set(_mdl_der_u_max, INF);
  v_set(_mdl_der_u_weight, 0.0);

  iv_resize(_mdl_y_active, _mdl_ny);
  v_resize(_mdl_y_nominal, _mdl_ny);
  v_resize(_mdl_y_bias, _mdl_ny);
  v_resize(_mdl_y_min, _mdl_ny);
  v_resize(_mdl_y_max, _mdl_ny);
  v_resize(_mdl_y_ref, _mdl_ny);
  v_resize(_mdl_y_weight, _mdl_ny);
  v_resize(_mdl_y_min_soft, _mdl_ny);
  v_resize(_mdl_y_min_weight1, _mdl_ny);
  v_resize(_mdl_y_min_weight2, _mdl_ny);
  v_resize(_mdl_y_max_soft, _mdl_ny);
  v_resize(_mdl_y_max_weight1, _mdl_ny);
  v_resize(_mdl_y_max_weight2, _mdl_ny);
  iv_zero(_mdl_y_active);
  v_set(_mdl_y_nominal, 1.0);
  v_set(_mdl_y_min, -INF);
  v_set(_mdl_y_bias, 0.0);
  v_set(_mdl_y_max, INF);
  v_set(_mdl_y_ref, 0.0);
  v_set(_mdl_y_weight, 0.0);
  v_set(_mdl_y_min_soft, -INF);
  v_set(_mdl_y_min_weight1, 0.0);
  v_set(_mdl_y_min_weight2, 0.0);
  v_set(_mdl_y_max_soft, INF);
  v_set(_mdl_y_max_weight1, 0.0);
  v_set(_mdl_y_max_weight2, 0.0);

  m_resize(_mdl_us, _KK+1, _mdl_nu);
  m_resize(_mdl_ys, _KK+1, _mdl_ny);
}

//--------------------------------------------------------------------------
void Prg_SFunctionOpt::setup(int k,
			     Omu_Vector &x, Omu_Vector &u, Omu_Vector &c)
{
  int i, idx, is, j;

  // complete general setup
  if (k == 0) {
    // obtain numbers of optimization variables (active model variables)
    _nu = 0;
    for (idx = 0; idx < _mdl_nu; idx++) {
      if (_mdl_u_active[idx])
	_nu++;
    }
    _nx = _nu + _mdl_nx;
    _nc = 0;
    _ns = 0;
    for (idx = 0; idx < _mdl_ny; idx++) {
      if (_mdl_y_min[idx] > -INF || _mdl_y_max[idx] < INF
	  || _mdl_y_weight[idx] != 0.0) {
	_mdl_y_active[idx] = 1;
	_nc++;
      }
      if (_mdl_y_min_soft[idx] > -INF) {
	//_mdl_y_active[idx] = 1;
	//_nc++;
	_ns++;
      }
      if (_mdl_y_max_soft[idx] < INF) {
	//_mdl_y_active[idx] = 1;
	//_nc++;
	_ns++;
      }
    }
  }

  // allocate states and controls for optimization
  int spsk;
  if (k < _K) {
    spsk = _sps;
    x.alloc(_nx);
    u.alloc(_nu + spsk*_ns);
  }
  else {
    spsk = 1;
    // at final time allocate additional final states for slack variables
    // as no control parameters allowed here
    x.alloc(_nx + spsk*_ns);
  }
  c.alloc(spsk*(_nc+_ns));

  // setup initial states
  if (k == 0) {
    for (i = 0, idx = 0; idx < _mdl_nu; idx++) {
      if (_mdl_u_active[idx]) {
	x.initial[i] = _mdl_us[0][idx] / _mdl_u_nominal[idx];
	x.min[i] = x.max[i] = x.initial[i];
	i++;
      }
    }
    for (i = _nu; i < _nx; i++) {
      x.initial[i] = _mdl_x0[i-_nu] / _mdl_x_nominal[i-_nu];
      x.min[i] = x.max[i] = x.initial[i];
    }
  }

  // setup control inputs
  for (i = 0, idx = 0; idx < _mdl_nu; idx++) {
    if (_mdl_u_active[idx]) {
      if (k > 0) {
	if (_mdl_u_min[idx] > -INF) {
	  x.min[i] = _mdl_u_min[idx] / _mdl_u_nominal[idx];
	}
	if (_mdl_u_max[idx] < INF) {
	  x.max[i] = _mdl_u_max[idx] / _mdl_u_nominal[idx];
	}
      }
      if (k < _K) {
	// rate of change bounds
	if (_mdl_der_u_min[idx] > -INF) {
	  u.min[i] = _mdl_der_u_min[idx] / _mdl_u_nominal[idx];
	}
	if (_mdl_der_u_max[idx] < INF) {
	  u.max[i] = _mdl_der_u_max[idx] / _mdl_u_nominal[idx];
	}
	// initial guess
	u.initial[i] =
	  (_mdl_us[ks(k+1)][idx] - _mdl_us[ks(k)][idx])
	  / (ts(ks(k+1)) - ts(ks(k)))
	  / _mdl_u_nominal[idx];
      }
      i++;
    }
  }

  // setup constraints on model outputs
  for (i = 0, is = 0, idx = 0; idx < _mdl_ny; idx++) {
    if (_mdl_y_min[idx] > -INF)
      for (j = 0; j < spsk; j++)
	c.min[i+j] = _mdl_y_min[idx] / _mdl_y_nominal[idx];
    if (_mdl_y_max[idx] < INF)
      for (j = 0; j < spsk; j++)
	c.max[i+j] = _mdl_y_max[idx] / _mdl_y_nominal[idx];
    if (_mdl_y_active[idx])
      i += spsk;
    if (_mdl_y_min_soft[idx] > -INF) {
      for (j = 0; j < spsk; j++)
	c.min[spsk*_nc + is + j] = _mdl_y_min_soft[idx] / _mdl_y_nominal[idx];
      is += spsk;
    }
    if (_mdl_y_max_soft[idx] < INF) {
      for (j = 0; j < spsk; j++)
	c.max[spsk*_nc + is + j] = _mdl_y_max_soft[idx] / _mdl_y_nominal[idx];
      is += spsk;
    }
  }

  // setup slack variables for soft constraints
  if (k < _K)
    for (i = _nu; i < _nu + spsk*_ns; i++) {
      u.min[i] = 0.0;
      u.initial[i] = 0.0;
    }
  else
    for (i = _nx; i < _nx + spsk*_ns; i++) {
      x.min[i] = 0.0;
      x.initial[i] = 0.0;
    }
}

//--------------------------------------------------------------------------
void Prg_SFunctionOpt::setup_struct(int k,
				    const Omu_Vector &x, const Omu_Vector &u,
				    Omu_DependentVec &xt, Omu_DependentVec &F,
				    Omu_DependentVec &f,
				    Omu_Dependent &f0, Omu_DependentVec &c)
{
  // consistic just takes states from optimizer
  m_ident(xt.Jx);
  m_zero(xt.Ju);
  xt.set_linear();

  // explicit ODE for continuous-time equations
  m_ident(F.Jxp);
  sm_mlt(-1.0, F.Jxp, F.Jxp);
  F.set_linear(Omu_Dependent::WRT_xp);
}

//--------------------------------------------------------------------------
void Prg_SFunctionOpt::init_simulation(int k,
				       Omu_Vector &x, Omu_Vector &u)
{
  // initial states
  if (k == 0)
    v_copy(x.initial, x);
}

//--------------------------------------------------------------------------
void Prg_SFunctionOpt::update(int kk,
			      const Omu_StateVec &x, const Omu_Vec &u,
			      const Omu_StateVec &xf,
			      Omu_DependentVec &f, Omu_Dependent &f0,
			      Omu_DependentVec &c)
{
  int i, idx, is, j;
  int spsk = kk < _KK? _sps: 1; // one sample at final time

  // junction conditions for continuous-time equations
  if (kk < _KK) {
    for (i = 0; i < _nx; i++)
      f[i] = xf[i];
  }

  // set simulation time
  ssSetT(_S, ts(kk));

  // initialize model inputs
  real_T **mdl_u = (real_T **)ssGetInputPortRealSignalPtrs(_S, 0);
  for (i = 0, idx = 0; idx < _mdl_nu; idx++) {
    if (_mdl_u_active[idx])
      _mdl_us[kk][idx] = x[i++] * _mdl_u_nominal[idx];
    *mdl_u[idx] = _mdl_us[kk][idx];
  }

  // initialize model (this might be required by the S-function
  // as time steps back, compared to the previous call from continuous)
  if (ssGetmdlInitializeConditions(_S) != NULL) {
    mdlInitializeConditions(_S);
    if (ssGetErrorStatus(_S)) {
      fprintf(stderr, "Error from mdlInitializeConditions: %s\n",
	      ssGetErrorStatus(_S));
      ssSetErrorStatus(_S, NULL);
      error(E_RANGE, "mdlInitializeConditions");
    }
  }

  // pass current states to model
  real_T *mdl_x = ssGetContStates(_S);
  for (i = 0; i < _mdl_nx; i++)
    mdl_x[i] = x[_nu + i] * _mdl_x_nominal[i];

  // call mdlOutputs
  mdlOutputs(_S, 0); 
  if (ssGetErrorStatus(_S)) {
    fprintf(stderr, "Error from mdlOutputs: %s\n", ssGetErrorStatus(_S));
    ssSetErrorStatus(_S, NULL);
    error(E_RANGE, "mdlOutputs");
  }

  // obtain model outputs
  real_T *mdl_y = ssGetOutputPortRealSignal(_S, 0);

  // correct model outputs with bias
  for (idx = 0; idx < _mdl_ny; idx++)
    mdl_y[idx] += _mdl_y_bias[idx];

  // utilize model outputs for constraints and optimization objective
  double help;
  f0 = 0.0;
  // contribution of controlled inputs
  for (i = 0, idx = 0; idx < _mdl_nu; idx++) {
    if (_mdl_u_active[idx]) {
      // control inputs
      help = (_mdl_us[kk][idx] - _mdl_u_ref[idx]) / _mdl_u_nominal[idx];
      f0 += _mdl_u_weight[idx]*help*help;
      // rates of change
      if (kk < _KK) {
	help = u[i];
	f0 += _mdl_der_u_weight[idx]*help*help;
      }
      i++;
    }
  }
  // contribution of active outputs
  for (i = 0, is = 0, idx = 0; idx < _mdl_ny; idx++) {
    // assign used outputs to constraints
    if (_mdl_y_active[idx]) {
      c[i + kk%spsk] = mdl_y[idx] / _mdl_y_nominal[idx];
      i += spsk;
    }
    // calculate objective term
    if (_mdl_y_weight[idx] != 0.0) {
      help = (mdl_y[idx] - _mdl_y_ref[idx]) / _mdl_y_nominal[idx];
      f0 += _mdl_y_weight[idx]*help*help;
    }
    // consider soft constraints
    if (_mdl_y_min_soft[idx] > -INF) {
      if (kk < _KK)
	help = u[_nu + is + kk%spsk];
      else
	help = x[_nx + is + kk%spsk];
      c[spsk*_nc + is + kk%spsk] = mdl_y[idx] / _mdl_y_nominal[idx] + help;
      f0 += _mdl_y_min_weight1[idx]*help + _mdl_y_min_weight2[idx]*help*help;
      is += spsk;
    }
    if (_mdl_y_max_soft[idx] < INF) {
      if (kk < _KK)
	help = u[_nu + is + kk%spsk];
      else
	help = x[_nx + is + kk%spsk];
      c[spsk*_nc + is + kk%spsk] = mdl_y[idx] / _mdl_y_nominal[idx] - help;
      f0 += _mdl_y_min_weight1[idx]*help + _mdl_y_min_weight2[idx]*help*help;
      is += spsk;
    }
  }
  //f0 *= x[0]; // extension for saving fuel
  // consider step length in objective (trapezoidal integration rule)
  if (kk == 0)
    f0 *= 0.5*(ts(kk+1) - ts(kk));
  else if (kk == _KK)
    f0 *= 0.5*(ts(kk) - ts(kk-1));
  else
    f0 *= 0.5*(ts(kk+1) - ts(kk-1));
#if 0
  // apply higher weight at final time 
  if (kk == _KK)
    f0 *= 100.0;
#endif
  // store values of model outputs
  for (i = 0; i < _mdl_ny; i++)
    _mdl_ys[kk][i] = value(mdl_y[i]);

  // obtain Jacobians if required
  if (f.is_required_J() || f0.is_required_g() || c.is_required_J()) {
    // call predefined update for numerical differentiation
    _S->mdlInfo->reservedForFutureInt[0] = 1;
    Omu_Program::update_grds(kk, x, u, xf, f, f0, c);
    _S->mdlInfo->reservedForFutureInt[0] = 0;
#if 1
    // recalculate gradient of f0 as complete numeric diff. gives worse results
    v_zero(f0.gx);
    v_zero(f0.gu);
    double dt;
    if (kk == 0)
      dt = 0.5*(ts(kk+1) - ts(kk));
    else if (kk == _KK)
      dt = 0.5*(ts(kk) - ts(kk-1));
    else
      dt = 0.5*(ts(kk+1) - ts(kk-1));
    // controlled inputs
    for (i = 0, idx = 0; idx < _mdl_nu; idx++) {
      // contribution of objective term
      if (_mdl_u_active[idx]) {
	// control inputs
	f0.gx[i] += dt * _mdl_u_weight[idx] *
	  2.0 * (x[i] - _mdl_u_ref[idx]/_mdl_u_nominal[idx]);
	// rates of change
	if (kk < _KK) {
	  f0.gu[i] += dt * _mdl_der_u_weight[idx] * 2.0 * u[i];
	}
	i++;
      }
    }
    // active outputs
    for (i = 0, is = 0, idx = 0; idx < _mdl_ny; idx++) {
      // contribution of objective term
      if (_mdl_y_weight[idx] != 0.0) {
	for (j = 0; j < _nx; j++)
	  f0.gx[j] += dt * _mdl_y_weight[idx] *
	    2.0 * (c[i + kk%spsk] - _mdl_y_ref[idx]/_mdl_y_nominal[idx]) *
	    c.Jx[i + kk%spsk][j];
      }
      if (_mdl_y_active[idx])
	i += spsk;

      // contribution of soft constraints
      if (_mdl_y_min_soft[idx] > -INF) {
	if (kk < _KK)
	  f0.gu[_nu + is + kk%spsk]
	    += dt*(_mdl_y_min_weight1[idx]
		   + 2.0*_mdl_y_min_weight2[idx]*u[_nu + is + kk%spsk]);
	else
	  f0.gx[_nx + is + kk%spsk]
	    += dt*(_mdl_y_min_weight1[idx]
		   + 2.0*_mdl_y_min_weight2[idx]*x[_nx + is + kk%spsk]);
	is += spsk;
      }
      if (_mdl_y_max_soft[idx] < INF) {
	if (kk < _KK)
	  f0.gu[_nu + is + kk%spsk]
	    += dt*(_mdl_y_max_weight1[idx]
		   + 2.0*_mdl_y_max_weight2[idx]*u[_nu + is + kk%spsk]);
	else
	  f0.gx[_nx + is + kk%spsk]
	    += dt*(_mdl_y_max_weight1[idx]
		   + 2.0*_mdl_y_max_weight2[idx]*x[_nx + is + kk%spsk]);
	is += spsk;
      }
    }
#endif
  }
}

//--------------------------------------------------------------------------
void Prg_SFunctionOpt::consistic(int kk, double t,
				 const Omu_StateVec &x, const Omu_Vec &u,
				 Omu_DependentVec &xt)
{
  int i, idx;

  // initialize model in first stage
  if (kk == 0 && ssGetmdlInitializeConditions(_S) != NULL) {

    // set simulation time
    ssSetT(_S, t);

    // initialize model inputs using linear interpolation over time
    double rt = (t - ts(kk)) / (ts(kk+1) - ts(kk));
    real_T **mdl_u = (real_T **)ssGetInputPortRealSignalPtrs(_S, 0);
    for (i = 0, idx = 0; idx < _mdl_nu; idx++) {
      if (_mdl_u_active[idx])
	*mdl_u[idx] = x[i++] * _mdl_u_nominal[idx];
      else
	*mdl_u[idx] = _mdl_us[kk][idx] * (1 - rt) + _mdl_us[kk+1][idx] * rt;
    }

    // initialize model
    mdlInitializeConditions(_S);
    if (ssGetErrorStatus(_S)) {
      fprintf(stderr, "Error from mdlInitializeConditions: %s\n",
	      ssGetErrorStatus(_S));
      ssSetErrorStatus(_S, NULL);
      error(E_RANGE, "mdlInitializeConditions");
    }
    // call mdlOutputs as mdlInitializeConditions might not really
    // perform the initialization
    mdlOutputs(_S, 0); 
    if (ssGetErrorStatus(_S)) {
      fprintf(stderr, "Error from mdlOutputs: %s\n", ssGetErrorStatus(_S));
      ssSetErrorStatus(_S, NULL);
      error(E_RANGE, "mdlOutputs");
    }
  }

  // take over states from optimizer
  v_copy(x, xt);
}

//--------------------------------------------------------------------------
void Prg_SFunctionOpt::continuous(int kk, double t,
				  const Omu_StateVec &x, const Omu_Vec &u,
				  const Omu_StateVec &xp, Omu_DependentVec &F)
{
  int i, idx;

  // set simulation time
  ssSetT(_S, t);

  // initialize model inputs using linear interpolation over time
  double rt = (t - ts(kk)) / (ts(kk+1) - ts(kk));
  real_T **mdl_u = (real_T **)ssGetInputPortRealSignalPtrs(_S, 0);
  for (i = 0, idx = 0; idx < _mdl_nu; idx++) {
    if (_mdl_u_active[idx])
      *mdl_u[idx] = x[i++] * _mdl_u_nominal[idx];
    else
      *mdl_u[idx] = _mdl_us[kk][idx] * (1 - rt) + _mdl_us[kk+1][idx] * rt;
  }

  // pass current states to model
  real_T *mdl_x = ssGetContStates(_S);
  for (i = 0; i < _mdl_nx; i++)
    mdl_x[i] = x[_nu + i] * _mdl_x_nominal[i];

  // evaluate continuous model equations
  mdlDerivatives(_S);
  if (ssGetErrorStatus(_S)) {
    fprintf(stderr, "Error from mdlDerivatives: %s\n", ssGetErrorStatus(_S));
    ssSetErrorStatus(_S, NULL);
    error(E_RANGE, "mdlDerivatives");
  }

  // get model derivatives and change to residual form
  real_T *mdl_xp = ssGetdX(_S);
  for (i = 0; i < _mdl_nx; i++)
    F[_nu + i] = mdl_xp[i]/_mdl_x_nominal[i] - xp[_nu + i];

  // model equations for piecewise linear control inputs
  for (i = 0; i < _nu; i++)
    F[i] = u[i] - xp[i];

  // obtain Jacobians if required
  if (F.is_required_J()) {
    // call predefined continuous for numerical differentiation
    _S->mdlInfo->reservedForFutureInt[0] = 1;
    Omu_Program::continuous_grds(kk, t, x, u, xp, F);
    _S->mdlInfo->reservedForFutureInt[0] = 0;
  }
}


//==========================================================================

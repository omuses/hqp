/*
 * Prg_SFunctionEst.C -- class definition
 *
 */

/*
    Copyright (C) 1997--2003  Ruediger Franke

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

#include "Prg_SFunctionEst.h"

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
#define assert(expr) if (!(expr)) m_error(E_INTERN, "assert(" #expr ")");

#define GET_SET_CB(vartype, prefix, name) \
  GET_CB(vartype, prefix, name), \
  IF_SET_CB(vartype, Prg_SFunctionEst, set_##name)

#define GET_CB(vartype, prefix, name) \
  prefix#name, \
  IF_GET_CB(vartype, Prg_SFunctionEst, name)

IF_CLASS_DEFINE("SFunctionEst", Prg_SFunctionEst, Omu_Program);

//--------------------------------------------------------------------------
Prg_SFunctionEst::Prg_SFunctionEst()
{
  _K = 0;
  _KK = 0;
  _multistage = true;

  _mdl_np = 0;

  _mdl_p_active = iv_get(_mdl_np);
  _mdl_x0_active = iv_get(_mdl_nx);
  _mdl_der_x0_min = v_get(_mdl_nx);
  _mdl_der_x0_max = v_get(_mdl_nx);
  _mdl_u_order = iv_get(_mdl_nu);
  _mdl_y_active = iv_get(_mdl_ny);
  _mdl_p_nominal = v_get(_mdl_np);
  _mdl_x_nominal = v_get(_mdl_nx);
  _mdl_y_nominal = v_get(_mdl_ny);
  iv_zero(_mdl_p_active);
  iv_zero(_mdl_x0_active);
  v_set(_mdl_der_x0_min, -Inf);
  v_set(_mdl_der_x0_max, Inf);
  iv_set(_mdl_u_order, 1);
  iv_zero(_mdl_y_active);
  v_ones(_mdl_p_nominal);
  v_ones(_mdl_x_nominal);
  v_ones(_mdl_y_nominal);

  _np = 0;
  _nx0 = 0;
  _ny = 0;
  _nx = _mdl_nx;

  _nex = 1;
  _mdl_x0s = m_get(_nex, _mdl_nx);
  _exs = iv_get(_KK+1);
  iv_zero(_exs);

  _mdl_us = m_get(_KK+1, _mdl_nu);
  _mdl_xs = m_get(_KK+1, _mdl_nx);
  _mdl_ys = m_get(_KK+1, _mdl_ny);
  _ys_ref = m_get(_KK+1, 1);

  _M = m_resize(m_get(1, 1), (_KK+1)*_ny, _np+_nx0);
  _dydx = m_resize(m_get(1, 1), _ny, _nx);
  _dxdpx0 = m_resize(m_get(1, 1), _nx, _np+_nx0);
  _dydpx0 = m_resize(m_get(1, 1), _ny, _np+_nx0);
  _dxfdx = m_resize(m_get(1, 1), _nx, _nx);
  _dxfdpx0 = m_resize(m_get(1, 1), _nx, _np+_nx0);

  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_p)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_p_min)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_p_max)));
  _ifList.append(new If_IntVec(GET_SET_CB(const IVECP, "", mdl_p_active)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_p_nominal)));

  // make mdl_x0 read only as mdl_x0s should be used
  // (mdl_x0s is needed to cover multiple experiments)
  _ifList.append(new If_RealVec(GET_CB(const VECP, "", mdl_x0)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_x0_min)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_x0_max)));
  _ifList.append(new If_IntVec(GET_SET_CB(const IVECP, "", mdl_x0_active)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_der_x0_min)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_der_x0_max)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_x_nominal)));

  _ifList.append(new If_IntVec(GET_SET_CB(const IVECP, "", mdl_u_order)));

  _ifList.append(new If_IntVec(GET_SET_CB(const IVECP, "", mdl_y_active)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_y_nominal)));

  _ifList.append(new If_RealMat(GET_SET_CB(const MATP, "", mdl_x0s)));
  _ifList.append(new If_RealMat(GET_SET_CB(const MATP, "", mdl_us)));
  _ifList.append(new If_RealMat(GET_SET_CB(const MATP, "", mdl_xs)));
  _ifList.append(new If_RealMat(GET_CB(const MATP, "", mdl_ys)));

  _ifList.append(new If_Int(GET_SET_CB(int, "prg_", nex)));
  _ifList.append(new If_Bool(GET_SET_CB(bool, "prg_", multistage)));
  _ifList.append(new If_RealMat(GET_SET_CB(const MATP, "prg_", ys_ref)));
  _ifList.append(new If_RealMat(GET_CB(const MATP, "prg_", M)));
}

//--------------------------------------------------------------------------
Prg_SFunctionEst::~Prg_SFunctionEst()
{
  m_free(_dxfdpx0);
  m_free(_dxfdx);
  m_free(_dydpx0);
  m_free(_dxdpx0);
  m_free(_dydx);
  m_free(_M);
  m_free(_ys_ref);
  m_free(_mdl_ys);
  m_free(_mdl_xs);
  m_free(_mdl_us);
  m_free(_mdl_x0s);
  iv_free(_exs);
  iv_free(_mdl_p_active);
  v_free(_mdl_der_x0_max);
  v_free(_mdl_der_x0_min);
  iv_free(_mdl_x0_active);
  iv_free(_mdl_u_order);
  iv_free(_mdl_y_active);
  v_free(_mdl_y_nominal);
  v_free(_mdl_x_nominal);
  v_free(_mdl_p_nominal);
}

//--------------------------------------------------------------------------
void Prg_SFunctionEst::read_mx_args(VECP p)
{
  mxArray *arg;
  int i, j, idx, nel;

  for (idx = 0, j = 0; j < _mdl_nargs; j++) {
    arg = _mx_args[j];
    if (mxIsDouble(arg)) {
      nel = mxGetNumberOfElements(arg);
      for (i = 0; i < nel; i++, idx++)
	p[idx] = mxGetPr(arg)[i];
    }
  }
  assert(idx == _mdl_np); // S-function parameters must not have changed
}

//--------------------------------------------------------------------------
void Prg_SFunctionEst::write_mx_args(VECP p)
{
  mxArray *arg;
  int i, j, idx, nel;

  for (idx = 0, j = 0; j < _mdl_nargs; j++) {
    arg = _mx_args[j];
    if (mxIsDouble(arg)) {
      nel = mxGetNumberOfElements(arg);
      for (i = 0; i < nel; i++, idx++)
	mxGetPr(arg)[i] = p[idx];
    }
  }
  assert(idx == _mdl_np); // S-function parameters must not have changed
}

//--------------------------------------------------------------------------
void Prg_SFunctionEst::write_active_mx_args(VECP p)
{
  mxArray *arg;
  int i, j, ip, idx, nel;

  for (ip = 0, idx = 0, j = 0; j < _mdl_nargs; j++) {
    arg = _mx_args[j];
    if (mxIsDouble(arg)) {
      nel = mxGetNumberOfElements(arg);
      for (i = 0; i < nel; i++, idx++)
	if (_mdl_p_active[idx])
	  mxGetPr(arg)[i] = p[ip++] * _mdl_p_nominal[idx];
    }
  }
  assert(idx == _mdl_np); // S-function parameters must not have changed
  assert(ip == _np);	  // _np must not have changed since setup_stages
}

//--------------------------------------------------------------------------
void Prg_SFunctionEst::setup_stages(IVECP ks, VECP ts)
{
  int i, j;

  // setup optimization problem
  if (_multistage) {
    _K = max(_K, _KK); // assume that one of both has been specified
    _KK = _K;
    stages_alloc(ks, ts, _K, 1);
  }
  else {
    if (_nex > 1) {
      m_error(E_FORMAT, "Prg_SFunctionEst::setup_stages: "
	      "prg_nex>1 requires prg_multistage=true");
    }
    _KK = max(_K, _KK); // assume that one of both has been specified
    _K = 1;
    stages_alloc(ks, ts, 1, _KK);
  }

  // setup S-function
  setup_sfun();

  // check for optional S-function methods that are required
  assert(ssGetmdlDerivatives(_S) != NULL);

  // determine number of parameters
  mxArray *arg;
  _mdl_np = 0;
  for (j = 0; j < _mdl_nargs; j++) {
    arg = _mx_args[j];
    // only consider parameters in double format for accessing via mxGetPr()
    if (mxIsDouble(arg))
      _mdl_np += mxGetNumberOfElements(arg);
  }
  
  // adapt sizes of model vectors
  _mdl_p.alloc(_mdl_np);
  v_zero(_mdl_p.min);

  _mdl_x0.alloc(_mdl_nx);
  v_copy(Prg_SFunction::_mdl_x0, _mdl_x0);

  iv_resize(_mdl_p_active, _mdl_np);
  iv_resize(_mdl_x0_active, _mdl_nx);
  v_resize(_mdl_der_x0_min, _mdl_nx);
  v_resize(_mdl_der_x0_max, _mdl_nx);
  iv_resize(_mdl_u_order, _mdl_nu);
  iv_resize(_mdl_y_active, _mdl_ny);
  v_resize(_mdl_p_nominal, _mdl_np);
  v_resize(_mdl_x_nominal, _mdl_nx);
  v_resize(_mdl_y_nominal, _mdl_ny);
  iv_zero(_mdl_p_active);
  iv_zero(_mdl_x0_active);
  v_set(_mdl_der_x0_min, -Inf);
  v_set(_mdl_der_x0_max, Inf);
  iv_set(_mdl_u_order, 1);
  iv_zero(_mdl_y_active);
  v_ones(_mdl_p_nominal);
  v_ones(_mdl_x_nominal);
  v_ones(_mdl_y_nominal);

  // initialize all _x0s with _x0
  m_resize(_mdl_x0s, _nex, _mdl_nx);
  for (i = 0; i < _nex; i++) {
    for (j = 0; j < _mdl_nx; j++)
      _mdl_x0s[i][j] = _mdl_x0[j];
  }
  iv_resize(_exs, _KK+1);
  iv_zero(_exs);

  m_resize(_mdl_us, _KK+1, _mdl_nu);
  m_resize(_mdl_xs, _KK+1, _mdl_nx);
  m_resize(_mdl_ys, _KK+1, _mdl_ny);

  // store parameters in _mdl_p
  read_mx_args(_mdl_p);
}

//--------------------------------------------------------------------------
void Prg_SFunctionEst::setup(int k,
			     Omu_VariableVec &x, Omu_VariableVec &u,
			     Omu_VariableVec &c)
{
  int i, idx, kk;

  // take over possibly modified _mdl_p
  write_mx_args(_mdl_p);

  // initialize number of experiment
  int ex = 0;
  bool new_experiment = k == 0 || ts(ks(k)) < ts(ks(k)-1);
  if (k > 0) {
    if (new_experiment)
      // new experiment as time steps back
      ex = _exs[ks(k)-1] + 1;
    else
      ex = _exs[ks(k)-1];
  }
  assert(ex < _nex);  // _nex and experiment data must fit
  // store number of experiment in all sample periods of stage
  _exs[ks(k)] = ex;
  if (k < _K) {
    for (kk = ks(k)+1; kk < ks(k+1); kk++)
      _exs[kk] = ex;
  }

  // complete general setup
  if (k == 0) {
    // obtain numbers of optimization variables (active model variables)
    _np = 0;
    for (i = 0; i < _mdl_np; i++) {
      if (_mdl_p_active[i])
	++_np;
    }
    _nx0 = 0;
    for (i = 0; i < _mdl_nx; i++) {
      if (_mdl_x0_active[i])
	++_nx0;
    }
    _ny = 0;
    for (i = 0; i < _mdl_ny; i++) {
      if (_mdl_y_active[i])
	++_ny;
    }
    _nx = _np + _mdl_nx;

    // adapt sizes of reference and nominal outputs
    m_resize(_ys_ref, _KK+1, _ny);

    m_resize(_M, (_KK+1)*_ny, _np+_nx0);
    m_resize(_dydx, _ny, _nx);
    m_resize(_dxdpx0, _nx, _np+_nx0);
    m_resize(_dydpx0, _ny, _np+_nx0);
    m_resize(_dxfdx, _nx, _nx);
    m_resize(_dxfdpx0, _nx, _np+_nx0);
  }

  if (!_multistage && _nx0 > 0)
    m_error(E_FORMAT, "Prg_SFunctionEst::setup: "
	    "estimation of initial states requires prg_multistage=true");

  // allocate optimization variables
  if (k > 0 && new_experiment) {
    x.alloc(_np, _nx);
    u.alloc(_mdl_nx);
  }
  else {
    if (_multistage || k == _K)
      x.alloc(_nx);
    else
      x.alloc(_np, _nx);
  }

  if (k == 0) {
    // setup parameters
    int ip = 0;
    for (i = 0; i < _mdl_np; i++) {
      if (_mdl_p_active[i]) {
	x.initial[ip] = _mdl_p[i] / _mdl_p_nominal[i];
	if (_mdl_p.min[i] > -Inf)
	  x.min[ip] = _mdl_p.min[i] / _mdl_p_nominal[i];
	if (_mdl_p.max[i] < Inf)
	  x.max[ip] = _mdl_p.max[i] / _mdl_p_nominal[i];
	++ip;
      }
    }
    assert(ip == _np);	// _np must not have changed since setup_stages

    // setup initial states
    for (i = _np; i < _nx; i++) {
      x.initial[i] = _mdl_x0s[0][i-_np] / _mdl_x_nominal[i-_np];
      if (!_mdl_x0_active[i-_np])
	x.min[i] = x.max[i] = x.initial[i];
      else {
	if (_mdl_x0.min[i-_np] > -Inf)
	  x.min[i] = _mdl_x0.min[i-_np] / _mdl_x_nominal[i-_np];
	if (_mdl_x0.max[i-_np] < Inf)
	  x.max[i] = _mdl_x0.max[i-_np] / _mdl_x_nominal[i-_np];
      }
    }
  }
  else if (new_experiment) {
    // setup initial states for subsequent experiments
    for (i = 0; i < _mdl_nx; i++) {
      u.initial[i] = _mdl_x0s[ex][i] / _mdl_x_nominal[i];
      if (!_mdl_x0_active[i])
	u.min[i] = u.max[i] = u.initial[i];
      else {
	if (_mdl_x0.min[i] > -Inf)
	  u.min[i] = _mdl_x0.min[i] / _mdl_x_nominal[i];
	if (_mdl_x0.max[i] < Inf)
	  u.max[i] = _mdl_x0.max[i] / _mdl_x_nominal[i];
      }
    }
  }
  // setup states of subsequent stages
  else {
    for (i = 0, idx = 0; idx < _mdl_np; idx++)
      if (_mdl_p_active[idx]) {
	x.initial[i] = _mdl_p[i] / _mdl_p_nominal[i];
	i++;
      }
    for (idx = 0; idx < _mdl_nx; idx++)
      x.initial[_np+idx] = _mdl_xs[ks(k)][idx] / _mdl_x_nominal[idx];
  }
  // constraints for assigning relevant model outputs to opt vars
  // and on time derivatives of initial states
  int nc = _ny;
  if (new_experiment) {
    for (i = 0; i < _mdl_nx; i++) {
      if (_mdl_der_x0_min[i] > -Inf || _mdl_der_x0_max[i] < Inf)
	++nc;
    }
  }
  c.alloc(nc);
  if (new_experiment) {
    for (i = _ny, idx = 0; idx < _mdl_nx; idx++) {
      if (_mdl_der_x0_min[idx] > -Inf)
	c.min[i] = _mdl_der_x0_min[idx] / _mdl_x_nominal[idx];
      if (_mdl_der_x0_max[idx] < Inf)
	c.max[i] = _mdl_der_x0_max[idx] / _mdl_x_nominal[idx];
      if (_mdl_der_x0_min[idx] > -Inf || _mdl_der_x0_max[idx] < Inf)
	++i;
    }
  }
}

//--------------------------------------------------------------------------
void Prg_SFunctionEst::setup_struct(int k,
				    const Omu_VariableVec &x,
				    const Omu_VariableVec &u,
				    Omu_DependentVec &xt, Omu_DependentVec &F,
				    Omu_DependentVec &f,
				    Omu_Dependent &f0, Omu_DependentVec &c)
{
  int i, j;
  int ex = _exs[ks(k)]; 	// experiment of current stage
  bool new_experiment = k == 0 || ex != _exs[ks(k)-1];

  if (k < _K) {
    int ex1 = _exs[ks(k+1)]; 	// experiment of subsequent stage

    // consistic just takes states from optimizer
    if (k == 0 || !new_experiment) {
      m_ident(xt.Jx);
      m_zero(xt.Ju);
    }
    else {
      m_zero(xt.Jx);
      for (i = 0; i < _np; i++)
	xt.Jx[i][i] = 1.0;
      m_zero(xt.Ju);
      for (i = _np; i < _nx; i++)
	xt.Ju[i][i-_np] = 1.0;
    }
    xt.set_linear();

    if (ex == ex1) {
      // explicit ODE for continuous-time equations
      for (i = 0; i < _np; i++) {
	for (j = 0; j < _nx; j++)
	  F.Jx[i][j] = 0.0;
      }

      m_zero(F.Jxp);
      for (i = _np; i < _nx; i++)
	F.Jxp[i][i] = -1.0;
      F.set_linear(Omu_Dependent::WRT_xp);
    }
    else {
      // no continuous-time equations between two experiments
      m_zero(F.Jx);
      m_zero(F.Jxp);
    }
    m_zero(F.Ju);

    // f in update depends linearly on x for discrete-time states (parameters)
    m_zero(f.Jx);
    for (i = 0; i < _np; i++)
      f.Jx[i][i] = 1.0;
    f.set_linear(Omu_Dependent::WRT_x);

    // f does not depend on u
    m_zero(f.Ju);

    // all results of update depend linearly (or not) on xf
    m_zero(f.Jxf);
    if (ex == ex1) {
      for (i = _np; i < _nx; i++)
	f.Jxf[i][i] = 1.0;
    }
    f.set_linear(Omu_Dependent::WRT_xf);
    v_zero(f0.gxf);
    f0.set_linear(Omu_Dependent::WRT_xf);
    m_zero(c.Jxf);
    c.set_linear(Omu_Dependent::WRT_xf);
  }
}

//--------------------------------------------------------------------------
void Prg_SFunctionEst::init_simulation(int k,
				       Omu_VariableVec &x, Omu_VariableVec &u)
{
  // initial states
  if (k == 0)
    v_copy(x.initial, x);
}

//--------------------------------------------------------------------------
void Prg_SFunctionEst::update(int kk,
			      const Omu_StateVec &x, const Omu_Vec &u,
			      const Omu_StateVec &xf,
			      Omu_DependentVec &f, Omu_Dependent &f0,
			      Omu_DependentVec &c)
{
  int i, j, idx;
  int ex = _exs[kk];
  bool new_experiment = kk == 0 || ex != _exs[kk-1];

  // set simulation time
  ssSetT(_S, ts(kk));

  // pass estimated parameters to model
  write_active_mx_args(x);

  // initialize model inputs
  real_T *mdl_u;
  if (ssGetInputPortRequiredContiguous(_S, 0))
    mdl_u = (real_T *)ssGetInputPortRealSignal(_S, 0);
  else
    mdl_u = (real_T *)*ssGetInputPortRealSignalPtrs(_S, 0);
  for (idx = 0; idx < _mdl_nu; idx++)
    mdl_u[idx] = _mdl_us[kk][idx];

  // initialize model (this is required here as parameters change)
  if (_np > 0 && ssGetmdlInitializeConditions(_S) != NULL) {
    mdlInitializeConditions(_S);
    if (ssGetErrorStatus(_S)) {
      fprintf(stderr, "Error from mdlInitializeConditions: %s\n",
	      ssGetErrorStatus(_S));
      ssSetErrorStatus(_S, NULL);
      m_error(E_RANGE, "mdlInitializeConditions");
    }
  }

  // pass current states to model
  real_T *mdl_x = ssGetContStates(_S);
  if (kk == 0 || !new_experiment) {
    for (idx = 0; idx < _mdl_nx; idx++)
      mdl_x[idx] = x[_np + idx] * _mdl_x_nominal[idx];
  }
  else {
    // control parameters contain initial states for new experiments
    for (idx = 0; idx < _mdl_nx; idx++)
      mdl_x[idx] = u[idx] * _mdl_x_nominal[idx];
  }

  // obtain model outputs
  mdlOutputs(_S, 0);
  if (ssGetErrorStatus(_S)) {
    fprintf(stderr, "Error from mdlOutputs: %s\n", ssGetErrorStatus(_S));
    ssSetErrorStatus(_S, NULL);
    m_error(E_RANGE, "mdlOutputs");
  }

  // store outputs in constraints
  real_T *mdl_y = ssGetOutputPortRealSignal(_S, 0);
  for (i = 0, idx = 0; idx < _mdl_ny; idx++) {
    if (_mdl_y_active[idx]) {
      c[i++] = mdl_y[idx] / _mdl_y_nominal[idx];
    }
  }
  assert(i == _ny);	// _ny must not have changed since setup_stages

  // calculate objective term
  f0 = 0.0;
  double help;
  for (i = 0, idx = 0; idx < _mdl_ny; idx++) {
    if (_mdl_y_active[idx]) {
      help = c[i] - _ys_ref[kk][i]/_mdl_y_nominal[idx];
      f0 += help*help;
      ++i;
    }
  }

  // constraints on time derivatives of initial states
  if (new_experiment) {
    real_T *mdl_xp = ssGetdX(_S);
    mdlDerivatives(_S);
    if (ssGetErrorStatus(_S)) {
      fprintf(stderr, "Error from mdlDerivatives: %s\n",
	      ssGetErrorStatus(_S));
      ssSetErrorStatus(_S, NULL);
      m_error(E_RANGE, "mdlDerivatives");
    }
    for (i = _ny, idx = 0; idx < _mdl_nx; idx++) {
      if (_mdl_der_x0_min[idx] > -Inf || _mdl_der_x0_max[idx] < Inf)
	c[i++] = mdl_xp[idx] / _mdl_x_nominal[idx];
    }
  }

  if (kk < _KK) {
    // discrete-time state equations for parameters
    for (i = 0; i < _np; i++)
      f[i] = x[i];

    // junction conditions for state equations
    if (ex == _exs[kk+1]) {
      for (i = _np; i < _nx; i++)
	f[i] = xf[i];
    }
  }

  // obtain Jacobians if required
  if (f.is_required_J() || f0.is_required_g() || c.is_required_J()) {

    // store current model parameters and states
    if (kk == 0) {
      // store model parameters
      read_mx_args(_mdl_p);
    }
    if (new_experiment) {
      // store model states
      real_T *mdl_x = ssGetContStates(_S);
      for (idx = 0; idx < _mdl_nx; idx++)
	_mdl_x0s[ex][idx] = mdl_x[idx];
    }

    // store values of model states
    for (idx = 0; idx < _mdl_nx; idx++)
      _mdl_xs[kk][idx] = value(mdl_x[idx]);

    // store values of model outputs
    for (idx = 0; idx < _mdl_ny; idx++)
      _mdl_ys[kk][idx] = value(mdl_y[idx]);

    // call predefined update for numerical differentiation
    _S->mdlInfo->reservedForFutureInt[0] = 1;
    Omu_Program::update_grds(kk, x, u, xf, f, f0, c);
    _S->mdlInfo->reservedForFutureInt[0] = 0;

    // correct gradient of f0 as complete nd gives bad results
    v_zero(f0.gx);
    for (i = 0, idx = 0; idx < _mdl_ny; idx++) {
      if (_mdl_y_active[idx]) {
	for (j = 0; j < _nx; j++)
	  f0.gx[j] += 
	    2.0 * (c[i] - _ys_ref[kk][i]/_mdl_y_nominal[idx]) * c.Jx[i][j];
	++i;
      }
    }

    //
    // compute Measurement matrix
    //

    // build dx/d(p,x0)
    if (kk == 0) {
      m_zero(_dxdpx0);
      for (i = 0; i < _np; i++)
	_dxdpx0[i][i] = 1.0;
      for (i = _np, idx = 0; idx < _mdl_nx; idx++)
	if (_mdl_x0_active[idx])
	  _dxdpx0[_np + idx][i++] = 1.0;
    }
    else {
      if (!new_experiment)
	// initial states are final states of previous sample period
	m_copy(_dxfdpx0, _dxdpx0);
      else {
	// re-initialize initial states for a new experiment
	for (i = _np, idx = 0; idx < _mdl_nx; idx++) {
	  for (j = 0; j < _np+_nx0; j++)
	    _dxdpx0[_np + idx][j] = 0.0;
	  if (_mdl_x0_active[idx])
	    _dxdpx0[_np + idx][i++] = 1.0;
	}
      }	
    }

    // build dy/d(p,x0)
    if (kk == 0 || !new_experiment)
      m_copy(c.Jx, _dydx);
    else {
      // collect dydx from c.Jx and c.Ju
      for (i = 0; i < _ny; i++) {
	for (j = 0; j < _np; j++)
	  _dydx[i][j] = c.Jx[i][j];
	for (j = _np; j < _nx; j++)
	  _dydx[i][j] = c.Ju[i][j-_np];
      }
    }
    m_mlt(_dydx, _dxdpx0, _dydpx0);

    // store results in M
    for (i = 0; i < _ny; i++)
      for (j = 0; j < _np + _nx0; j++)
	_M[kk*_ny + i][j] = _dydpx0[i][j];

    // consider contribution of continuous-time equations
    if (kk < _KK && ex == _exs[kk+1]) {
      if (kk == 0 || !new_experiment)
	m_copy(xf.Sx, _dxfdx);
      else {
	// collect dxfdx from xf.Sx and xf.Su
	for (i = 0; i < _nx; i++) {
	  for (j = 0; j < _np; j++)
	    _dxfdx[i][j] = xf.Sx[i][j];
	  for (j = _np; j < _nx; j++)
	    _dxfdx[i][j] = xf.Su[i][j-_np];
	}
      }
      if (_multistage)
	m_mlt(_dxfdx, _dxdpx0, _dxfdpx0);
      else
	// if not multistage, then dxf/dx == dxf/d(p,x0) as x == (p,x0)
	m_copy(_dxfdx, _dxfdpx0);
    }
  }
}

//--------------------------------------------------------------------------
void Prg_SFunctionEst::consistic(int kk, double t,
				 const Omu_StateVec &x, const Omu_Vec &u,
				 Omu_DependentVec &xt)
{
  int i;
  int ex = _exs[kk];
  bool new_experiment = kk == 0 || ex != _exs[kk-1];

  // initialize model in first stage
  if (kk == 0 && ssGetmdlInitializeConditions(_S) != NULL) {

    // set simulation time
    ssSetT(_S, t);

    // pass estimated parameters to model
    write_active_mx_args(x);

    // initialize model inputs
    real_T *mdl_u;
    if (ssGetInputPortRequiredContiguous(_S, 0))
      mdl_u = (real_T *)ssGetInputPortRealSignal(_S, 0);
    else
      mdl_u = (real_T *)*ssGetInputPortRealSignalPtrs(_S, 0);
    for (i = 0; i < _mdl_nu; i++)
      mdl_u[i] = _mdl_us[kk][i];

    // initialize model
    mdlInitializeConditions(_S);
    if (ssGetErrorStatus(_S)) {
      fprintf(stderr, "Error from mdlInitializeConditions: %s\n",
	      ssGetErrorStatus(_S));
      ssSetErrorStatus(_S, NULL);
      m_error(E_RANGE, "mdlInitializeConditions");
    }
  }

  // take over states from optimizer
  if (kk == 0 || !new_experiment)
    v_copy(x, xt);
  else {
    for (i = 0; i < _np; i++)
      xt[i] = x[i];
    for (i = _np; i < _nx; i++)
      xt[i] = u[i-_np];
  }
}

//--------------------------------------------------------------------------
void Prg_SFunctionEst::continuous(int kk, double t,
				  const Omu_StateVec &x, const Omu_Vec &u,
				  const Omu_StateVec &xp, Omu_DependentVec &F)
{
  int i;

  // set simulation time
  ssSetT(_S, t);

  // pass estimated parameters to model
  write_active_mx_args(x);

  // initialize model inputs using linear interpolation over time
  double rt = (t - ts(kk)) / (ts(kk+1) - ts(kk));
  real_T *mdl_u;
  if (ssGetInputPortRequiredContiguous(_S, 0))
    mdl_u = (real_T *)ssGetInputPortRealSignal(_S, 0);
  else
    mdl_u = (real_T *)*ssGetInputPortRealSignalPtrs(_S, 0);
  for (i = 0; i < _mdl_nu; i++) {
    if (_mdl_u_order[i] == 0)
      mdl_u[i] = _mdl_us[kk][i];
    else
      mdl_u[i] = _mdl_us[kk][i] * (1 - rt) + _mdl_us[kk+1][i] * rt;
  }

  if (_np > 0 && ssGetmdlInitializeConditions(_S) != NULL) {
    // initialize model (this is required as parameters change)
    mdlInitializeConditions(_S);
    if (ssGetErrorStatus(_S)) {
      fprintf(stderr, "Error from mdlInitializeConditions: %s\n",
	      ssGetErrorStatus(_S));
      ssSetErrorStatus(_S, NULL);
      m_error(E_RANGE, "mdlInitializeConditions");
    }
  }

  // pass current states to model
  real_T *mdl_x = ssGetContStates(_S);
  for (i = 0; i < _mdl_nx; i++)
    mdl_x[i] = x[_np + i] * _mdl_x_nominal[i];

  // set model outputs before calculating derivatives
  // (note: this is required for sub-blocks providing inputs other blocks)
  mdlOutputs(_S, 0); 
  if (ssGetErrorStatus(_S)) {
    fprintf(stderr, "Error from mdlOutputs: %s\n", ssGetErrorStatus(_S));
    ssSetErrorStatus(_S, NULL);
    m_error(E_RANGE, "mdlOutputs");
  }

  // evaluate continuous model equations
  mdlDerivatives(_S);
  if (ssGetErrorStatus(_S)) {
    fprintf(stderr, "Error from mdlDerivatives: %s\n", ssGetErrorStatus(_S));
    ssSetErrorStatus(_S, NULL);
    m_error(E_RANGE, "mdlDerivatives");
  }

  // get model derivatives and change to residual form
  real_T *mdl_xp = ssGetdX(_S);
  for (i = 0; i < _mdl_nx; i++)
    F[_np + i] = mdl_xp[i]/_mdl_x_nominal[i] - xp[_np + i];

  // obtain Jacobians if required
  if (F.is_required_J()) {
    // call predefined continuous for numerical differentiation
    _S->mdlInfo->reservedForFutureInt[0] = 1;
    Omu_Program::continuous_grds(kk, t, x, u, xp, F);
    _S->mdlInfo->reservedForFutureInt[0] = 0;
  }
}


//==========================================================================

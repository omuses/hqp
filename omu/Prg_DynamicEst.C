/*
 * Prg_DynamicEst.C -- class definition
 *
 */

/*
    Copyright (C) 1997--2017  Ruediger Franke

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

#include "Prg_DynamicEst.h"

#include <stdlib.h>

#include <If_Int.h>
#include <If_Real.h>
#include <If_RealVec.h>
#include <If_RealMat.h>
#include <If_Int.h>
#include <If_IntVec.h>
#include <If_Method.h>

#include <meschach/matrix2.h> // for confidence intervals

// redefine assert to throw an error instead of aborting
#undef assert
#define assert(expr) if (!(expr)) m_error(E_INTERN, "assert(" #expr ")");

#define GET_SET_CB(vartype, prefix, name) \
  GET_CB(vartype, prefix, name), \
  IF_SET_CB(vartype, Prg_DynamicEst, set_##name)

#define GET_CB(vartype, prefix, name) \
  prefix#name, \
  IF_GET_CB(vartype, Prg_DynamicEst, name)

// Throw E_CONV as errors occuring during solution are generally
// due to bad values and need to be treated (esp. during stepsize check).
#undef SMETHOD_ERROR
#define SMETHOD_ERROR E_CONV

typedef If_Method<Prg_DynamicEst> If_Cmd;

IF_CLASS_DEFINE("DynamicEst", Prg_DynamicEst, Omu_Program);
IF_CLASS_DEFINE("SFunctionEst", Prg_SFunctionEst, Omu_Program);

//--------------------------------------------------------------------------
Prg_DynamicEst::Prg_DynamicEst()
{
  _multistage = 1;

  _mdl_p_active = iv_get(_mdl_np);
  _mdl_p_confidence = v_get(_mdl_np);
  _mdl_x0_active = iv_get(_mdl_nx);
  _mdl_x0_confidence = v_get(_mdl_nx);
  _mdl_der_x0_min = v_get(_mdl_nx);
  _mdl_der_x0_max = v_get(_mdl_nx);
  _mdl_u_order = iv_get(_mdl_nu);
  _mdl_y_active = iv_get(_mdl_ny);
  _mdl_p_nominal = v_get(_mdl_np);
  iv_zero(_mdl_p_active);
  v_zero(_mdl_p_confidence);
  iv_zero(_mdl_x0_active);
  v_zero(_mdl_x0_confidence);
  v_set(_mdl_der_x0_min, -Inf);
  v_set(_mdl_der_x0_max, Inf);
  iv_set(_mdl_u_order, 1);
  iv_zero(_mdl_y_active);
  v_ones(_mdl_p_nominal);

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
  _ys_ref = m_get(_KK+1, _ny);

  _ssr = 0.0;
  _M2 = m_resize(m_get(1, 1), (_KK+1)*_ny, _np+_nx0);
  _P2 = m_resize(m_get(1, 1), _np+_nx0, _np+_nx0);
  _COV = m_resize(m_get(1, 1), _np+_nx0, _np+_nx0);
  _dydx = m_resize(m_get(1, 1), _ny, _nx);
  _dxdpx0 = m_resize(m_get(1, 1), _nx, _np+_nx0);
  _dydpx0 = m_resize(m_get(1, 1), _ny, _np+_nx0);
  _dxfdx = m_resize(m_get(1, 1), _nx, _nx);
  _dxfdpx0 = m_resize(m_get(1, 1), _nx, _np+_nx0);

  _ifList.append(new If_Cmd("prg_setup_model",
			    &Prg_DynamicEst::setup_model, this));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_p)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_p_min)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_p_max)));
  _ifList.append(new If_IntVec(GET_SET_CB(const IVECP, "", mdl_p_active)));
  _ifList.append(new If_RealVec(GET_CB(const VECP, "", mdl_p_confidence)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_p_nominal)));

  // make mdl_x0 read only as mdl_x0s should be used
  // (mdl_x0s is needed to cover multiple experiments)
  _ifList.append(new If_RealVec(GET_CB(const VECP, "", mdl_x0)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_x0_min)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_x0_max)));
  _ifList.append(new If_IntVec(GET_SET_CB(const IVECP, "", mdl_x0_active)));
  _ifList.append(new If_RealVec(GET_CB(const VECP, "", mdl_x0_confidence)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_der_x0_min)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_der_x0_max)));

  _ifList.append(new If_IntVec(GET_SET_CB(const IVECP, "", mdl_u_order)));

  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_x_min)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_x_max)));

  _ifList.append(new If_IntVec(GET_SET_CB(const IVECP, "", mdl_y_active)));

  _ifList.append(new If_RealMat(GET_SET_CB(const MATP, "", mdl_x0s)));
  _ifList.append(new If_RealMat(GET_SET_CB(const MATP, "", mdl_us)));
  _ifList.append(new If_RealMat(GET_SET_CB(const MATP, "", mdl_xs)));
  _ifList.append(new If_RealMat(GET_CB(const MATP, "", mdl_ys)));

  _ifList.append(new If_Int(GET_SET_CB(int, "prg_", nex)));
  _ifList.append(new If_Int(GET_SET_CB(int, "prg_", multistage)));
  _ifList.append(new If_RealMat(GET_SET_CB(const MATP, "prg_", ys_ref)));
  _ifList.append(new If_RealMat(GET_CB(const MATP, "prg_", M)));
  _ifList.append(new If_RealMat(GET_CB(const MATP, "prg_", P)));
  _ifList.append(new If_RealMat(GET_CB(const MATP, "prg_", COV)));
}

//--------------------------------------------------------------------------
Prg_DynamicEst::~Prg_DynamicEst()
{
  m_free(_dxfdpx0);
  m_free(_dxfdx);
  m_free(_dydpx0);
  m_free(_dxdpx0);
  m_free(_dydx);
  m_free(_COV);
  m_free(_P2);
  m_free(_M2);
  m_free(_ys_ref);
  m_free(_mdl_ys);
  m_free(_mdl_xs);
  m_free(_mdl_us);
  m_free(_mdl_x0s);
  iv_free(_exs);
  v_free(_mdl_p_nominal);
  v_free(_mdl_p_confidence);
  iv_free(_mdl_p_active);
  v_free(_mdl_der_x0_max);
  v_free(_mdl_der_x0_min);
  v_free(_mdl_x0_confidence);
  iv_free(_mdl_x0_active);
  iv_free(_mdl_u_order);
  iv_free(_mdl_y_active);
}

//--------------------------------------------------------------------------
void Prg_DynamicEst::write_active_mx_args(VECP p)
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
void Prg_DynamicEst::setup_model()
{
  if (_mdl_logging >= If_LogInfo)
    If_Log("Info", "Prg_DynamicEst::setup_model");

  // load S-function
  Omu_Model::setup_model(_t0);

  // check for optional S-function methods that are required
  if (_mdl_nx > _mdl_nd)
    assert(ssGetmdlDerivatives(_SS) != NULL);

  // adapt sizes of model vectors
  _mdl_p.alloc(_mdl_np);
  v_zero(_mdl_p.min);

  _mdl_x0.alloc(_mdl_nx);
  _mdl_x.alloc(_mdl_nx);

  iv_resize(_mdl_p_active, _mdl_np);
  v_resize(_mdl_p_confidence, _mdl_np);
  iv_resize(_mdl_x0_active, _mdl_nx);
  v_resize(_mdl_x0_confidence, _mdl_nx);
  v_resize(_mdl_der_x0_min, _mdl_nx);
  v_resize(_mdl_der_x0_max, _mdl_nx);
  iv_resize(_mdl_u_order, _mdl_nu);
  iv_resize(_mdl_y_active, _mdl_ny);
  v_resize(_mdl_p_nominal, _mdl_np);
  iv_zero(_mdl_p_active);
  v_zero(_mdl_p_confidence);
  iv_zero(_mdl_x0_active);
  v_zero(_mdl_x0_confidence);
  v_set(_mdl_der_x0_min, -Inf);
  v_set(_mdl_der_x0_max, Inf);
  iv_set(_mdl_u_order, 1);
  iv_zero(_mdl_y_active);
  if (_mdl_is_fmu) {
    // take over default values from model description
    if(Tcl_VarEval(theInterp, "mdl_p_nominal ${::fmu::", _mdl_name,
                   "::parameterNominalValues}", NULL) != TCL_OK)
      m_error(E_INTERN, "can't obtain nominal parameter values");
  }
  else {
    v_ones(_mdl_p_nominal);
  }
}

//--------------------------------------------------------------------------
void Prg_DynamicEst::setup_stages(IVECP ks, VECP ts)
{
  int kk, i, j;

  // setup S-function
  if (_mdl_needs_setup)
    setup_model();

  // setup optimization problem
  if (_multistage) {
    _K = max(_K, _KK); // assume that one of both has been specified
    _KK = _K;
    stages_alloc(ks, ts, _K, 1);
  }
  else {
    if (_nex > 1) {
      m_error(E_FORMAT, "Prg_DynamicEst::setup_stages: "
	      "prg_nex>1 requires prg_multistage=true");
    }
    _KK = max(_K, _KK); // assume that one of both has been specified
    _K = 1;
    stages_alloc(ks, ts, 1, _KK);
  }

  // initialize _x0 with start values from model and _x0s with _x0
  v_copy(Omu_Model::_mdl_x_start, _mdl_x0);
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
  v_copy(Omu_Model::_mdl_p, _mdl_p);

  // setup _mdl_xs with start values from model
  for (kk = 0; kk < _KK; kk++) {
    for (j = 0; j < _mdl_nx; j++)
      _mdl_xs[kk][j] = _mdl_x0[j];
  }

  // setup _mdl_us with start values from model
  for (kk = 0; kk < _KK; kk++) {
    for (j = 0; j < _mdl_nu; j++)
      _mdl_us[kk][j] = Omu_Model::_mdl_u_start[j];
  }
}

//--------------------------------------------------------------------------
void Prg_DynamicEst::setup(int k,
			   Omu_VariableVec &x, Omu_VariableVec &u,
			   Omu_VariableVec &c)
{
  int i, idx, kk;

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
    // take over possibly modified _mdl_p
    write_mx_args(_mdl_p);

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

    // adapt sizes that depend on active variables
    m_resize(_ys_ref, _KK+1, _ny);

    m_resize(_M2, (_KK+1)*_ny, _np+_nx0);
    m_resize(_P2, _np+_nx0, _np+_nx0);
    m_resize(_COV, _np+_nx0, _np+_nx0);
    m_resize(_dydx, _ny, _nx);
    m_resize(_dxdpx0, _nx, _np+_nx0);
    m_resize(_dydpx0, _ny, _np+_nx0);
    m_resize(_dxfdx, _nx, _nx);
    m_resize(_dxfdpx0, _nx, _np+_nx0);
    m_zero(_M2);
    m_zero(_P2);
    m_zero(_COV);
    _ssr = 0.0;
  }

  if (!_multistage && _nx0 > 0)
    m_error(E_FORMAT, "Prg_DynamicEst::setup: "
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

  if (k == 0) {
    // setup initial states
    for (i = _np; i < _nx; i++) {
      x.initial[i] = _mdl_xs[ks(0)][i-_np] / _mdl_x_nominal[i-_np];
      if (!_mdl_x0_active[i-_np])
	x.min[i] = x.max[i] = x.initial[i];
      else {
	if (_mdl_x0.min[i-_np] > _mdl_x.min[i-_np])
	  x.min[i] = _mdl_x0.min[i-_np] / _mdl_x_nominal[i-_np];
	else
	  x.min[i] = _mdl_x.min[i-_np] / _mdl_x_nominal[i-_np];
	if (_mdl_x0.max[i-_np] < _mdl_x.max[i-_np])
	  x.max[i] = _mdl_x0.max[i-_np] / _mdl_x_nominal[i-_np];
	else
	  x.max[i] = _mdl_x.max[i-_np] / _mdl_x_nominal[i-_np];
      }
    }
  }
  else if (new_experiment) {
    // setup initial states for subsequent experiments
    for (i = 0; i < _mdl_nx; i++) {
      u.initial[i] = _mdl_xs[ks(k)][i] / _mdl_x_nominal[i];
      if (!_mdl_x0_active[i])
	u.min[i] = u.max[i] = u.initial[i];
      else {
	if (_mdl_x0.min[i] > _mdl_x.min[i])
	  u.min[i] = _mdl_x0.min[i] / _mdl_x_nominal[i];
	else
	  u.min[i] = _mdl_x.min[i] / _mdl_x_nominal[i];
	if (_mdl_x0.max[i] < _mdl_x.max[i])
	  u.max[i] = _mdl_x0.max[i] / _mdl_x_nominal[i];
	else
	  u.max[i] = _mdl_x.max[i] / _mdl_x_nominal[i];
      }
    }
  }
  // setup states of subsequent stages
  else {
    for (idx = 0; idx < _mdl_nx; idx++) {
      x.initial[_np+idx] = _mdl_xs[ks(k)][idx] / _mdl_x_nominal[idx];
      x.min[_np+idx] = _mdl_x.min[idx] / _mdl_x_nominal[idx];
      x.max[_np+idx] = _mdl_x.max[idx] / _mdl_x_nominal[idx];
    }
  }
  // constraints for assigning relevant model outputs to opt vars
  // and on time derivatives of initial states
  int nc = _ny;
  if (new_experiment) {
    for (idx = 0; idx < _mdl_nx; idx++) {
      if (_mdl_x0_active[idx]
          && (_mdl_der_x0_min[idx] > -Inf || _mdl_der_x0_max[idx] < Inf))
	++nc;
    }
  }
  c.alloc(nc);
  if (new_experiment) {
    for (i = _ny, idx = 0; idx < _mdl_nx; idx++) {
      if (_mdl_x0_active[idx]
          && (_mdl_der_x0_min[idx] > -Inf || _mdl_der_x0_max[idx] < Inf)) {
        if (_mdl_der_x0_min[idx] > -Inf)
          c.min[i] = _mdl_der_x0_min[idx] / _mdl_x_nominal[idx];
        if (_mdl_der_x0_max[idx] < Inf)
          c.max[i] = _mdl_der_x0_max[idx] / _mdl_x_nominal[idx];
        ++i;
      }
    }
  }
}

//--------------------------------------------------------------------------
void Prg_DynamicEst::setup_struct(int k,
				  const Omu_VariableVec &x,
				  const Omu_VariableVec &u,
				  Omu_DependentVec &xt, Omu_DependentVec &F,
				  Omu_DependentVec &f,
				  Omu_Dependent &f0, Omu_DependentVec &c)
{
  int i, j;
  int ex = _exs[ks(k)]; 	// experiment of current stage
  bool new_experiment = k == 0 || ex != _exs[ks(k)-1];

  // consistic just takes states from optimizer
  // note: possible changes due to discrete events (mdlUpdate) are neglected
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

  if (k < _K) {
    int ex1 = _exs[ks(k+1)]; 	// experiment of subsequent stage

    if (ex == ex1) {
      // explicit ODE for continuous-time equations
      for (i = 0; i < _np + _mdl_nd; i++) {
	for (j = 0; j < _nx; j++)
	  F.Jx[i][j] = 0.0;
      }

      m_zero(F.Jdx);
      for (i = _np + _mdl_nd; i < _nx; i++)
	F.Jdx[i][i] = -1.0;
      F.set_linear(Omu_Dependent::WRT_dx);
    }
    else {
      // no continuous-time equations between two experiments
      m_zero(F.Jx);
      m_zero(F.Jdx);
    }
    m_zero(F.Ju);

    if (_mdl_nd == 0) {
      // f in update depends linearly on x for parameters
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
    }
    v_zero(f0.gxf);
    f0.set_linear(Omu_Dependent::WRT_xf);
    m_zero(c.Jxf);
    c.set_linear(Omu_Dependent::WRT_xf);
  }
}

//--------------------------------------------------------------------------
void Prg_DynamicEst::init_simulation(int k,
				     Omu_VariableVec &x, Omu_VariableVec &u)
{
  if (_mdl_logging >= If_LogInfo)
    If_Log("Info", "Prg_DynamicEst::init_simulation at k = %d", k);

  // initial states
  if (k == 0 || _multistage > 1)
    v_copy(x.initial, x);
}

//--------------------------------------------------------------------------
void Prg_DynamicEst::update(int kk,
			    const Omu_StateVec &x, const Omu_Vec &u,
			    const Omu_StateVec &xf,
			    Omu_DependentVec &f, Omu_Dependent &f0,
			    Omu_DependentVec &c)
{
  if (_mdl_logging >= If_LogInfo)
    If_Log("Info", "Prg_DynamicEst::update at kk = %d", kk);

  int i, j, idx;
  int ex = _exs[kk];
  bool new_experiment = kk == 0 || ex != _exs[kk-1];

  // set simulation time
  ssSetT(_SS, ts(kk));

  // initialize model inputs
  real_T *mdl_u = NULL;
  if (ssGetNumInputPorts(_SS) > 0) {
    if (ssGetInputPortRequiredContiguous(_SS, 0))
      mdl_u = (real_T *)ssGetInputPortRealSignal(_SS, 0);
    else
      mdl_u = (real_T *)*ssGetInputPortRealSignalPtrs(_SS, 0);
  }
  for (idx = 0; idx < _mdl_nu; idx++)
    mdl_u[idx] = _mdl_us[kk][idx];

  // pass optimized parameters to model
  // Note: this is done in any call for numerical approximation of Jacobian
  // Note: initialization is required if model copies mx_args to private memory
  if (_np > 0) {
    write_active_mx_args(x);
    //if (ssGetmdlInitializeConditions(_SS) != NULL)
    //  SMETHOD_CALL(mdlInitializeConditions, _SS);
  }

  // pass current states to model
  real_T *mdl_xd = ssGetDiscStates(_SS);
  real_T *mdl_xc = ssGetContStates(_SS);
  if (kk == 0 || !new_experiment) {
    for (i = 0; i < _mdl_nd; i++)
      mdl_xd[i] = x[_np + i] * _mdl_x_nominal[i];
    for (i = _mdl_nd; i < _mdl_nx; i++)
      mdl_xc[i - _mdl_nd] = x[_np + i] * _mdl_x_nominal[i];
  }
  else {
    for (i = 0; i < _mdl_nd; i++)
      mdl_xd[i] = u[i] * _mdl_x_nominal[i];
    for (i = _mdl_nd; i < _mdl_nx; i++)
      mdl_xc[i - _mdl_nd] = u[i] * _mdl_x_nominal[i];
  }

  // obtain model outputs
  SMETHOD_CALL2(mdlOutputs, _SS, 0);

  // store outputs in constraints
  real_T *mdl_y = NULL;
  if (ssGetNumOutputPorts(_SS) > 0)
    mdl_y = ssGetOutputPortRealSignal(_SS, 0);
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
  if (new_experiment && _mdl_nx > _mdl_nd) {
    real_T *mdl_dx = ssGetdX(_SS);
    SMETHOD_CALL(mdlDerivatives, _SS);
    for (i = _ny, idx = _mdl_nd; idx < _mdl_nx; idx++) {
      if (_mdl_x0_active[idx]
          && (_mdl_der_x0_min[idx] > -Inf || _mdl_der_x0_max[idx] < Inf))
	c[i++] = mdl_dx[idx - _mdl_nd] / _mdl_x_nominal[idx];
    }
  }

  // store values of model states
  for (idx = 0; idx < _mdl_nd; idx++)
    _mdl_xs[kk][idx] = value(mdl_xd[idx]);
  for (idx = _mdl_nd; idx < _mdl_nx; idx++)
    _mdl_xs[kk][idx] = value(mdl_xc[idx - _mdl_nd]);

  // store values of model outputs
  for (idx = 0; idx < _mdl_ny; idx++)
    _mdl_ys[kk][idx] = value(mdl_y[idx]);

  if (new_experiment) {
    // store initial model states
    for (idx = 0; idx < _mdl_nd; idx++)
      _mdl_x0s[ex][idx] = mdl_xd[idx];
    for (idx = _mdl_nd; idx < _mdl_nx; idx++)
      _mdl_x0s[ex][idx] = mdl_xc[idx - _mdl_nd];
  }

  if (kk < _KK) {
    // discrete-time state equations for parameters
    for (i = 0; i < _np; i++)
      f[i] = x[i];
    if (ex == _exs[kk+1]) {
      if (_mdl_nd > 0) {
        // call mdlUpdate to get discrete events processed
        if (ssGetmdlUpdate(_SS) != NULL) {
          setContinuousTask(false);
          setSampleHit(true);
          if (_mdl_is_fmu) {
            // obtain discrete states at end of sample interval
            ssSetT(_SS, ts(kk + 1));
            SMETHOD_CALL2(mdlOutputs, _SS, 0);
          }
          SMETHOD_CALL2(mdlUpdate, _SS, 0);
          setSampleHit(false);
          setContinuousTask(true);
        }
        // read discrete states from model
        for (i = 0; i < _mdl_nd; i++) {
          f[_np + i] = mdl_xd[i] / _mdl_x_nominal[i];
        }
      }
      // junction conditions for continuous-time state equations
      for (i = _np + _mdl_nd; i < _nx; i++)
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

    // call predefined update for numerical differentiation
    Omu_Program::update_grds(kk, x, u, xf, f, f0, c);

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
    // compute Measurement, Precision and Covariance matrix
    //

    if (_mdl_nd > 0)
      // no covariance matrix supported yet
      return;

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
	_M2[kk*_ny + i][j] = _dydpx0[i][j];

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

    // store sum of residuals
    if (kk == 0)
      _ssr = 0.0;
    _ssr += f0;

    // obtain precision matrix, covariance matrix, and confidence intervals
    // Note: skip if nothing is estimated as m_inverse crashes for dim=0
    if (kk == _KK && _M2->n > 0) {
      static double tn[] = {1, 5, 10, 15, 20, 30, 40,
			    50, 60, 80, 100, 200, 500, 100000};
      static double tv[] = {12.706, 2.571, 2.228, 2.131, 2.086, 2.042, 2.021,
			    2.009, 2.000, 1.990,  1.984, 1.972, 1.965, 1.960};
      double n = _ny*(_KK+1) - _np - _nex*_nx0 - 1;
      double t_f2_95;
      if (n < 1)
	t_f2_95 = Inf;
      else if (n >= tn[sizeof(tn)/sizeof(double)-1])
	t_f2_95 = tv[sizeof(tv)/sizeof(double)-1];
      else {
	for (i = 0; i < (int)(sizeof(tn)/sizeof(double)-1); i++) {
	  if (tn[i] <= n && n <= tn[i+1])
	    break;
	}
	double rn = (n - tn[i]) / (tn[i+1] - tn[i]);
	t_f2_95 = tv[i] * (1.0 - rn) + tv[i+1] * rn;
      }
      mtrm_mlt(_M2, _M2, _COV);
      m_catchall(// try
                 m_inverse(_COV, _P2),
                 // catch
                 m_error(E_CONV, 
                         "Prg_DynamicEst::update: singular Covariance matrix")
                 );
      sm_mlt(_ssr/n, _P2, _COV);
      for (i = 0, idx = 0; idx < _mdl_np; idx++) {
	if (_mdl_p_active[idx]) {
	  _mdl_p_confidence[idx]
	    = t_f2_95 * sqrt(_COV[i][i]) * _mdl_p_nominal[idx];
	  i++;
	}
      }
      for (idx = 0; idx < _mdl_nx; idx++) {
	if (_mdl_x0_active[idx]) {
	  _mdl_x0_confidence[idx]
	    = t_f2_95 * sqrt(_COV[i][i]) * _mdl_x_nominal[idx];
	  i++;
	}
      }
    }
  }
}

//--------------------------------------------------------------------------
void Prg_DynamicEst::consistic(int kk, double t,
			       const Omu_StateVec &x, const Omu_Vec &u,
			       Omu_DependentVec &xt)
{
  if (_mdl_logging >= If_LogInfo)
    If_Log("Info", "Prg_DynamicEst::consistic at kk = %d, t = %.3f", kk, t);

  int i;
  int ex = _exs[kk];
  bool new_experiment = kk == 0 || ex != _exs[kk-1];

  // set simulation time
  ssSetT(_SS, t);

  // initialize model inputs
  real_T *mdl_u = NULL;
  if (ssGetNumInputPorts(_SS) > 0) {
    if (ssGetInputPortRequiredContiguous(_SS, 0))
      mdl_u = (real_T *)ssGetInputPortRealSignal(_SS, 0);
    else
      mdl_u = (real_T *)*ssGetInputPortRealSignalPtrs(_SS, 0);
  }
  for (i = 0; i < _mdl_nu; i++)
    mdl_u[i] = _mdl_us[kk][i];

  // pass optimized parameters to model
  // Note: this is done in any call for numerical approximation of Jacobian
  // Note: initialization is required if model copies mx_args to private memory
  if (_np > 0) {
    write_active_mx_args(x);
    //if (ssGetmdlInitializeConditions(_SS) != NULL)
    //  SMETHOD_CALL(mdlInitializeConditions, _SS);
  }
  // initialize model in first stage
  // This way model initial conditions get applied as functions of parameters.
  // Don't initialize if time changed, e.g. for subsequent simulation calls
  if ((_mdl_needs_init || kk == 0 && t == _t0_setup_model || kk != 0 && new_experiment)
      && ssGetmdlInitializeConditions(_SS) != NULL) {
    // initialize model
    if (kk < _KK)
      setSampleTime(ts(kk+1) - ts(kk));
    SMETHOD_CALL(mdlInitializeConditions, _SS);
    _mdl_needs_init = false;
  }

  // disable discrete sample times for continuous evaluation
  setSampleHit(false);
  // mark continuous sample time and process mdlUpdate in case of success
  // or after initialization of FMU
  if (setContinuousTask(true) || _mdl_is_fmu && new_experiment) {
    // pass states from optimizer to model
    real_T *mdl_xd = ssGetDiscStates(_SS);
    real_T *mdl_xc = ssGetContStates(_SS);
    if (kk == 0 || !new_experiment) {
      for (i = 0; i < _mdl_nd; i++)
        mdl_xd[i] = x[_np + i] * _mdl_x_nominal[i];
      for (i = _mdl_nd; i < _mdl_nx; i++)
        mdl_xc[i - _mdl_nd] = x[_np + i] * _mdl_x_nominal[i];
    }
    else {
      for (i = 0; i < _mdl_nd; i++)
        mdl_xd[i] = u[i] * _mdl_x_nominal[i];
      for (i = _mdl_nd; i < _mdl_nx; i++)
        mdl_xc[i - _mdl_nd] = u[i] * _mdl_x_nominal[i];
    }

    // call mdlUpdate to get discrete events processed
    // Note: this is done once at the beginning of a sample interval;
    // no event processing takes place during the integration.
    if (ssGetmdlUpdate(_SS) != NULL) {
      if (_mdl_is_fmu)
        setSampleHit(true);
      // also call mdlOutputs as done by Simulink before each mdlUpdate
      SMETHOD_CALL2(mdlOutputs, _SS, 0); 
      SMETHOD_CALL2(mdlUpdate, _SS, 0);
      if (_mdl_is_fmu)
        setSampleHit(false);
    }

    // take over estimated parameters from optimizer
    for (i = 0; i < _np; i++)
      xt[i] = x[i];

    // read back states from model
    // Note: take estimated initial states from the optimizer
    // to prevent their overriding by the model, e.g. from parameters;
    // they are read once in Omu_Model::setup_model() instead.
    for (i = 0; i < _mdl_nx; i++) {
      if (new_experiment && _mdl_x0_active[i]) {
        if (kk == 0)
          xt[_np + i] = x[_np + i];
        else
          xt[_np + i] = u[i];
      }
      else {
        if (i < _mdl_nd)
          xt[_np + i] = mdl_xd[i] / _mdl_x_nominal[i];
        else
          xt[_np + i] = mdl_xc[i - _mdl_nd] / _mdl_x_nominal[i];
      }
    }
  }
  else {
    // take over estimated parameters from optimizer
    for (i = 0; i < _np; i++)
      xt[i] = x[i];
    // take over regular and estimated states from optimizer
    for (; i < _nx; i++) {
      if (kk == 0 || !new_experiment)
        xt[i] = x[i];
      else
        xt[i] = u[i - _np];
    }
  }
}

//--------------------------------------------------------------------------
void Prg_DynamicEst::continuous(int kk, double t,
				const Omu_StateVec &x, const Omu_Vec &u,
				const Omu_StateVec &dx, Omu_DependentVec &F)
{
  if (_mdl_logging >= If_LogInfo)
    If_Log("Info", "Prg_DynamicEst::continuous at kk = %d, t = %.3f", kk, t);

  int i;

  // set simulation time
  ssSetT(_SS, t);

  // initialize model inputs using linear interpolation over time
  double rt = (t - ts(kk)) / (ts(kk+1) - ts(kk));
  real_T *mdl_u = NULL;
  if (ssGetNumInputPorts(_SS) > 0) {
    if (ssGetInputPortRequiredContiguous(_SS, 0))
      mdl_u = (real_T *)ssGetInputPortRealSignal(_SS, 0);
    else
      mdl_u = (real_T *)*ssGetInputPortRealSignalPtrs(_SS, 0);
  }
  for (i = 0; i < _mdl_nu; i++) {
    if (_mdl_u_order[i] == 0)
      mdl_u[i] = _mdl_us[kk][i];
    else
      mdl_u[i] = _mdl_us[kk][i] * (1 - rt) + _mdl_us[kk+1][i] * rt;
  }

  // pass optimized parameters to model
  // Note: this is done in any call for numerical approximation of Jacobian
  // Note: initialization is required if model copies mx_args to private memory
  if (_np > 0) {
    write_active_mx_args(x);
    //if (ssGetmdlInitializeConditions(_SS) != NULL)
    //  SMETHOD_CALL(mdlInitializeConditions, _SS);
  }

  // pass current states to model
  real_T *mdl_xd = ssGetDiscStates(_SS);
  real_T *mdl_xc = ssGetContStates(_SS);
  for (i = 0; i < _mdl_nd; i++)
    mdl_xd[i] = x[_np + i] * _mdl_x_nominal[i];
  for (i = _mdl_nd; i < _mdl_nx; i++)
    mdl_xc[i - _mdl_nd] = x[_np + i] * _mdl_x_nominal[i];

  // set model outputs before calculating derivatives
  // (note: this is required for sub-blocks providing inputs other blocks)
  SMETHOD_CALL2(mdlOutputs, _SS, 0); 

  // evaluate continuous model equations
  if (_mdl_nx > _mdl_nd)
    SMETHOD_CALL(mdlDerivatives, _SS);

  // get model derivatives and change to residual form
  real_T *mdl_dx = ssGetdX(_SS);
  for (i = _mdl_nd; i < _mdl_nx; i++)
    F[_np + i] = mdl_dx[i - _mdl_nd]/_mdl_x_nominal[i] - dx[_np + i];

  // obtain Jacobians if required
  if (F.is_required_J()) {
    // call predefined continuous for numerical differentiation
    Omu_Program::continuous_grds(kk, t, x, u, dx, F);
  }
}


//==========================================================================

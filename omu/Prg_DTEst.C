/*
 * Prg_DTEst.C -- class definition
 *
 */

/*
    Copyright (C) 1997--2019  Ruediger Franke

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

#include "Prg_DTEst.h"
#include "Hqp_omp.h"

#include <stdlib.h>

#include <If_Real.h>
#include <If_RealVec.h>
#include <If_RealMat.h>
#include <If_Int.h>
#include <If_IntVec.h>
#include <If_Bool.h>
#include <If_Method.h>

#include <meschach/matrix2.h> // for confidence intervals

// redefine assert to throw an error instead of aborting
#undef assert
#define assert(expr) if (!(expr)) m_error(E_INTERN, "assert(" #expr ")");

#define GET_SET_CB(vartype, prefix, name) \
  GET_CB(vartype, prefix, name), \
  IF_SET_CB(vartype, Prg_DTEst, set_##name)

#define GET_CB(vartype, prefix, name) \
  prefix#name, \
  IF_GET_CB(vartype, Prg_DTEst, name)

// Throw E_CONV as errors occuring during solution are generally
// due to bad values and need to be treated (esp. during stepsize check).
#undef SMETHOD_ERROR
#define SMETHOD_ERROR E_CONV

typedef If_Method<Prg_DTEst> If_Cmd;

IF_CLASS_DEFINE("DTEst", Prg_DTEst, Hqp_SqpProgram);

//--------------------------------------------------------------------------
Prg_DTEst::Prg_DTEst()
/* Note: multiple threads not supported my _mx_args and _COV
  : Hqp_Docp(omp_get_max_threads())
  , Omu_Model(ncpu())
*/
{
  _stages_ok = false;
  _within_grds = false;
  _ad = true;
  _fscale = 1.0;

  _mdl_p_active = iv_get(_mdl_np);
  _mdl_p_confidence = v_get(_mdl_np);
  _mdl_x0_active = iv_get(_mdl_nx);
  _mdl_x0_confidence = v_get(_mdl_nx);
  _mdl_u_order = iv_get(_mdl_nu);
  _mdl_y_active = iv_get(_mdl_ny);
  _mdl_p_nominal = v_get(_mdl_np);
  iv_zero(_mdl_p_active);
  v_zero(_mdl_p_confidence);
  iv_zero(_mdl_x0_active);
  v_zero(_mdl_x0_confidence);
  iv_set(_mdl_u_order, 1);
  iv_zero(_mdl_y_active);
  v_ones(_mdl_p_nominal);

  _np = 0;
  _nx0 = 0;
  _ny = 0;
  _nx = _mdl_nx;

  _K = 0;
  _nex = 1;
  _mdl_x0s = m_get(_nex, _mdl_nx);
  _exs = iv_get(_K+1);
  iv_zero(_exs);

  _mdl_us = m_get(_K+1, _mdl_nu);
  _mdl_xs = m_get(_K+1, _mdl_nx);
  _mdl_ys = m_get(_K+1, _mdl_ny);
  _ys_ref = m_get(_K+1, _ny);
  _t0 = 0.0;
  _tf = 1.0;
  _ts = v_get(_K+1);

  _ssr = 0.0;
  _M2 = m_resize(m_get(1, 1), (_K+1)*_ny, _np+_nx0);
  _P2 = m_resize(m_get(1, 1), _np+_nx0, _np+_nx0);
  _COV = m_resize(m_get(1, 1), _np+_nx0, _np+_nx0);
  _dydx = m_resize(m_get(1, 1), _ny, _nx);
  _dxdpx0 = m_resize(m_get(1, 1), _nx, _np+_nx0);
  _dydpx0 = m_resize(m_get(1, 1), _ny, _np+_nx0);
  _dfdx = m_resize(m_get(1, 1), _nx, _nx);
  _dfdpx0 = m_resize(m_get(1, 1), _nx, _np+_nx0);

  _ifList.append(new If_Cmd("prg_setup_stages",
                            &Prg_DTEst::setup_stages, this));
  _ifList.append(new If_Int(GET_SET_CB(int, "prg_", K)));
  _ifList.append(new If_Real(GET_SET_CB(double, "prg_", t0)));
  _ifList.append(new If_Real(GET_SET_CB(double, "prg_", tf)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "prg_", ts)));
  _ifList.append(new If_Bool(GET_SET_CB(bool, "prg_", ad)));
  _ifList.append(new If_Real(GET_SET_CB(double, "prg_", fscale)));

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

  _ifList.append(new If_IntVec(GET_SET_CB(const IVECP, "", mdl_u_order)));

  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_x_min)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_x_max)));

  _ifList.append(new If_IntVec(GET_SET_CB(const IVECP, "", mdl_y_active)));

  _ifList.append(new If_RealMat(GET_SET_CB(const MATP, "", mdl_x0s)));
  _ifList.append(new If_RealMat(GET_SET_CB(const MATP, "", mdl_us)));
  _ifList.append(new If_RealMat(GET_SET_CB(const MATP, "", mdl_xs)));
  _ifList.append(new If_RealMat(GET_CB(const MATP, "", mdl_ys)));

  _ifList.append(new If_Int(GET_SET_CB(int, "prg_", nex)));
  _ifList.append(new If_RealMat(GET_SET_CB(const MATP, "prg_", ys_ref)));
  _ifList.append(new If_RealMat(GET_CB(const MATP, "prg_", M)));
  _ifList.append(new If_RealMat(GET_CB(const MATP, "prg_", P)));
  _ifList.append(new If_RealMat(GET_CB(const MATP, "prg_", COV)));
}

//--------------------------------------------------------------------------
Prg_DTEst::~Prg_DTEst()
{
  m_free(_dfdpx0);
  m_free(_dfdx);
  m_free(_dydpx0);
  m_free(_dxdpx0);
  m_free(_dydx);
  m_free(_COV);
  m_free(_P2);
  m_free(_M2);
  v_free(_ts);
  m_free(_ys_ref);
  m_free(_mdl_ys);
  m_free(_mdl_xs);
  m_free(_mdl_us);
  m_free(_mdl_x0s);
  iv_free(_exs);
  v_free(_mdl_p_nominal);
  v_free(_mdl_p_confidence);
  iv_free(_mdl_p_active);
  v_free(_mdl_x0_confidence);
  iv_free(_mdl_x0_active);
  iv_free(_mdl_u_order);
  iv_free(_mdl_y_active);
}

//--------------------------------------------------------------------------
void Prg_DTEst::write_active_mx_args(VECP p)
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
void Prg_DTEst::setup_model()
{
  if (_mdl_logging >= If_LogInfo)
    If_Log("Info", "Prg_DTEst::setup_model");

  // limit number of CPUs to number of stages
  if (ncpu() > _K + 1)
    set_ncpu(_K + 1);

  // load FMU or S-function and adapt ncpu
  // re-allocate model vectors if new model was loaded
  if (Omu_Model::setup_model(_t0, ncpu())) {

    // no continuous states
    assert(_mdl_nd == _mdl_nx);

    // adapt sizes of model vectors
    _mdl_p.alloc(_mdl_np);
    v_zero(_mdl_p.min);

    _mdl_x0.alloc(_mdl_nx);
    _mdl_x.alloc(_mdl_nx);

    iv_resize(_mdl_p_active, _mdl_np);
    v_resize(_mdl_p_confidence, _mdl_np);
    iv_resize(_mdl_x0_active, _mdl_nx);
    v_resize(_mdl_x0_confidence, _mdl_nx);
    iv_resize(_mdl_u_order, _mdl_nu);
    iv_resize(_mdl_y_active, _mdl_ny);
    v_resize(_mdl_p_nominal, _mdl_np);
    iv_zero(_mdl_p_active);
    v_zero(_mdl_p_confidence);
    iv_zero(_mdl_x0_active);
    v_zero(_mdl_x0_confidence);
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
}

//--------------------------------------------------------------------------
void Prg_DTEst::setup_horizon(int &k0, int &kf)
{
  // call setup_stages if this wasn't already done through the interface
  if (!_stages_ok)
    setup_stages();

  _stages_ok = false;	// for subsequent modifications

  k0 = 0;
  kf = _K;
}

//--------------------------------------------------------------------------
void Prg_DTEst::setup_stages()
{
  int k, i, j;

  if (_mdl_logging >= If_LogInfo)
    If_Log("Info", "Prg_DTEst::setup_stages");

  // setup model
  setup_model();

  // get model structure
  if (_mdl_nx > _mdl_nd)
    m_error(E_FORMAT, "Prg_DTEst::setup_stages: "
            "prg_name DTEst requires model without continuous states");

  // setup optimization program
  v_resize(_ts, _K + 1);
  _ts[0] = _t0; // explicitly assign first value to avoid div. by zero for _K=0
  for (k = 1; k <= _K; k++)
    _ts[k] = _t0 + (double)k * (_tf - _t0) / (double)_K;

  // initialize _x0 with start values from model and _x0s with _x0
  v_copy(Omu_Model::_mdl_x_start, _mdl_x0);
  m_resize(_mdl_x0s, _nex, _mdl_nx);
  for (i = 0; i < _nex; i++) {
    for (j = 0; j < _mdl_nx; j++)
      _mdl_x0s[i][j] = _mdl_x0[j];
  }
  iv_resize(_exs, _K+1);
  iv_zero(_exs);

  m_resize(_mdl_us, _K+1, _mdl_nu);
  m_resize(_mdl_xs, _K+1, _mdl_nx);
  m_resize(_mdl_ys, _K+1, _mdl_ny);

  // store parameters in _mdl_p
  v_copy(Omu_Model::_mdl_p, _mdl_p);

  // setup _mdl_xs with start values from model
  for (k = 0; k < _K; k++) {
    for (j = 0; j < _mdl_nx; j++)
      _mdl_xs[k][j] = _mdl_x0[j];
  }

  // setup _mdl_us with start values from model
  for (k = 0; k < _K; k++) {
    for (j = 0; j < _mdl_nu; j++)
      _mdl_us[k][j] = Omu_Model::_mdl_u_start[j];
  }

  _stages_ok = true;
}

//--------------------------------------------------------------------------
void Prg_DTEst::setup_vars(int k,
			   VECP x, VECP x_min, VECP x_max, IVECP x_int,
			   VECP u, VECP u_min, VECP u_max, IVECP u_int,
			   VECP c, VECP c_min, VECP c_max)
{
  int i, idx;

  // initialize number of experiment
  int ex = 0;
  bool new_experiment = k == 0 || _ts[k] < _ts[k-1];
  if (k > 0) {
    if (new_experiment)
      // new experiment as time steps back
      ex = _exs[k-1] + 1;
    else
      ex = _exs[k-1];
  }
  assert(ex < _nex);  // _nex and experiment data must fit
  // store number of experiment for stage
  _exs[k] = ex;

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
    m_resize(_ys_ref, _K+1, _ny);

    m_resize(_M2, (_K+1)*_ny, _np+_nx0);
    m_resize(_P2, _np+_nx0, _np+_nx0);
    m_resize(_COV, _np+_nx0, _np+_nx0);
    m_resize(_dydx, _ny, _nx);
    m_resize(_dxdpx0, _nx, _np+_nx0);
    m_resize(_dydpx0, _ny, _np+_nx0);
    m_resize(_dfdx, _nx, _nx);
    m_resize(_dfdpx0, _nx, _np+_nx0);
    m_zero(_M2);
    m_zero(_P2);
    m_zero(_COV);
    _ssr = 0.0;
  }

  // allocate optimization variables
  if (k > 0 && new_experiment) {
    alloc_vars(x, x_min, x_max, x_int, _np);
    alloc_vars(u, u_min, u_max, u_int, _mdl_nx);
  }
  else {
    alloc_vars(x, x_min, x_max, x_int, _nx);
  }

  // setup parameters
  int ip = 0;
  for (i = 0; i < _mdl_np; i++) {
    if (_mdl_p_active[i]) {
      x[ip] = _mdl_p[i] / _mdl_p_nominal[i];
      if (_mdl_p.min[i] > -Inf)
        x_min[ip] = _mdl_p.min[i] / _mdl_p_nominal[i];
      if (_mdl_p.max[i] < Inf)
        x_max[ip] = _mdl_p.max[i] / _mdl_p_nominal[i];
      ++ip;
    }
  }
  assert(ip == _np);	// _np must not have changed since setup_stages

  if (k == 0) {
    // setup initial states
    for (i = _np; i < _nx; i++) {
      x[i] = _mdl_xs[0][i-_np] / _mdl_x_nominal[i-_np];
      if (!_mdl_x0_active[i-_np])
	x_min[i] = x_max[i] = x[i];
      else {
	if (_mdl_x0.min[i-_np] > _mdl_x.min[i-_np])
	  x_min[i] = _mdl_x0.min[i-_np] / _mdl_x_nominal[i-_np];
	else
	  x_min[i] = _mdl_x.min[i-_np] / _mdl_x_nominal[i-_np];
	if (_mdl_x0.max[i-_np] < _mdl_x.max[i-_np])
	  x_max[i] = _mdl_x0.max[i-_np] / _mdl_x_nominal[i-_np];
	else
	  x_max[i] = _mdl_x.max[i-_np] / _mdl_x_nominal[i-_np];
      }
    }
  }
  else if (new_experiment) {
    // setup initial states for subsequent experiments
    for (i = 0; i < _mdl_nx; i++) {
      u[i] = _mdl_xs[k][i] / _mdl_x_nominal[i];
      if (!_mdl_x0_active[i])
	u_min[i] = u_max[i] = u[i];
      else {
	if (_mdl_x0.min[i] > _mdl_x.min[i])
	  u_min[i] = _mdl_x0.min[i] / _mdl_x_nominal[i];
	else
	  u_min[i] = _mdl_x.min[i] / _mdl_x_nominal[i];
	if (_mdl_x0.max[i] < _mdl_x.max[i])
	  u_max[i] = _mdl_x0.max[i] / _mdl_x_nominal[i];
	else
	  u_max[i] = _mdl_x.max[i] / _mdl_x_nominal[i];
      }
    }
  }
  // setup states of subsequent stages
  else {
    for (idx = 0; idx < _mdl_nx; idx++) {
      x[_np+idx] = _mdl_xs[k][idx] / _mdl_x_nominal[idx];
      x_min[_np+idx] = _mdl_x.min[idx] / _mdl_x_nominal[idx];
      x_max[_np+idx] = _mdl_x.max[idx] / _mdl_x_nominal[idx];
    }
  }
  // constraints for assigning relevant model outputs to opt vars
  alloc_vars(c, c_min, c_max, IVNULL, _ny);
}

//--------------------------------------------------------------------------
void Prg_DTEst::setup_struct(int k, const VECP x, const VECP u,
			     MATP fx, MATP fu, IVECP f_lin,
			     VECP f0x, VECP f0u, int &f0_lin,
			     MATP cx, MATP cu, IVECP c_lin,
			     MATP Lxx, MATP Luu, MATP Lxu)
{
  int i, j, idx, jdx, offs;
  int ex = _exs[k];
  SimStruct *S = _SS[0];

  if (_mdl_is_fmu && _mdl_jac && ssGetmdlJacobian(S) != NULL) {
    if (k == 0) {
      for (idx = 0; idx < _mdl_ny; idx++) {
        _mdl_jac_y_active[idx] = _mdl_y_active[idx];
      }
      setup_jac();
    }

    // obtain sparse structure of Jacobian
    m_zero(fx);
    m_zero(fu);
    m_zero(cx);
    m_zero(cu);
    // assume full Jacobian wrt parameters as not available from model
    if (k < _K) {
      for (j = 0; j < _np; j++) {
        fx[j][j] = 1.0;
        if (ex == _exs[k+1]) {
          for (i = _np; i < _nx; i++)
            fx[i][j] = 1.0;
        }
      }
    }
    // obtain structure wrt model states and inputs
    real_T *pr = ssGetJacobianPr(S);
    int i_end = ssGetJacobianJc(S)[_mdl_nx + _mdl_nu];
    for (i = 0; i < i_end; i++)
      *pr++ = 1.0;
    fetch_jacxu(S, k, x, u, fx, fu, cx, cu);
  }
}

//--------------------------------------------------------------------------
void Prg_DTEst::init_simulation(int k, VECP x, VECP u)
{
  if (_mdl_logging >= If_LogInfo)
    If_Log("Info", "Prg_DTEst::init_simulation at k = %d", k);

  // initial states
  if (k == 0) {
    int i, ip = 0;
    for (i = 0; i < _mdl_np; i++) {
      if (_mdl_p_active[i]) {
        x[ip] = _mdl_p[i] / _mdl_p_nominal[i];
        ++ip;
      }
    }
    for (i = _np; i < _nx; i++) {
      x[i] = _mdl_xs[0][i-_np] / _mdl_x_nominal[i-_np];
    }
  }
}

//--------------------------------------------------------------------------
void Prg_DTEst::update_vals(int k, const VECP x, const VECP u,
			    VECP f, Real &f0, VECP c)
{
  int i, j, idx;
  int ex = _exs[k];
  bool new_experiment = k == 0 || ex != _exs[k-1];
  int tn = omp_get_thread_num();
  SimStruct *S = _SS[tn];

  if (_mdl_logging >= If_LogInfo)
    If_Log("Info", "Prg_DTEst::update_vals at k = %d, tn = %d", k, tn);

  // initialize model inputs
  real_T *mdl_u = NULL;
  if (ssGetNumInputPorts(S) > 0) {
    if (ssGetInputPortRequiredContiguous(S, 0))
      mdl_u = (real_T *)ssGetInputPortRealSignal(S, 0);
    else
      mdl_u = (real_T *)*ssGetInputPortRealSignalPtrs(S, 0);
  }
  for (idx = 0; idx < _mdl_nu; idx++)
    mdl_u[idx] = _mdl_us[k][idx];

  // pass optimized parameters to model
  // Note: this is done in any call for numerical approximation of Jacobian
  // Note: initialization is required if model copies mx_args to private memory
  if (_np > 0) {
    write_active_mx_args(x);
    //_mdl_needs_init[tn] = 1;
  }

  // set simulation time
  ssSetT(S, _ts[k]);
  if (k < _K)
    setSampleTime(S, ts(k+1) - ts(k));

  // pass current states to model
  real_T *mdl_xd = ssGetDiscStates(S);
  real_T *mdl_xc = ssGetContStates(S);
  bool needs_update = true;
  bool after_init = false;

  if (k == 0 || !new_experiment) {
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

  // initialize model in first stage
  // Don't initialize if time changed, e.g. for subsequent simulation calls
  if ((_mdl_needs_init[tn] || k == 0 && ts(k) == _t0_setup_model)
      && ssGetmdlInitializeConditions(S) != NULL) {
    // initialize model
    SMETHOD_CALL(mdlInitializeConditions, S);
    _mdl_needs_init[tn] = 0;
    // call mdlOutputs because implementation may be delayed
    after_init = true;
  }

  if (after_init || _mdl_is_fmu && k == 0) {
    // call mdlUpdate and disable continuous task to trigger initial clock
    if (_mdl_is_fmu) {
      setSampleHit(S, true);
      setContinuousTask(S, k != 0);
    }
    SMETHOD_CALL2(mdlOutputs, S, 0);
    if (ssGetmdlUpdate(S) != NULL) {
      SMETHOD_CALL2(mdlUpdate, S, 0);
    }
    if (_mdl_is_fmu) {
      setSampleHit(S, false);
      setContinuousTask(S, true);
    }
    needs_update = false;
    // take over regular and estimated states from optimizer
    for (i = 0; i < _mdl_nd; i++) {
      if (k == 0 || !new_experiment) {
        if (_mdl_x0_active[i] || !new_experiment) {
          mdl_xd[i] = x[_np + i] * _mdl_x_nominal[i];
          needs_update = true;
        }
      }
      else {
        if (_mdl_x0_active[i]) {
          mdl_xd[i] = u[i] * _mdl_x_nominal[i];
          needs_update = true;
        }
      }
    }
  }

  // call mdlOutputs/mdlUpdate to get outputs for current states
  // enable continuous task to not update states here
  if (needs_update) {
    if (_mdl_is_fmu)
      setSampleHit(S, true);
    // also call mdlOutputs as done by Simulink before each mdlUpdate
    SMETHOD_CALL2(mdlOutputs, S, 0);
    if (ssGetmdlUpdate(S) != NULL) {
      SMETHOD_CALL2(mdlUpdate, S, 0);
    }
    if (_mdl_is_fmu)
      setSampleHit(S, false);
  }

  // obtain model outputs
  real_T *mdl_y = NULL;
  if (ssGetNumOutputPorts(S) > 0)
    mdl_y = ssGetOutputPortRealSignal(S, 0);

  // store outputs in constraints
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
      help = c[i] - _ys_ref[k][i]/_mdl_y_nominal[idx];
      f0 += help*help;
      ++i;
    }
  }

  if (!_within_grds) {
    // store values of model states
    for (idx = 0; idx < _mdl_nd; idx++)
      _mdl_xs[k][idx] = value(mdl_xd[idx]);
    for (idx = _mdl_nd; idx < _mdl_nx; idx++)
      _mdl_xs[k][idx] = value(mdl_xc[idx - _mdl_nd]);

    // store values of model outputs
    for (idx = 0; idx < _mdl_ny; idx++)
      _mdl_ys[k][idx] = value(mdl_y[idx]);

    if (new_experiment) {
      // store initial model states
      for (idx = 0; idx < _mdl_nd; idx++)
        _mdl_x0s[ex][idx] = mdl_xd[idx];
      for (idx = _mdl_nd; idx < _mdl_nx; idx++)
        _mdl_x0s[ex][idx] = mdl_xc[idx - _mdl_nd];
    }
  }

  if (k < _K) {
    // discrete-time state equations for parameters
    for (i = 0; i < _np; i++)
      f[i] = x[i];
    if (ex == _exs[k+1]) {
      if (_mdl_nd > 0) {
        // call mdlUpdate to get discrete events processed
        if (ssGetmdlUpdate(S) != NULL) {
          setContinuousTask(S, false);
          setSampleHit(S, true);
          if (_mdl_is_fmu) {
            // obtain discrete states at end of sample interval
            ssSetT(S, _ts[k + 1]);
            SMETHOD_CALL2(mdlOutputs, S, 0);
          }
          SMETHOD_CALL2(mdlUpdate, S, 0);
          setSampleHit(S, false);
          setContinuousTask(S, true);
        }
        // read discrete states from model
        for (i = 0; i < _mdl_nd; i++) {
          f[_np + i] = mdl_xd[i] / _mdl_x_nominal[i];
        }
      }
    }
  }
  f0 *= _fscale;
}

//--------------------------------------------------------------------------
void Prg_DTEst::update_stage(int k, const VECP x, const VECP u,
			     VECP f, Real &f0, VECP c,
			     MATP fx, MATP fu, VECP f0x, VECP f0u,
			     MATP cx, MATP cu,
			     const VECP rf, const VECP rc,
			     MATP Lxx, MATP Luu, MATP Lxu)
{
  int i, j, idx;
  int ex = _exs[k];
  bool new_experiment = k == 0 || ex != _exs[k-1];
  int tn = omp_get_thread_num();
  SimStruct *S = _SS[tn];

  if (_mdl_logging >= If_LogInfo)
    If_Log("Info", "Prg_DTEst::update_stage at k = %d, tn = %d", k, tn);

  // store current model parameters and states
  if (k == 0) {
    // store model parameters
    read_mx_args(_mdl_p);
  }

  if (!_ad || !_mdl_jac || ssGetmdlJacobian(S) == NULL) {
    // call predefined update for numerical differentiation
    update_vals(k, x, u, f, f0, c);
    _within_grds = true;
    update_grds(k, x, u, fx, fu, f0x, f0u, cx, cu);
    update_hela(k, x, u, rf, rc, Lxx, Luu, Lxu);
    _within_grds = false;
    if (!_ad)
      return;
  }
  else {
    // update values
    update_vals(k, x, u, f, f0, c);
    // obtain Jacobians
    m_zero(fx);
    m_zero(fu);
    m_zero(cx);
    m_zero(cu);
    // obtain Jacobian wrt parameters numerically
    Real vj_bak, dvj, df0;
    int nf = f->dim;
    int nc = c->dim;
    VECP df = v_get(nf);
    VECP dc = v_get(nc);
    for (j = 0; j < _np; j++) {
      vj_bak = x[j];
      dvj = 1e-4 * fabs(vj_bak) + 1e-6;
      x->ve[j] += dvj;

      v_zero(df);
      df0 = 0.0;
      v_zero(dc);

      update_vals(k, x, u, df, df0, dc);

      v_sub(df, f, df);
      for (i = 0; i < nf; i++)
        fx[i][j] = df[i] / dvj;

      f0x[j] = (df0 - f0) / dvj;

      v_sub(dc, c, dc);
      for (i = 0; i < nc; i++)
        cx[i][j] = dc[i] / dvj;

      x->ve[j] = vj_bak;
    }
    v_free(dc);
    v_free(df);

    // obtain Jacobian wrt model states and inputs from model
    // first restore states after mdlUpdate has been processed
    // and call mdlOutputs
    real_T *mdl_xd = ssGetDiscStates(S);
    real_T *mdl_xc = ssGetContStates(S);
    if (k == 0 || !new_experiment) {
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
    SMETHOD_CALL2(mdlOutputs, S, 0); 

    SMETHOD_CALL(mdlJacobian, S);

    fetch_jacxu(S, k, x, u, fx, fu, cx, cu);
  }

  // (re-)obtain gradient of f0
  v_zero(f0x);
  v_zero(f0u);
  for (i = 0, idx = 0; idx < _mdl_ny; idx++) {
    if (_mdl_y_active[idx]) {
      for (j = 0; j < _np; j++)
        f0x[j] += 
          2.0 * (c[i] - _ys_ref[k][i]/_mdl_y_nominal[idx]) * cx[i][j];
      if (k == 0 || !new_experiment) {
        for (; j < _nx; j++)
          f0x[j] += 
            2.0 * (c[i] - _ys_ref[k][i]/_mdl_y_nominal[idx]) * cx[i][j];
      }
      else {
        for (; j < _nx; j++)
          f0u[j-_np] += 
            2.0 * (c[i] - _ys_ref[k][i]/_mdl_y_nominal[idx]) * cu[i][j-_np];
      }
      ++i;
    }
  }
  sv_mlt(_fscale, f0x, f0x);
  sv_mlt(_fscale, f0u, f0u);

  //
  // compute Measurement, Precision and Covariance matrix
  //

  // build dx/d(p,x0)
  if (k == 0) {
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
      m_copy(_dfdpx0, _dxdpx0);
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
  if (k == 0 || !new_experiment)
    m_copy(cx, _dydx);
  else {
    // collect dydx from c.Jx and c.Ju
    for (i = 0; i < _ny; i++) {
      for (j = 0; j < _np; j++)
        _dydx[i][j] = cx[i][j];
      for (j = _np; j < _nx; j++)
        _dydx[i][j] = cu[i][j-_np];
    }
  }
  m_mlt(_dydx, _dxdpx0, _dydpx0);

  // store results in M
  for (i = 0; i < _ny; i++)
    for (j = 0; j < _np + _nx0; j++)
      _M2[k*_ny + i][j] = _dydpx0[i][j];

  // consider contribution of continuous-time equations for next k
  if (k < _K && ex == _exs[k+1]) {
    if (k == 0 || !new_experiment)
      m_copy(fx, _dfdx);
    else {
      // collect dfdx from fx and fu
      for (i = 0; i < _nx; i++) {
        for (j = 0; j < _np; j++)
          _dfdx[i][j] = fx[i][j];
        for (j = _np; j < _nx; j++)
          _dfdx[i][j] = fu[i][j-_np];
      }
    }
    m_mlt(_dfdx, _dxdpx0, _dfdpx0);
  }

  // store sum of residuals
  if (k == 0)
    _ssr = 0.0;
  _ssr += f0;

  // obtain precision matrix, covariance matrix, and confidence intervals
  // Note: skip if nothing is estimated as m_inverse crashes for dim=0
  if (k == _K && _M2->n > 0) {
    static double tn[] = {1, 5, 10, 15, 20, 30, 40,
                          50, 60, 80, 100, 200, 500, 100000};
    static double tv[] = {12.706, 2.571, 2.228, 2.131, 2.086, 2.042, 2.021,
                          2.009, 2.000, 1.990,  1.984, 1.972, 1.965, 1.960};
    double n = _ny*(_K+1) - _np - _nex*_nx0 - 1;
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
                       "Prg_DTEst::update_stage: singular Covariance matrix")
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

//--------------------------------------------------------------------------
void Prg_DTEst::fetch_jacxu(SimStruct *S,
                            int k, const VECP x, const VECP u,
                            MATP fx, MATP fu, MATP cx, MATP cu)
{
  int ex = _exs[k];
  bool new_experiment = k == 0 || ex != _exs[k-1];
  int i, idx, j, jdx, rdx;
  int ixf;
  int spsk = 1;
  int upsk = 1;
  int mdl_x_idx, mdl_y_idx;
  int mdl_nc = _mdl_nx - _mdl_nd; // number of continuous states
  real_T *pr = ssGetJacobianPr(S);
  int_T *ir = ssGetJacobianIr(S);
  int_T *jc = ssGetJacobianJc(S);

  // Jacobian wrt S-function states (ddxdy/dx)
  for (jdx = 0; jdx < _mdl_nx; jdx++) {
    j = jdx < mdl_nc? _mdl_nd + jdx: jdx - mdl_nc;
    mdl_x_idx = jdx < mdl_nc? _mdl_nd + jdx: jdx - mdl_nc;
    for (i = 0, idx = _mdl_nx, rdx = jc[jdx]; rdx < jc[jdx+1]; rdx++) {
      if (ir[rdx] < _mdl_nx) {
        if (k < _K && ex == _exs[k+1]) {
          // junction conditions for discrete states
          if (k == 0 || !new_experiment) {
            ixf = ir[rdx] - mdl_nc;
            fx[_np + ixf][_np + j] = pr[rdx] /
              _mdl_x_nominal[ixf] * _mdl_x_nominal[mdl_x_idx];
          }
          else {
            ixf = ir[rdx];
            fu[_np + ixf][j] = pr[rdx] /
              _mdl_x_nominal[ixf] * _mdl_x_nominal[mdl_x_idx];
          }
        }
      }
      else {
        mdl_y_idx = ir[rdx] - _mdl_nx;
        if (_mdl_y_active[mdl_y_idx]) {
          // need to loop through idx to obtain i considering active outputs
          for (; idx < ir[rdx]; idx++) {
            if (_mdl_y_active[idx - _mdl_nx])
              i += spsk;
          }
          if (k == 0 || !new_experiment) {
            cx[i + k%spsk][_np + j] = pr[rdx] /
              _mdl_y_nominal[mdl_y_idx] * _mdl_x_nominal[mdl_x_idx];
          }
          else {
            cu[i + k%spsk][j] = pr[rdx] /
              _mdl_y_nominal[mdl_y_idx] * _mdl_x_nominal[mdl_x_idx];
          }
        }
      }
    }
  }
}

//==========================================================================

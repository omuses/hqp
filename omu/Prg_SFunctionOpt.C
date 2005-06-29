/*
 * Prg_SFunctionOpt.C -- class definition
 *
 */

/*
    Copyright (C) 1997--2005  Ruediger Franke

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

#include <If_Int.h>
#include <If_Bool.h>
#include <If_Real.h>
#include <If_RealVec.h>
#include <If_RealMat.h>
#include <If_IntVec.h>
#include <If_Int.h>
#include <If_IntVec.h>

// redefine assert to throw an error instead of aborting
#undef assert
#define assert(expr) if (!(expr)) m_error(E_INTERN, "assert(" #expr ")");

#define GET_SET_CB(vartype, prefix, name) \
  GET_CB(vartype, prefix, name), \
  IF_SET_CB(vartype, Prg_SFunctionOpt, set_##name)

#define GET_CB(vartype, prefix, name) \
  prefix#name, \
  IF_GET_CB(vartype, Prg_SFunctionOpt, name)

// Call an S-function method and check for errors.
// Throw E_CONV as errors occuring during solution are generally
// due to bad values and need to be treated (esp. during stepsize check).
#define SMETHOD_CALL(method, S) \
  ssSetErrorStatus(S, NULL); \
  method(S); \
  if (ssGetErrorStatus(S)) { \
    fprintf(stderr, "Error from " #method ": %s\n", \
	    ssGetErrorStatus(S)); \
    m_error(E_CONV, ssGetErrorStatus(S)); \
  }

#define SMETHOD_CALL2(method, S, tid) \
  ssSetErrorStatus(S, NULL); \
  method(S, tid); \
  if (ssGetErrorStatus(S)) { \
    fprintf(stderr, "Error from " #method ": %s\n", \
	    ssGetErrorStatus(S)); \
    m_error(E_CONV, ssGetErrorStatus(S)); \
  }

IF_CLASS_DEFINE("SFunctionOpt", Prg_SFunctionOpt, Omu_Program);

//--------------------------------------------------------------------------
Omu_OptVarVec::Omu_OptVarVec()
{
  weight1 = v_resize(v_get(1), 0);
  weight2 = v_resize(v_get(1), 0);
  ref = v_resize(v_get(1), 0);
  active = iv_resize(iv_get(1), 0);
}

//--------------------------------------------------------------------------
Omu_OptVarVec::~Omu_OptVarVec()
{
  v_free(ref);
  v_free(weight2);
  v_free(weight1);
  iv_free(active);
}

//--------------------------------------------------------------------------
void Omu_OptVarVec::resize(int n)
{
  Omu_VariableVec::alloc(n);
  v_resize(weight1, n);
  v_resize(weight2, n);
  v_resize(ref, n);
  iv_resize(active, n);

  v_set(weight1, 0.0);
  v_set(weight2, 0.0);
  v_set(ref, 0.0);
  iv_zero(active);
}

//==========================================================================

//--------------------------------------------------------------------------
Prg_SFunctionOpt::Prg_SFunctionOpt()
{
  _sps = 1;
  _multistage = true;

  _mdl_x0_active = iv_get(_mdl_nx);
  _mdl_der_x0_min = v_get(_mdl_nx);
  _mdl_der_x0_max = v_get(_mdl_nx);
  _mdl_u_order = iv_get(_mdl_nu);
  _mdl_u0_nfixed = iv_get(_mdl_nu);
  _mdl_u_decimation = iv_get(_mdl_nu);
  _mdl_u_nominal = v_get(_mdl_nu);
  _mdl_x_nominal = v_get(_mdl_nx);
  _mdl_y_nominal = v_get(_mdl_ny);
  _t_nominal = 1.0;
  _mdl_y_bias = v_get(_mdl_ny);
  v_set(_mdl_x0, 0.0);
  iv_set(_mdl_x0_active, 0);
  v_set(_mdl_der_x0_min, -Inf);
  v_set(_mdl_der_x0_max, Inf);
  iv_set(_mdl_u_order, 1);
  iv_set(_mdl_u0_nfixed, 0);
  iv_set(_mdl_u_decimation, 1);
  v_set(_mdl_u_nominal, 1.0);
  v_set(_mdl_x_nominal, 1.0);
  v_set(_mdl_y_nominal, 1.0);
  v_set(_mdl_y_bias, 0.0);

  // numbers of optimization variables
  _nu = 0;
  _nx = _nu + _mdl_nx;
  _nc = 0;
  _ns = 0;
  _nsc = 0;
  _ncf = 0;

  _mdl_us = m_get(_KK+1, _mdl_nu);
  _mdl_xs = m_get(_KK+1, _mdl_nx);
  _mdl_ys = m_get(_KK+1, _mdl_ny);

  _ifList.append(new If_Int("prg_sps", &_sps));
  _ifList.append(new If_Bool(GET_SET_CB(bool, "prg_", multistage)));

  // redefine mdl_x0 in order to also consider mdl_xs
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_x0)));
  _ifList.append(new If_IntVec(GET_SET_CB(const IVECP, "", mdl_x0_active)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_x0_min)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_x0_max)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_der_x0_min)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_der_x0_max)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_u0_min)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_u0_max)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_y0_min)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_y0_max)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_y0_weight1)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_y0_weight2)));

  _ifList.append(new If_IntVec(GET_SET_CB(const IVECP, "", mdl_u_order)));
  _ifList.append(new If_IntVec(GET_SET_CB(const IVECP, "", mdl_u_active)));
  _ifList.append(new If_IntVec(GET_SET_CB(const IVECP, "", mdl_u0_nfixed)));
  _ifList.append(new If_IntVec(GET_SET_CB(const IVECP, "", mdl_u_decimation)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_u_nominal)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_u_min)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_u_max)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_u_ref)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_u_weight1)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_u_weight2)));

  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_der_u_min)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_der_u_max)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_der_u_ref)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "",
					   mdl_der_u_weight1)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "",
					   mdl_der_u_weight2)));

  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_x_nominal)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_x_min)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_x_max)));

  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_y_bias)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_y_nominal)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_y_min)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_y_max)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_y_ref)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_y_weight1)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_y_weight2)));

  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_y_soft_min)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_y_soft_max)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "",
					   mdl_y_soft_weight1)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "",
					   mdl_y_soft_weight2)));

  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_yf_min)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_yf_max)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_yf_weight1)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_yf_weight2)));

  _ifList.append(new If_RealMat(GET_SET_CB(const MATP, "", mdl_us)));
  _ifList.append(new If_RealMat(GET_SET_CB(const MATP, "", mdl_xs)));
  _ifList.append(new If_RealMat(GET_CB(const MATP, "", mdl_ys)));
}

//--------------------------------------------------------------------------
Prg_SFunctionOpt::~Prg_SFunctionOpt()
{
  m_free(_mdl_ys);
  m_free(_mdl_xs);
  m_free(_mdl_us);
  v_free(_mdl_y_bias);
  v_free(_mdl_y_nominal);
  v_free(_mdl_x_nominal);
  v_free(_mdl_u_nominal);
  iv_free(_mdl_u_decimation);
  iv_free(_mdl_u0_nfixed);
  iv_free(_mdl_u_order);
  iv_free(_mdl_x0_active);
  v_free(_mdl_der_x0_max);
  v_free(_mdl_der_x0_min);
}

//--------------------------------------------------------------------------
void Prg_SFunctionOpt::setup_model()
{
  // load S-function
  Prg_SFunction::setup_model();

  // check for optional S-function methods that are required
  assert(ssGetmdlDerivatives(_S) != NULL);

  // adapt sizes of model vectors
  _mdl_x0.alloc(_mdl_nx);
  _mdl_u0.alloc(_mdl_nu);
  _mdl_y0.resize(_mdl_ny);
  _mdl_u.resize(_mdl_nu);
  _mdl_der_u.resize(_mdl_nu);
  _mdl_x.alloc(_mdl_nx);
  _mdl_y.resize(_mdl_ny);
  _mdl_y_soft.resize(_mdl_ny);
  _mdl_yf.resize(_mdl_ny);

  iv_resize(_mdl_x0_active, _mdl_nx);
  v_resize(_mdl_der_x0_min, _mdl_nx);
  v_resize(_mdl_der_x0_max, _mdl_nx);
  iv_resize(_mdl_u_order, _mdl_nu);
  iv_resize(_mdl_u0_nfixed, _mdl_nu);
  iv_resize(_mdl_u_decimation, _mdl_nu);
  v_resize(_mdl_u_nominal, _mdl_nu);
  v_resize(_mdl_x_nominal, _mdl_nx);
  v_resize(_mdl_y_nominal, _mdl_ny);
  v_resize(_mdl_y_bias, _mdl_ny);
  iv_set(_mdl_x0_active, 0);
  v_set(_mdl_der_x0_min, -Inf);
  v_set(_mdl_der_x0_max, Inf);
  iv_set(_mdl_u_order, 1);
  iv_set(_mdl_u0_nfixed, 0);
  iv_set(_mdl_u_decimation, 1);
  v_set(_mdl_u_nominal, 1.0);
  v_set(_mdl_x_nominal, 1.0);
  v_set(_mdl_y_nominal, 1.0);
  v_set(_mdl_y_bias, 0.0);
}

//--------------------------------------------------------------------------
void Prg_SFunctionOpt::setup_stages(IVECP ks, VECP ts)
{
  int kk, j;

  // setup S-function
  if (_mdl_needs_setup)
    setup_model();

  // setup optimization problem
  if (_sps < 1) {
    m_error(E_FORMAT, "Prg_SFunctionOpt::setup_stages: "
	    "prg_sps must at least be one");
  }
  if (_multistage) {
    _K = _KK/_sps;
    stages_alloc(ks, ts, _K, _sps);
  }
  else {
    if (_sps > 1) {
      m_error(E_FORMAT, "Prg_SFunctionOpt::setup_stages: "
	      "prg_sps>1 requires prg_multistage=true");
    }
    _K = 1;
    stages_alloc(ks, ts, 1, _KK);
  }

  m_resize(_mdl_us, _KK+1, _mdl_nu);
  m_resize(_mdl_xs, _KK+1, _mdl_nx);
  m_resize(_mdl_ys, _KK+1, _mdl_ny);

  // setup _mdl_xs with initial states from model
  v_copy(Prg_SFunction::_mdl_x0, _mdl_x0);
  for (kk = 0; kk < _KK; kk++) {
    for (j = 0; j < _mdl_nx; j++)
      _mdl_xs[kk][j] = _mdl_x0[j];
  }
}

//--------------------------------------------------------------------------
void Prg_SFunctionOpt::setup(int k,
			     Omu_VariableVec &x, Omu_VariableVec &u,
			     Omu_VariableVec &c)
{
  int i, idx, isc, j;

  // complete general setup
  if (k == 0) {
    // take over possibly modified _mdl_p
    write_mx_args(_mdl_p);

    // obtain numbers of optimization variables (active model variables)
    _nu = 0;
    for (idx = 0; idx < _mdl_nu; idx++) {
      if (_mdl_u.active[idx])
	_nu++;
    }
    _nx = _nu + _mdl_nx;
    _nc = 0;
    _ns = 0;
    _nsc = 0;
    _nc0 = 0;
    _ncf = 0;
    for (idx = 0; idx < _mdl_nx; idx++) {
      if (_mdl_x0_active[idx]
          && (_mdl_der_x0_min[idx] > -Inf || _mdl_der_x0_max[idx] < Inf)) {
	_nc0++;
      }
    }
    for (idx = 0; idx < _mdl_ny; idx++) {
      if (_mdl_y.min[idx] > -Inf || _mdl_y.max[idx] < Inf
	  || _mdl_y.weight1[idx] != 0.0 || _mdl_y.weight2[idx] != 0.0) {
	_mdl_y.active[idx] = 1;
	_nc++;
      }
      if (_mdl_y_soft.min[idx] > -Inf) {
	_nsc++;
      }
      if (_mdl_y_soft.max[idx] < Inf) {
	_nsc++;
      }
      if (_mdl_y_soft.min[idx] > -Inf || _mdl_y_soft.max[idx] < Inf) {
	_mdl_y_soft.active[idx] = 1;
	_ns++;
      }
      if (_mdl_y0.min[idx] > -Inf || _mdl_y0.max[idx] < Inf
	  || _mdl_y0.weight1[idx] != 0.0 || _mdl_y0.weight2[idx] != 0.0) {
	if (!_multistage)
	  m_error(E_FORMAT, "Prg_SFunctionOpt::setup: "
		  "use of mdl_y0 requires prg_multistage=true");
	_mdl_y0.active[idx] = 1;
	_nc0++;
      }
      if (_mdl_yf.min[idx] > -Inf || _mdl_yf.max[idx] < Inf
	  || _mdl_yf.weight1[idx] != 0.0 || _mdl_yf.weight2[idx] != 0.0) {
	_mdl_yf.active[idx] = 1;
	_ncf++;
      }
    }
    // initialize nominal time
    _t_nominal = (ts(_KK) - ts(0)) / _KK;
  }

  // allocate states and controls for optimization
  int spsk, upsk;
  if (k < _K) {
    if (_multistage) {
      spsk = _sps;
      upsk = 1;
      x.alloc(_nx);
      u.alloc(_nu + spsk*_ns);
      if (k == 0)
	c.alloc(spsk*(_nc+_nsc) + _nc0);
      else
	c.alloc(spsk*(_nc+_nsc));
    }
    else {
      spsk = _KK;
      upsk = _KK;
      x.alloc(_nu, _nx);
      u.alloc(upsk*_nu + spsk*_ns);
      c.alloc(spsk*(_nc+_nsc) + upsk*_nu);
    }
  }
  else {
    spsk = 1;
    upsk = 1;
    // at final time allocate additional final states for slack variables
    // as no control parameters allowed here
    // (propagate us and xs to final stage as they are accessed in update)
    x.alloc(_nx + spsk*_ns);
    if (_K > 0)
      c.alloc(spsk*(_nc+_nsc) + _ncf);
    else
      c.alloc(spsk*(_nc+_nsc) + _nc0 + _ncf);
  }

  // setup initial states
  if (k == 0) {
    for (i = 0, idx = 0; idx < _mdl_nu; idx++) {
      if (_mdl_u_order[idx] < 0 || _mdl_u_order[idx] > 1)
	m_error(E_FORMAT, "Prg_SFunctionOpt::setup: "
		"mdl_u_order must be 0 or 1");
      if (_mdl_u.active[idx]) {
	x.initial[i] = _mdl_us[0][idx] / _mdl_u_nominal[idx];
	if (_mdl_u0_nfixed[idx] > 1 - _mdl_u_order[idx])
	  x.min[i] = x.max[i] = x.initial[i];
	else if (_mdl_u0_nfixed[idx] > 0 && _mdl_u_order[idx] == 0) {
	  // apply rate of change bounds to initial step for zero order hold
	  if (_mdl_der_u.min[idx] > -Inf)
	    x.min[i] = x.initial[i]
	      + _mdl_der_u.min[idx] / _mdl_u_nominal[idx]*_t_nominal;
	  if (_mdl_der_u.max[idx] < Inf)
	    x.max[i] = x.initial[i]
	      + _mdl_der_u.max[idx] / _mdl_u_nominal[idx]*_t_nominal;
	}
        i++;
      }
    }
    for (i = _nu; i < _nx; i++) {
      x.initial[i] = _mdl_xs[ks(0)][i-_nu] / _mdl_x_nominal[i-_nu];
      if (_mdl_x0_active[i-_nu]) {
	if (!_multistage)
	  m_error(E_FORMAT, "Prg_SFunctionOpt::setup: "
		  "mdl_x0_active=1 requires prg_multistage=true");
        // use the more restrictive bound of x0_min/max and x_min/max
        if (_mdl_x0.min[i-_nu] > _mdl_x.min[i-_nu])
          x.min[i] = _mdl_x0.min[i-_nu] / _mdl_x_nominal[i-_nu];
        else
          x.min[i] = _mdl_x.min[i-_nu] / _mdl_x_nominal[i-_nu];
        if (_mdl_x0.max[i-_nu] < _mdl_x.max[i-_nu])
          x.max[i] = _mdl_x0.max[i-_nu] / _mdl_x_nominal[i-_nu];
        else
          x.max[i] = _mdl_x.max[i-_nu] / _mdl_x_nominal[i-_nu];
      }
      else
	x.min[i] = x.max[i] = x.initial[i];
    }
  }
  // setup states of subsequent stages
  else {
    for (i = 0, idx = 0; idx < _mdl_nu; idx++)
      if (_mdl_u.active[idx]) {
	x.initial[i] = _mdl_us[ks(k)][idx] / _mdl_u_nominal[idx];
	i++;
      }
    for (i = _nu; i < _nx; i++) {
      x.initial[i] = _mdl_xs[ks(k)][i-_nu] / _mdl_x_nominal[i-_nu];
      x.min[i] = _mdl_x.min[i-_nu] / _mdl_x_nominal[i-_nu];
      x.max[i] = _mdl_x.max[i-_nu] / _mdl_x_nominal[i-_nu];
    }
  }

  // setup control inputs
  for (i = 0, idx = 0; idx < _mdl_nu; idx++) {
    if (_mdl_u.active[idx]) {
      if (_multistage && k >= _mdl_u0_nfixed[idx] + _mdl_u_order[idx] - 1
          || !_multistage && k == 0 && _mdl_u0_nfixed[idx] == 0) {
	// control bounds
        // at k=0, use the more restrictive bound of u0_min/max and u_min/max
        if (k == 0 && _mdl_u0.min[idx] > _mdl_u.min[idx])
          x.min[i] = _mdl_u0.min[idx] / _mdl_u_nominal[idx];
	else if (_mdl_u.min[idx] > -Inf)
	  x.min[i] = _mdl_u.min[idx] / _mdl_u_nominal[idx];
        if (k == 0 && _mdl_u0.max[idx] < _mdl_u.max[idx])
          x.max[i] = _mdl_u0.max[idx] / _mdl_u_nominal[idx];
	else if (_mdl_u.max[idx] < Inf)
	  x.max[i] = _mdl_u.max[idx] / _mdl_u_nominal[idx];
      }
      if (!_multistage && k < _K) {
	// treat control bounds via general constraints
	if (_mdl_u.min[idx] > -Inf) {
	  for (j = _mdl_u0_nfixed[idx] - _mdl_u_order[idx]; j < upsk; j++)
	    c.min[spsk*(_nc+_nsc) + i + j] =
	      _mdl_u.min[idx] / _mdl_u_nominal[idx];
	}
	if (_mdl_u.max[idx] < Inf) {
	  for (j = _mdl_u0_nfixed[idx] - _mdl_u_order[idx]; j < upsk; j++)
	    c.max[spsk*(_nc+_nsc) + i + j] =
	      _mdl_u.max[idx] / _mdl_u_nominal[idx];
	}
      }
      if (k < _K) {
	// rate of change bounds
	if (_mdl_der_u.min[idx] > -Inf) {
	  for (j = 0; j < upsk; j++)
	    u.min[i+j] =
	      _mdl_der_u.min[idx] / _mdl_u_nominal[idx]*_t_nominal;
	}
	if (_mdl_der_u.max[idx] < Inf) {
	  for (j = 0; j < upsk; j++)
	    u.max[i+j] =
	      _mdl_der_u.max[idx] / _mdl_u_nominal[idx]*_t_nominal;
	}
	// initial guess
	if (_multistage) {
	  u.initial[i] =
	    (_mdl_us[ks(k+1)][idx] - _mdl_us[ks(k)][idx])
	    / (ts(ks(k+1)) - ts(ks(k)))
	    / _mdl_u_nominal[idx]*_t_nominal;
	}
	else {
	  for (j = 0; j < upsk; j++) {
	    u.initial[i+j] =
	      (_mdl_us[j+1][idx] - _mdl_us[j][idx])
	      / (ts(j+1) - ts(j))
	      / _mdl_u_nominal[idx]*_t_nominal;
	  }
	}
	// override rate of change bounds for fixed controls
	if (_multistage) {
	  if (k < _mdl_u0_nfixed[idx] + _mdl_u_order[idx] - 2) {
	    u.min[i] = u.max[i] = u.initial[i];
	  }
	}
	else {
	  int jend = min(_mdl_u0_nfixed[idx] + _mdl_u_order[idx] - 2, upsk);
	  for (j = 0; j < jend; j++)
	    u.min[i+j] = u.max[i+j] = u.initial[i+j];
	}
	// override rate of change bounds for joined controls
	if (_multistage) {
	  if ((k + 1) % _mdl_u_decimation[idx] != 0) {
	    u.min[i] = u.max[i] = 0.0;
	  }
	}
	else {
	  for (j = 0; j < upsk; j++)
	    if ((j + 1) % _mdl_u_decimation[idx] != 0)
	      u.min[i+j] = u.max[i+j] = 0.0;
	}
      }
      // for zero order hold: constraint u parameter of last interval to zero
      // as control must not change at final time, compared to last interval
      if (k == _K-1 && _mdl_u_order[idx] == 0) {
	if (_multistage)
	  u.min[i] = u.max[i] = u.initial[i] = 0.0;
	else {
	  j = upsk - 1;
	  u.min[i+j] = u.max[i+j] = u.initial[i+j] = 0.0;
	}
      }
      i += upsk;
    }
  }

  // setup constraints on model outputs
  for (i = 0, isc = 0, idx = 0; idx < _mdl_ny; idx++) {
    if (_mdl_y.min[idx] > -Inf)
      for (j = 0; j < spsk; j++)
	c.min[i+j] = _mdl_y.min[idx] / _mdl_y_nominal[idx];
    if (_mdl_y.max[idx] < Inf)
      for (j = 0; j < spsk; j++)
	c.max[i+j] = _mdl_y.max[idx] / _mdl_y_nominal[idx];
    if (_mdl_y.active[idx])
      i += spsk;
    if (_mdl_y_soft.min[idx] > -Inf) {
      for (j = 0; j < spsk; j++)
	c.min[spsk*_nc + isc + j] = _mdl_y_soft.min[idx] / _mdl_y_nominal[idx];
      isc += spsk;
    }
    if (_mdl_y_soft.max[idx] < Inf) {
      for (j = 0; j < spsk; j++)
	c.max[spsk*_nc + isc + j] = _mdl_y_soft.max[idx] / _mdl_y_nominal[idx];
      isc += spsk;
    }
  }

  // setup slack variables for soft constraints
  if (k < _K)
    for (i = upsk*_nu; i < upsk*_nu + spsk*_ns; i++) {
      u.min[i] = 0.0;
      u.initial[i] = 0.0;
    }
  else
    for (i = _nx; i < _nx + spsk*_ns; i++) {
      x.min[i] = 0.0;
      x.initial[i] = 0.0;
    }

  // setup state and output constraints at initial time
  if (k == 0) {
    for (i = spsk*(_nc+_nsc), idx = 0; idx < _mdl_nx; idx++) {
      if (_mdl_x0_active[idx]
          && (_mdl_der_x0_min[idx] > -Inf || _mdl_der_x0_max[idx] < Inf)) {
        if (_mdl_der_x0_min[idx] > -Inf)
          c.min[i] = _mdl_der_x0_min[idx] / _mdl_x_nominal[idx];
        if (_mdl_der_x0_max[idx] < Inf)
          c.max[i] = _mdl_der_x0_max[idx] / _mdl_x_nominal[idx];
	i++;
      }
    }
    for (idx = 0; idx < _mdl_ny; idx++) {
      if (_mdl_y0.active[idx]) {
        if (_mdl_y0.min[idx] > -Inf)
          c.min[i] = _mdl_y0.min[idx] / _mdl_y_nominal[idx];
        if (_mdl_y0.max[idx] < Inf)
          c.max[i] = _mdl_y0.max[idx] / _mdl_y_nominal[idx];
	i++;
      }
    }
  }

  // setup output constraints at final time
  if (k == _K) {
    for (i = spsk*(_nc+_nsc), idx = 0; idx < _mdl_ny; idx++) {
      if (_mdl_yf.active[idx]) {
        if (_mdl_yf.min[idx] > -Inf)
          c.min[i] = _mdl_yf.min[idx] / _mdl_y_nominal[idx];
        if (_mdl_yf.max[idx] < Inf)
          c.max[i] = _mdl_yf.max[idx] / _mdl_y_nominal[idx];
        i++;
      }
    }
  }
}

//--------------------------------------------------------------------------
void Prg_SFunctionOpt::setup_struct(int k,
				    const Omu_VariableVec &x,
				    const Omu_VariableVec &u,
				    Omu_DependentVec &xt, Omu_DependentVec &F,
				    Omu_DependentVec &f,
				    Omu_Dependent &f0, Omu_DependentVec &c)
{
  int i, idx, j;

  // consistic just takes states from optimizer
  // note: possible changes due to discrete events (mdlUpdate) are neglected
  m_ident(xt.Jx);
  m_zero(xt.Ju);
  xt.set_linear();

  // explicit ODE for continuous-time equations
  m_ident(F.Jdx);
  sm_mlt(-1.0, F.Jdx, F.Jdx);
  F.set_linear(Omu_Dependent::WRT_dx);

  // setup struct of Jacobians for wrt. active inputs
  if (k < _K) {
    m_zero(F.Ju);
    for (i = 0, idx = 0; idx < _mdl_nu; idx++) {
      if (_mdl_u.active[idx]) {
	// no dependency on states for active inputs
	for (j = 0; j < _nx; j++)
	  F.Jx[i][j] = 0.0;
	if (_mdl_u_order[idx] == 0) {
	  // no dependecy on time derivative for zero order hold
	  F.Jdx[i][i] = 0.0;
	  // F.Ju is zero for zero order hold
	  for (j = 0; j < _nu; j++)
	    F.Ju[i][j] = 0.0;
	} else {
	  if (_multistage) {
	    // F.Ju is constant for multistage
	    for (j = 0; j < _nu; j++)
	      F.Ju[i][j] = 0.0;
	    F.Ju[i][i] = 1.0/_t_nominal;
	  }
	  else {
	    // just allocate dense structure
	    for (j = 0; j < _nu; j++)
	      F.Ju[i][j] = 1.0;
	  }
	}
	i++;
      }
    }
    if (_multistage)
      F.set_linear(Omu_Dependent::WRT_u);
  }

  // constraints c in update do not depend on xf
  m_zero(c.Jxf);
  c.set_linear(Omu_Dependent::WRT_xf);
  if (_multistage && _ns == 0) {
    m_zero(c.Ju);
    c.set_linear(Omu_Dependent::WRT_u);
  }

  // objective does not depend on xf
  v_zero(f0.gxf);
  f0.set_linear(Omu_Dependent::WRT_xf);
}

//--------------------------------------------------------------------------
void Prg_SFunctionOpt::init_simulation(int k,
				       Omu_VariableVec &x, Omu_VariableVec &u)
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
  int i, idx, is, isc;
  int spsk = kk < _KK? (_multistage? _sps: _KK): 1; // one sample at final time
  int upsk = _multistage? 1: _KK;

  // junction conditions for continuous-time equations
  if (kk < _KK) {
    // controlled inputs
    for (i = 0, idx = 0; idx < _mdl_nu; idx++) {
      if (_mdl_u.active[idx]) {
	if (_mdl_u_order[idx] == 0 && (!_multistage || (kk+1)%spsk == 0))
	  // zero order hold at end of stage:
	  // apply step in u for subsequent stage
	  f[i] = x[i] + ((ts(kk+1) - ts(kk-spsk/upsk+1)) * u[i*upsk + kk%upsk]
			 /_t_nominal);
	else
	  // piecewise linear interpolation or within stage with zoh
	  f[i] = xf[i];
	i++;
      }
    }
    assert(i == _nu); // problem structure must not have changed
    // states
    for (; i < _nx; i++)
      f[i] = xf[i];
    // re-use slack variables of last but one sample period at final time
    // note: this is only because no control parameters allowed at fineal time
    // and because states need to be defined with state equations
    if (kk == _KK-1) {
      for (; i < _nx + _ns; i++)
	f[i] = u[upsk*_nu + (i-_nx+1)*spsk-1];
    }
  }

  // set simulation time
  ssSetT(_S, ts(kk));

  // initialize model inputs
  real_T *mdl_u = NULL;
  if (ssGetNumInputPorts(_S) > 0) {
    if (ssGetInputPortRequiredContiguous(_S, 0))
      mdl_u = (real_T *)ssGetInputPortRealSignal(_S, 0);
    else
      mdl_u = (real_T *)*ssGetInputPortRealSignalPtrs(_S, 0);
  }
  for (i = 0, idx = 0; idx < _mdl_nu; idx++) {
    if (_mdl_u.active[idx])
      _mdl_us[kk][idx] = x[i++] * _mdl_u_nominal[idx];
    mdl_u[idx] = _mdl_us[kk][idx];
  }

  // pass current states to model
  real_T *mdl_x = ssGetContStates(_S);
  for (i = 0; i < _mdl_nx; i++)
    mdl_x[i] = x[_nu + i] * _mdl_x_nominal[i];

  // call mdlOutputs
  SMETHOD_CALL2(mdlOutputs, _S, 0); 

  // obtain model outputs
  real_T *mdl_y = NULL;
  if (ssGetNumOutputPorts(_S) > 0)
    mdl_y = ssGetOutputPortRealSignal(_S, 0);

  // correct model outputs with bias
  for (idx = 0; idx < _mdl_ny; idx++)
    mdl_y[idx] += _mdl_y_bias[idx];

  // utilize model outputs for constraints and optimization objective
  double help;
  f0 = 0.0;

  double dt; 	// factor for linear interpol. (trapezoidal integration rule)
  double dt0; 	// factor for zero order hold
  double dtu; 	// factor used for controlled inputs in objective function
  if (kk == 0) {
    if (_KK > 0) {
      dt = 0.5*(ts(kk+1) - ts(kk));
      dt0 = ts(kk+1) - ts(kk);
    } else {
      // steady-state problem
      dt = dt0 = 1.0;
    }
  } else if (kk == _KK) {
    dt = 0.5*(ts(kk) - ts(kk-1));
    dt0 = 0.0;
  } else {
    dt = 0.5*(ts(kk+1) - ts(kk-1));
    dt0 = ts(kk+1) - ts(kk);
  }

  // contribution of controlled inputs
  for (i = 0, idx = 0; idx < _mdl_nu; idx++) {
    if (_mdl_u.active[idx]) {
      // control inputs
      dtu = (_mdl_u_order[idx] == 0)? dt0: dt;
      f0 += dtu * _mdl_u.weight1[idx] * _mdl_us[kk][idx] / _mdl_u_nominal[idx];
      help = (_mdl_us[kk][idx] - _mdl_u.ref[idx]) / _mdl_u_nominal[idx];
      f0 += dtu * _mdl_u.weight2[idx]*help*help;
      // rates of change
      if (kk < _KK) {
	f0 += dt0 * _mdl_der_u.weight1[idx] * u[i*upsk + kk%upsk]/_t_nominal;
	help = u[i*upsk + kk%upsk]/_t_nominal - _mdl_der_u.ref[idx]
	  / _mdl_u_nominal[idx];
	f0 += dt0 * _mdl_der_u.weight2[idx]*help*help;
      }
      i++;
    }
  }
  // if not multistage, treat control bounds via general constraints
  // (note: do not access xf to avoid non-linear dependency in update)
  if (!_multistage && kk < _KK) {
    for (i = 0; i < _nu; i++)
      c[spsk*(_nc+_nsc) + i*spsk + kk%spsk] = x[i] +
	(ts(kk+1)-ts(kk)) * u[i*upsk + kk%upsk]/_t_nominal; // == xf[i]
  }
  // contribution of active outputs
  for (i = 0, is = 0, isc = 0, idx = 0; idx < _mdl_ny; idx++) {
    // assign used outputs to constraints
    if (_mdl_y.active[idx]) {
      c[i + kk%spsk] = mdl_y[idx] / _mdl_y_nominal[idx];
      i += spsk;
    }
    // calculate objective term
    if (_mdl_y.weight1[idx] != 0.0)
      f0 += dt * _mdl_y.weight1[idx] * mdl_y[idx] / _mdl_y_nominal[idx];
    if (_mdl_y.weight2[idx] != 0.0) {
      help = (mdl_y[idx] - _mdl_y.ref[idx]) / _mdl_y_nominal[idx];
      f0 += dt * _mdl_y.weight2[idx]*help*help;
    }
    // consider soft constraints
    if (_mdl_y_soft.active[idx] == 1) {
      if (kk < _KK)
	help = u[upsk*_nu + is + kk%spsk];
      else
	help = x[_nx + is + kk%spsk];

      f0 += dt * (_mdl_y_soft.weight1[idx]*help
		  + _mdl_y_soft.weight2[idx]*help*help);

      if (_mdl_y_soft.min[idx] > -Inf) {
	c[spsk*_nc + isc + kk%spsk] = mdl_y[idx] / _mdl_y_nominal[idx] + help;
	isc += spsk;
      }
      if (_mdl_y_soft.max[idx] < Inf) {
	c[spsk*_nc + isc + kk%spsk] = mdl_y[idx] / _mdl_y_nominal[idx] - help;
	isc += spsk;
      }
      is += spsk;
    }
  }

  // additional terms at initial time
  if (kk == 0) {
    i = spsk*(_nc+_nsc);
    // constraints on derivatives of initial states
    SMETHOD_CALL(mdlDerivatives, _S);
    real_T *mdl_dx = ssGetdX(_S);
    for (idx = 0; idx < _mdl_nx; idx++) {
      if (_mdl_x0_active[idx]
          && (_mdl_der_x0_min[idx] > -Inf || _mdl_der_x0_max[idx] < Inf)) {
        c[i] = mdl_dx[idx] / _mdl_x_nominal[idx];
        i++;
      }
    }
    for (idx = 0; idx < _mdl_ny; idx++) {
      // assign used outputs to constraints
      if (_mdl_y0.active[idx]) {
	c[i] = mdl_y[idx] / _mdl_y_nominal[idx];
	i++;
      }
      // initial objective terms
      if (_mdl_y0.weight1[idx] != 0.0)
	f0 += _mdl_y0.weight1[idx] * mdl_y[idx] / _mdl_y_nominal[idx];
      if (_mdl_y0.weight2[idx] != 0.0) {
	help = (mdl_y[idx] - _mdl_y.ref[idx]) / _mdl_y_nominal[idx];
	f0 += _mdl_y0.weight2[idx]*help*help;
      }
    }
  }

  // additional terms at final time
  if (kk == _KK) {
    for (i = _nc+_nsc, idx = 0; idx < _mdl_ny; idx++) {
      // assign used outputs to constraints
      if (_mdl_yf.active[idx]) {
	c[i] = mdl_y[idx] / _mdl_y_nominal[idx];
	i++;
      }
      // final objective terms
      if (_mdl_yf.weight1[idx] != 0.0)
	f0 += _mdl_yf.weight1[idx] * mdl_y[idx] / _mdl_y_nominal[idx];
      if (_mdl_yf.weight2[idx] != 0.0) {
	help = (mdl_y[idx] - _mdl_y.ref[idx]) / _mdl_y_nominal[idx];
	f0 += _mdl_yf.weight2[idx]*help*help;
      }
    }
  }

  // store values of model states
  if (kk == 0)
    for (i = 0; i < _mdl_nx; i++)
      _mdl_x0[i] = value(mdl_x[i]);
  for (i = 0; i < _mdl_nx; i++)
    _mdl_xs[kk][i] = value(mdl_x[i]);

  // store values of model outputs
  for (i = 0; i < _mdl_ny; i++)
    _mdl_ys[kk][i] = value(mdl_y[i]);

  // obtain Jacobians if required
  if (f.is_required_J() || f0.is_required_g() || c.is_required_J())
    update_grds(kk, x, u, xf, f, f0, c);
}

//--------------------------------------------------------------------------
void Prg_SFunctionOpt::update_grds(int kk, 
				   const Omu_StateVec &x, const Omu_Vec &u,
				   const Omu_StateVec &xf,
				   Omu_DependentVec &f, Omu_Dependent &f0,
				   Omu_DependentVec  &c)
{
  int i, idx, ii, iidx0, iidx, isc, iscdx, is, j, jdx, rdx;
  int spsk = kk < _KK? (_multistage? _sps: _KK): 1; // one sample at final time
  int upsk = _multistage? 1: _KK;

  if (ssGetmdlJacobian(_S) == NULL) {
    // call predefined update for numerical differentiation
    Omu_Program::update_grds(kk, x, u, xf, f, f0, c);
  }
  else {
    // exploit mdlJacobian to obtain c.Jx
    int mdl_u_idx, mdl_x_idx, mdl_dx_idx, mdl_y_idx;
    real_T *pr = ssGetJacobianPr(_S);
    int_T *ir = ssGetJacobianIr(_S);
    int_T *jc = ssGetJacobianJc(_S);

    mdlJacobian(_S);

    m_zero(c.Jx);
    m_zero(c.Ju);
    // Jacobian wrt S-function states (ddxdy/dx)
    for (j = _nu, jdx = 0; jdx < _mdl_nx; jdx++, j++) {
      mdl_x_idx = jdx;
      for (i = 0, idx = _mdl_nx, isc = spsk*_nc, iscdx = _mdl_nx,
	     ii = _nc+_nsc, iidx0 = 0, iidx = _mdl_nx,
	     rdx = jc[jdx]; rdx < jc[jdx+1]; rdx++) {
        if (ir[rdx] < _mdl_nx) {
          mdl_dx_idx = ir[rdx];
	  // constraints on derivatives for initial states
	  if (kk == 0 && _mdl_x0_active[mdl_dx_idx]
              && (_mdl_der_x0_min[mdl_dx_idx] > -Inf
                  || _mdl_der_x0_max[mdl_dx_idx] < Inf)) {
	    // need to loop to obtain ii considering active initial states
	    for (; iidx0 < mdl_dx_idx; iidx0++) {
	      if (_mdl_x0_active[iidx0]
                  && (_mdl_der_x0_min[iidx0] > -Inf
                      || _mdl_der_x0_max[iidx0] < Inf))
		ii++;
	    }
	    c.Jx[ii][j] = pr[rdx] /
	      _mdl_x_nominal[mdl_dx_idx] * _mdl_x_nominal[mdl_x_idx];
	  }
        }
        else {
	  mdl_y_idx = ir[rdx] - _mdl_nx;
	  if (_mdl_y.active[mdl_y_idx]) {
	    // need to loop through idx to obtain i considering active outputs
	    for (; idx < ir[rdx]; idx++) {
	      if (_mdl_y.active[idx - _mdl_nx])
		i += spsk;
	    }
	    c.Jx[i + kk%spsk][j] = pr[rdx] /
	      _mdl_y_nominal[mdl_y_idx] * _mdl_x_nominal[mdl_x_idx];
	  }
	  // soft constraints
	  if (_mdl_y_soft.active[mdl_y_idx] == 1) {
	    // need to loop to obtain isc considering active outputs
	    for (; iscdx < ir[rdx]; iscdx++) {
	      if (_mdl_y_soft.min[iscdx - _mdl_nx] > -Inf)
		isc += spsk;
	      if (_mdl_y_soft.max[iscdx - _mdl_nx] < Inf)
		isc += spsk;
	    }
	    if (_mdl_y_soft.min[mdl_y_idx] > -Inf) {
	      c.Jx[isc + kk%spsk][j] = pr[rdx] /
		_mdl_y_nominal[mdl_y_idx] * _mdl_x_nominal[mdl_x_idx];
	      isc += spsk;
	    }
	    if (_mdl_y_soft.max[mdl_y_idx] < Inf) {
	      c.Jx[isc + kk%spsk][j] = pr[rdx] /
		_mdl_y_nominal[mdl_y_idx] * _mdl_x_nominal[mdl_x_idx];
	    }
	    if (_mdl_y_soft.min[mdl_y_idx] > -Inf) {
	      isc -= spsk;
	    }
	  }
	  // constraints (used model outputs) at initial time
	  if (kk == 0 && _mdl_y0.active[mdl_y_idx]) {
	    // need to loop to obtain ii considering active outputs
	    for (; iidx < ir[rdx]; iidx++) {
	      if (_mdl_y0.active[iidx - _mdl_nx])
		ii++;
	    }
	    c.Jx[ii][j] = pr[rdx] /
	      _mdl_y_nominal[mdl_y_idx] * _mdl_x_nominal[mdl_x_idx];
	  }
	  // constraints (used model outputs) at final time
	  if (kk == _KK && _mdl_yf.active[mdl_y_idx]) {
	    // need to loop to obtain ii considering active outputs
	    for (; iidx < ir[rdx]; iidx++) {
	      if (_mdl_yf.active[iidx - _mdl_nx])
		ii++;
	    }
	    c.Jx[ii][j] = pr[rdx] /
	      _mdl_y_nominal[mdl_y_idx] * _mdl_x_nominal[mdl_x_idx];
	  }
	}
      }
    }
    // Jacobian wrt S-function inputs (ddxdy/du)
    // (note: these are states for the optimizer)
    for (j = 0, jdx = _mdl_nx; jdx < _mdl_nx + _mdl_nu; jdx++) {
      mdl_u_idx = jdx - _mdl_nx;
      if (_mdl_u.active[mdl_u_idx]) {
	for (i = 0, idx = _mdl_nx, isc = spsk*_nc, iscdx = _mdl_nx,
	       ii = _nc+_nsc, iidx0 = 0, iidx = _mdl_nx,
	       rdx = jc[jdx]; rdx < jc[jdx+1]; rdx++) {
          if (ir[rdx] < _mdl_nx) {
            mdl_dx_idx = ir[rdx];
            // constraints on derivatives for initial states
            if (kk == 0 && _mdl_x0_active[mdl_dx_idx]
                && (_mdl_der_x0_min[mdl_dx_idx] > -Inf
                    || _mdl_der_x0_max[mdl_dx_idx] < Inf)) {
              // need to loop to obtain ii considering active initial states
              for (; iidx0 < mdl_dx_idx; iidx0++) {
                if (_mdl_x0_active[iidx0]
                    && (_mdl_der_x0_min[iidx0] > -Inf
                        || _mdl_der_x0_max[iidx0] < Inf))
                  ii++;
              }
              c.Jx[ii][j] = pr[rdx] /
                _mdl_x_nominal[mdl_dx_idx] * _mdl_u_nominal[mdl_u_idx];
            }
          }
	  else {
	    mdl_y_idx = ir[rdx] - _mdl_nx;
	    if (_mdl_y.active[mdl_y_idx]) {
	      // need to loop to obtain i considering active outputs
	      for (; idx < ir[rdx]; idx++) {
		if (_mdl_y.active[idx - _mdl_nx])
		  i += spsk;
	      }
	      c.Jx[i + kk%spsk][j] = pr[rdx] /
		_mdl_y_nominal[mdl_y_idx] * _mdl_u_nominal[mdl_u_idx];
	    }
	    // soft constraints
	    if (_mdl_y_soft.active[mdl_y_idx] == 1) {
	      // need to loop to obtain isc considering active outputs
	      for (; iscdx < ir[rdx]; iscdx++) {
		if (_mdl_y_soft.min[iscdx - _mdl_nx] > -Inf)
		  isc += spsk;
		if (_mdl_y_soft.max[iscdx - _mdl_nx] < Inf)
		  isc += spsk;
	      }
	      if (_mdl_y_soft.min[mdl_y_idx] > -Inf) {
		c.Jx[isc + kk%spsk][j] = pr[rdx] /
		  _mdl_y_nominal[mdl_y_idx] * _mdl_u_nominal[mdl_u_idx];
		isc += spsk;
	      }
	      if (_mdl_y_soft.max[mdl_y_idx] < Inf) {
		c.Jx[isc + kk%spsk][j] = pr[rdx] /
		  _mdl_y_nominal[mdl_y_idx] * _mdl_u_nominal[mdl_u_idx];
	      }
	      if (_mdl_y_soft.min[mdl_y_idx] > -Inf) {
		isc -= spsk;
	      }
	    }
	    // constraints at initial time
	    if (kk == 0 && _mdl_y0.active[mdl_y_idx]) {
	      // need to loop to obtain ii considering active outputs
	      for (; iidx < ir[rdx]; iidx++) {
		if (_mdl_y0.active[iidx - _mdl_nx])
		  ii++;
	      }
	      c.Jx[ii][j] = pr[rdx] /
		_mdl_y_nominal[mdl_y_idx] * _mdl_u_nominal[mdl_u_idx];
	    }
	    // constraints at final time
	    if (kk == _KK && _mdl_yf.active[mdl_y_idx]) {
	      // need to loop to obtain ii considering active outputs
	      for (; iidx < ir[rdx]; iidx++) {
		if (_mdl_yf.active[iidx - _mdl_nx])
		  ii++;
	      }
	      c.Jx[ii][j] = pr[rdx] /
		_mdl_y_nominal[mdl_y_idx] * _mdl_u_nominal[mdl_u_idx];
	    }
	  }
	}
	j++;
      }
    }
    // contribution of slack variables for soft constraints
    for (is = 0, isc = spsk*_nc, idx = 0; idx < _mdl_ny; idx++) {
      if (_mdl_y_soft.active[idx] == 1) {
	if (_mdl_y_soft.min[idx] > -Inf) {
	  if (kk < _KK)
	    c.Ju[isc + kk%spsk][upsk*_nu + is + kk%spsk] += 1.0;
	  else
	    c.Jx[isc + kk%spsk][_nx + is + kk%spsk] += 1.0;
	  isc += spsk;
	}
	if (_mdl_y_soft.max[idx] < Inf) {
	  if (kk < _KK)
	    c.Ju[isc + kk%spsk][upsk*_nu + is + kk%spsk] -= 1.0;
	  else
	    c.Jx[isc + kk%spsk][_nx + is + kk%spsk] -= 1.0;
	  isc += spsk;
	}
	is += spsk;
      }
    }
    // control bounds if not multistage
    if (!_multistage && kk < _KK) {
      for (i = 0; i < _nu; i++) {
	c.Jx[spsk*(_nc+_nsc) + i*spsk + kk%spsk][i] = 1.0;
	c.Ju[spsk*(_nc+_nsc) + i*spsk + kk%spsk][i*upsk + kk%upsk] =
	  (ts(kk+1)-ts(kk))/_t_nominal;
      }
    }
  }

  // junction conditions for continuous-time equations
  if (kk < _KK) {
    // default values
    m_zero(f.Ju);
    m_zero(f.Jx);
    m_ident(f.Jxf);
    // modifications for controlled inputs with zero order hold
    for (i = 0, idx = 0; idx < _mdl_nu; idx++) {
      if (_mdl_u.active[idx]) {
	if (_mdl_u_order[idx] == 0 && (!_multistage || (kk+1)%spsk == 0)) {
	  // zero order hold at end of stage:
	  // apply step in u for subsequent stage
	  f.Jxf[i][i] = 0.0;
	  f.Jx[i][i] = 1.0;
	  f.Ju[i][i*upsk + kk%upsk] = ((ts(kk+1) - ts(kk-spsk/upsk+1))
				       /_t_nominal);
	}
	i++;
      }
    }
    if (kk == _KK-1) {
      // slack variables at final time
      for (i = _nx; i < _nx + _ns; i++)
	// note: f.Jxf has only _nx columns, but _nx+_ns rows
	// that is why don't set f.Jxf[i][i] = 0.0;
	f.Ju[i][upsk*_nu + (i-_nx+1)*spsk-1] = 1.0;
    }
  }

  // apply chain rule to calculate gradient of f0 
  // (do this as well if Omu_Program::update_grds was used, as complete
  //  numeric differentiation gives bad results for quadratic terms)
  double dt; 	// factor for linear interpol. (trapezoidal integration rule)
  double dt0; 	// factor for zero order hold
  double dtu; 	// factor used for controlled inputs in objective function
  if (kk == 0) {
    if (_KK > 0) {
      dt = 0.5*(ts(kk+1) - ts(kk));
      dt0 = ts(kk+1) - ts(kk);
    } else {
      // steady-state problem
      dt = dt0 = 1.0;
    }
  } else if (kk == _KK) {
    dt = 0.5*(ts(kk) - ts(kk-1));
    dt0 = 0.0;
  } else {
    dt = 0.5*(ts(kk+1) - ts(kk-1));
    dt0 = ts(kk+1) - ts(kk);
  }

  v_zero(f0.gx);
  v_zero(f0.gu);
  // controlled inputs
  for (i = 0, idx = 0; idx < _mdl_nu; idx++) {
    // contribution of objective term
    if (_mdl_u.active[idx]) {
      // controlled inputs
      dtu = (_mdl_u_order[idx] == 0)? dt0: dt;
      f0.gx[i] += dtu * _mdl_u.weight1[idx];
      f0.gx[i] += dtu * _mdl_u.weight2[idx] *
	2.0 * (x[i] - _mdl_u.ref[idx]/_mdl_u_nominal[idx]);
      // rates of change
      if (kk < _KK) {
	f0.gu[i] += dt0 * _mdl_der_u.weight1[idx]/_t_nominal;
	f0.gu[i] += dt0 * _mdl_der_u.weight2[idx]/_t_nominal *
	  2.0 * (u[i]/_t_nominal - _mdl_der_u.ref[idx]/_mdl_u_nominal[idx]);
      }
      i++;
    }
  }
  // active outputs
  for (i = 0, is = 0, idx = 0; idx < _mdl_ny; idx++) {
    // contribution of objective term
    if (_mdl_y.weight1[idx] != 0.0) {
      for (j = 0; j < _nx; j++)
	f0.gx[j] += dt * _mdl_y.weight1[idx] * c.Jx[i + kk%spsk][j];
    }
    if (_mdl_y.weight2[idx] != 0.0) {
      for (j = 0; j < _nx; j++)
	f0.gx[j] += dt * _mdl_y.weight2[idx] *
	  2.0 * (c[i + kk%spsk] - _mdl_y.ref[idx]/_mdl_y_nominal[idx]) *
	  c.Jx[i + kk%spsk][j];
    }
    if (_mdl_y.active[idx])
      i += spsk;

    // contribution of soft constraints
    if (_mdl_y_soft.active[idx] == 1) {
      if (kk < _KK)
	f0.gu[upsk*_nu + is + kk%spsk]
	  += dt * (_mdl_y_soft.weight1[idx]
		   + 2.0*_mdl_y_soft.weight2[idx]*u[upsk*_nu + is + kk%spsk]);
      else
	f0.gx[_nx + is + kk%spsk]
	  += dt * (_mdl_y_soft.weight1[idx]
		   + 2.0*_mdl_y_soft.weight2[idx]*x[_nx + is + kk%spsk]);
      is += spsk;
    }
  }
  // additional terms at initial time
  if (kk == 0) {
    for (i = spsk*(_nc+_nsc), idx = 0; idx < _mdl_ny; idx++) {
      // contributions of final objective terms
      if (_mdl_y0.weight1[idx] != 0.0) {
	for (j = 0; j < _nx; j++)
	  f0.gx[j] += _mdl_y0.weight1[idx] * c.Jx[i][j];
      }
      if (_mdl_y0.weight2[idx] != 0.0) {
	for (j = 0; j < _nx; j++)
	  f0.gx[j] += _mdl_y0.weight2[idx] *
	    2.0 * (c[i] - _mdl_y.ref[idx]/_mdl_y_nominal[idx]) *
	    c.Jx[i][j];
      }
      if (_mdl_y0.active[idx])
	i++;
    }
  }
  // additional terms at final time
  if (kk == _KK) {
    for (i = _nc+_nsc, idx = 0; idx < _mdl_ny; idx++) {
      // contributions of final objective terms
      if (_mdl_yf.weight1[idx] != 0.0) {
	for (j = 0; j < _nx; j++)
	  f0.gx[j] += _mdl_yf.weight1[idx] * c.Jx[i][j];
      }
      if (_mdl_yf.weight2[idx] != 0.0) {
	for (j = 0; j < _nx; j++)
	  f0.gx[j] += _mdl_yf.weight2[idx] *
	    2.0 * (c[i] - _mdl_y.ref[idx]/_mdl_y_nominal[idx]) *
	    c.Jx[i][j];
      }
      if (_mdl_yf.active[idx])
	i++;
    }
  }
}

//--------------------------------------------------------------------------
void Prg_SFunctionOpt::consistic(int kk, double t,
				 const Omu_StateVec &x, const Omu_Vec &u,
				 Omu_DependentVec &xt)
{
  int i, idx;

  // set simulation time
  ssSetT(_S, t);

  // initialize model inputs
  real_T *mdl_u = NULL;
  if (ssGetNumInputPorts(_S) > 0) {
    if (ssGetInputPortRequiredContiguous(_S, 0))
      mdl_u = (real_T *)ssGetInputPortRealSignal(_S, 0);
    else
      mdl_u = (real_T *)*ssGetInputPortRealSignalPtrs(_S, 0);
  }
  for (i = 0, idx = 0; idx < _mdl_nu; idx++) {
    if (_mdl_u.active[idx])
      mdl_u[idx] = x[i++] * _mdl_u_nominal[idx];
    else
      mdl_u[idx] = _mdl_us[kk][idx];
  }

  // initialize model in first stage
  if (kk == 0 && ssGetmdlInitializeConditions(_S) != NULL) {
    // initialize model
    SMETHOD_CALL(mdlInitializeConditions, _S);
  }

  // pass states from optimizer to model
  real_T *mdl_x = ssGetContStates(_S);
  for (i = 0; i < _mdl_nx; i++)
    mdl_x[i] = x[_nu + i] * _mdl_x_nominal[i];

  // call mdlUpdate to get discrete events processed
  // Note: this is done once at the beginning of a sample interval;
  // no event processing takes place during the integration.
  // mdlUpdate is not called at initial time to prevent initialization
  // of potentially optimized initial states, e.g. from parameters;
  // it is called once in Prg_SFunction::setup_model() instead.
  //if (kk > 0 && ssGetmdlUpdate(_S) != NULL) {
  if (ssGetmdlUpdate(_S) != NULL) {
    // also call mdlOutputs as done by Simulink before each mdlUpdate
    SMETHOD_CALL2(mdlOutputs, _S, 0); 
    SMETHOD_CALL2(mdlUpdate, _S, 0);
  }

  // take over optimized control inputs from optimizer
  for (i = 0; i < _nu; i++)
    xt[i] = x[i];
  // read back states from model
  // Note: take optimized initial states from the optimizer
  // to prevent their overriding by the model, e.g. from parameters;
  // they are read once in Prg_SFunction::setup_model() instead.
  for (i = 0; i < _mdl_nx; i++) {
    if (kk == 0 && _mdl_x0_active[i])
      xt[_nu + i] = x[_nu + i];
    else
      xt[_nu + i] = mdl_x[i] / _mdl_x_nominal[i];
  }
}

//--------------------------------------------------------------------------
void Prg_SFunctionOpt::continuous(int kk, double t,
				  const Omu_StateVec &x, const Omu_Vec &u,
				  const Omu_StateVec &dx, Omu_DependentVec &F)
{
  int i, idx;
  int upsk = _multistage? 1: _KK;

  // set simulation time
  ssSetT(_S, t);

  // initialize model inputs
  real_T *mdl_u = NULL;
  if (ssGetNumInputPorts(_S) > 0) {
    if (ssGetInputPortRequiredContiguous(_S, 0))
      mdl_u = (real_T *)ssGetInputPortRealSignal(_S, 0);
    else
      mdl_u = (real_T *)*ssGetInputPortRealSignalPtrs(_S, 0);
  }
  for (i = 0, idx = 0; idx < _mdl_nu; idx++) {
    if (_mdl_u.active[idx])
      mdl_u[idx] = x[i++] * _mdl_u_nominal[idx];
    else
      mdl_u[idx] = _mdl_us[kk][idx];
  }

  // pass current states to model
  real_T *mdl_x = ssGetContStates(_S);
  for (i = 0; i < _mdl_nx; i++)
    mdl_x[i] = x[_nu + i] * _mdl_x_nominal[i];

  // set model outputs before calculating derivatives
  // (note: this is required for sub-blocks providing inputs other blocks)
  // (furthermore a model may check for discrete changes in mdlOutputs)
  SMETHOD_CALL2(mdlOutputs, _S, 0); 

  // evaluate continuous model equations
  SMETHOD_CALL(mdlDerivatives, _S);

  // get model derivatives and change to residual form
  real_T *mdl_dx = ssGetdX(_S);
  for (i = 0; i < _mdl_nx; i++)
    F[_nu + i] = mdl_dx[i]/_mdl_x_nominal[i] - dx[_nu + i];

  // model equations for controlled inputs
  for (i = 0, idx = 0; idx < _mdl_nu; idx++) {
    if (_mdl_u.active[idx]) {
      if (_mdl_u_order[idx] == 0)
	// zero order hold
	F[i] = 0.0 - dx[i];
      else
	// piecewise linear interpolation
	F[i] = u[i*upsk + kk%upsk]/_t_nominal - dx[i];
      i++;
    }
  }
  assert(i == _nu); // problem structure must not have changed

  // obtain Jacobians if required
  if (F.is_required_J())
    continuous_grds(kk, t, x, u, dx, F);
}

//--------------------------------------------------------------------------
void Prg_SFunctionOpt::continuous_grds(int kk, double t,
				       const Omu_StateVec &x,
				       const Omu_Vec &u,
				       const Omu_StateVec &dx,
				       Omu_DependentVec &F)
{
  if (ssGetmdlJacobian(_S) == NULL) {
    // call predefined continuous_grds for numerical differentiation
    Omu_Program::continuous_grds(kk, t, x, u, dx, F);
  }
  else {
    // exploit mdlJacobian
    int i, idx, j, jdx, rdx;
    int mdl_u_idx, mdl_x_idx;
    real_T *pr = ssGetJacobianPr(_S);
    int_T *ir = ssGetJacobianIr(_S);
    int_T *jc = ssGetJacobianJc(_S);

    mdlJacobian(_S);

    // obtain Jx
    m_zero(F.Jx);
    // Jacobian wrt S-function states (ddx/dx)
    for (jdx = 0; jdx < _mdl_nx; jdx++) {
      mdl_x_idx = jdx;
      for (rdx = jc[jdx]; rdx < jc[jdx+1]; rdx++) {
	idx = ir[rdx];
	if (idx >= _mdl_nx)
	  break;
	i = _nu + idx;
	j = _nu + jdx;
	F.Jx[i][j] = pr[rdx] /
	  _mdl_x_nominal[idx] * _mdl_x_nominal[mdl_x_idx];
      }
    }
    // Jacobian wrt S-function inputs (ddx/du)
    // (note: these are states for the optimizer)
    for (j = 0, jdx = _mdl_nx; jdx < _mdl_nx + _mdl_nu; jdx++) {
      mdl_u_idx = jdx - _mdl_nx;
      if (_mdl_u.active[mdl_u_idx]) {
	for (rdx = jc[jdx]; rdx < jc[jdx+1]; rdx++) {
	  idx = ir[rdx];
	  if (idx >= _mdl_nx)
	    break;
	  i = _nu + idx;
	  F.Jx[i][j] = pr[rdx] /
	    _mdl_x_nominal[i-_nu] * _mdl_u_nominal[mdl_u_idx];
	}
	j++;
      }
    }

    // set F.Ju if not multistage
    // note: this can't be done in setup_struct as F.Ju changes within stage
    if (!_multistage) {
      int upsk = _KK;
      m_zero(F.Ju);
      for (i = 0, idx = 0; idx < _mdl_nu; idx++) {
	if (_mdl_u.active[idx]) {
	  if (_mdl_u_order[idx] > 0)
	    // piecewise linear interpolation
	    F.Ju[i][i*upsk + kk%upsk] = 1.0/_t_nominal;
	  i++;
	}
      }
    }
  }
}


//==========================================================================

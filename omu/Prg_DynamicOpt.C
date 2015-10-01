/*
 * Prg_DynamicOpt.C -- class definition
 *
 */

/*
    Copyright (C) 1997--2014  Ruediger Franke

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

#include "Prg_DynamicOpt.h"

#include <stdlib.h>

#include <If_Int.h>
#include <If_Real.h>
#include <If_RealVec.h>
#include <If_RealMat.h>
#include <If_IntVec.h>
#include <If_Int.h>
#include <If_IntVec.h>
#include <If_Method.h>

// redefine assert to throw an error instead of aborting
#undef assert
#define assert(expr) if (!(expr)) m_error(E_INTERN, "assert(" #expr ")");

#define GET_SET_CB(vartype, prefix, name) \
  GET_CB(vartype, prefix, name), \
  IF_SET_CB(vartype, Prg_DynamicOpt, set_##name)

#define GET_CB(vartype, prefix, name) \
  prefix#name, \
  IF_GET_CB(vartype, Prg_DynamicOpt, name)

// Throw E_CONV as errors occuring during solution are generally
// due to bad values and need to be treated (esp. during stepsize check).
#undef SMETHOD_ERROR
#define SMETHOD_ERROR E_CONV

typedef If_Method<Prg_DynamicOpt> If_Cmd;

IF_CLASS_DEFINE("DynamicOpt", Prg_DynamicOpt, Omu_Program);
IF_CLASS_DEFINE("SFunctionOpt", Prg_SFunctionOpt, Omu_Program);

//--------------------------------------------------------------------------
Prg_DynamicOpt::Prg_DynamicOpt()
{
  _sps = 1;
  _multistage = 1;
  _t_scale_idx = -1;
  _t_active = 0;
  _t_scale_i = -1;
  _t_scale_nominal = 1.0;

  _mdl_x0_active = iv_get(_mdl_nx);
  _mdl_der_x0_min = v_get(_mdl_nx);
  _mdl_der_x0_max = v_get(_mdl_nx);
  _mdl_u_order = iv_get(_mdl_nu);
  _mdl_u0_nfixed = iv_get(_mdl_nu);
  _mdl_u_decimation = iv_get(_mdl_nu);
  _mdl_u_periodic = iv_get(_mdl_nu);
  _mdl_x_periodic = iv_get(_mdl_nx);
  _t_nominal = 1.0;
  _mdl_y_order = iv_get(_mdl_ny);
  _mdl_y_bias = v_get(_mdl_ny);
  v_set(_mdl_x0, 0.0);
  iv_set(_mdl_x0_active, 0);
  v_set(_mdl_der_x0_min, -Inf);
  v_set(_mdl_der_x0_max, Inf);
  iv_set(_mdl_u_order, 1);
  iv_set(_mdl_u0_nfixed, 0);
  iv_set(_mdl_u_decimation, 1);
  iv_set(_mdl_u_periodic, 0);
  iv_set(_mdl_x_periodic, 0);
  iv_set(_mdl_y_order, 1);
  v_set(_mdl_y_bias, 0.0);

  // numbers of optimization variables
  _nu = 0;
  _nx = _nu + _mdl_nx;
  _nc = 0;
  _ns = 0;
  _nsc = 0;
  _nsu = 0;
  _nsuc = 0;
  _nc0 = 0;
  _ncf = 0;
  _nsf = 0;
  _nscf = 0;

  _mdl_us = m_get(_KK+1, _mdl_nu);
  _mdl_xs = m_get(_KK+1, _mdl_nx);
  _mdl_ys = m_get(_KK+1, _mdl_ny);
  _taus = v_get(_KK+1);

  _ifList.append(new If_Cmd("prg_setup_model",
			    &Prg_DynamicOpt::setup_model, this));
  _ifList.append(new If_Int("prg_sps", &_sps));
  _ifList.append(new If_Int(GET_SET_CB(int, "prg_", multistage)));
  _ifList.append(new If_Int(GET_SET_CB(int, "", mdl_t_scale_idx)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "prg_", taus)));

  // redefine mdl_x0 in order to also consider mdl_xs
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_x0)));
  _ifList.append(new If_IntVec(GET_SET_CB(const IVECP, "", mdl_x0_active)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_x0_min)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_x0_max)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_der_x0_min)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_der_x0_max)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_u0)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_u0_min)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_u0_max)));
  _ifList.append(new If_RealVec(GET_CB(const VECP, "", mdl_y0)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_y0_min)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_y0_max)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_y0_weight1)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_y0_weight2)));

  _ifList.append(new If_IntVec(GET_SET_CB(const IVECP, "", mdl_u_order)));
  _ifList.append(new If_IntVec(GET_SET_CB(const IVECP, "", mdl_u_active)));
  _ifList.append(new If_IntVec(GET_SET_CB(const IVECP, "", mdl_u_integer)));
  _ifList.append(new If_IntVec(GET_SET_CB(const IVECP, "", mdl_u0_nfixed)));
  _ifList.append(new If_IntVec(GET_SET_CB(const IVECP, "", mdl_u_decimation)));
  _ifList.append(new If_IntVec(GET_SET_CB(const IVECP, "", mdl_u_periodic)));
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

  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "",
                                           mdl_der_u_soft_min)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "",
                                           mdl_der_u_soft_max)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "",
					   mdl_der_u_soft_weight1)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "",
					   mdl_der_u_soft_weight2)));

  _ifList.append(new If_IntVec(GET_SET_CB(const IVECP, "", mdl_x_integer)));
  _ifList.append(new If_IntVec(GET_SET_CB(const IVECP, "", mdl_x_periodic)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_x_min)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_x_max)));

  _ifList.append(new If_IntVec(GET_SET_CB(const IVECP, "", mdl_y_order)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_y_bias)));
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

  _ifList.append(new If_RealVec(GET_CB(const VECP, "", mdl_uf)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_uf_min)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_uf_max)));
  _ifList.append(new If_RealVec(GET_CB(const VECP, "", mdl_yf)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_yf_min)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_yf_max)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_yf_weight1)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_yf_weight2)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_yf_soft_min)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_yf_soft_max)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_yf_soft_weight1)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_yf_soft_weight2)));

  _ifList.append(new If_RealMat(GET_SET_CB(const MATP, "", mdl_us)));
  _ifList.append(new If_RealMat(GET_SET_CB(const MATP, "", mdl_xs)));
  _ifList.append(new If_RealMat(GET_CB(const MATP, "", mdl_ys)));
}

//--------------------------------------------------------------------------
Prg_DynamicOpt::~Prg_DynamicOpt()
{
  v_free(_taus);
  m_free(_mdl_ys);
  m_free(_mdl_xs);
  m_free(_mdl_us);
  v_free(_mdl_y_bias);
  iv_free(_mdl_y_order);
  iv_free(_mdl_x_periodic);
  iv_free(_mdl_u_periodic);
  iv_free(_mdl_u_decimation);
  iv_free(_mdl_u0_nfixed);
  iv_free(_mdl_u_order);
  iv_free(_mdl_x0_active);
  v_free(_mdl_der_x0_max);
  v_free(_mdl_der_x0_min);
}

//--------------------------------------------------------------------------
void Prg_DynamicOpt::setup_model()
{
  // load S-function
  Omu_Model::setup_model(_t0);

  // check for optional S-function methods that are required
  if (_mdl_nx > _mdl_nd)
    assert(ssGetmdlDerivatives(_SS) != NULL);

  // adapt sizes of model vectors
  _mdl_x0.alloc(_mdl_nx);
  _mdl_u0.alloc(_mdl_nu);
  _mdl_y0.resize(_mdl_ny);
  _mdl_u.resize(_mdl_nu);
  _mdl_der_u.resize(_mdl_nu);
  _mdl_der_u_soft.resize(_mdl_nu);
  _mdl_x.alloc(_mdl_nx);
  _mdl_y.resize(_mdl_ny);
  _mdl_y_soft.resize(_mdl_ny);
  _mdl_uf.resize(_mdl_nu);
  _mdl_yf.resize(_mdl_ny);
  _mdl_yf_soft.resize(_mdl_ny);

  iv_resize(_mdl_x0_active, _mdl_nx);
  v_resize(_mdl_der_x0_min, _mdl_nx);
  v_resize(_mdl_der_x0_max, _mdl_nx);
  iv_resize(_mdl_u_order, _mdl_nu);
  iv_resize(_mdl_u0_nfixed, _mdl_nu);
  iv_resize(_mdl_u_decimation, _mdl_nu);
  iv_resize(_mdl_u_periodic, _mdl_nu);
  iv_resize(_mdl_x_periodic, _mdl_nx);
  iv_resize(_mdl_y_order, _mdl_ny);
  v_resize(_mdl_y_bias, _mdl_ny);
  iv_set(_mdl_x0_active, 0);
  v_set(_mdl_der_x0_min, -Inf);
  v_set(_mdl_der_x0_max, Inf);
  iv_set(_mdl_u_order, 1);
  iv_set(_mdl_u0_nfixed, 0);
  iv_set(_mdl_u_decimation, 1);
  iv_set(_mdl_u_periodic, 0);
  iv_set(_mdl_x_periodic, 0);
  iv_set(_mdl_y_order, 1);
  v_set(_mdl_y_bias, 0.0);
}

//--------------------------------------------------------------------------
void Prg_DynamicOpt::setup_stages(IVECP ks, VECP ts)
{
  int kk, j;

  // setup S-function
  if (_mdl_needs_setup)
    setup_model();

  // setup optimization problem
  if (_sps < 1) {
    m_error(E_FORMAT, "Prg_DynamicOpt::setup_stages: "
	    "prg_sps must at least be one");
  }
  if (_multistage) {
    _K = _KK/_sps;
    stages_alloc(ks, ts, _K, _sps);
  }
  else {
    if (_sps > 1) {
      m_error(E_FORMAT, "Prg_DynamicOpt::setup_stages: "
	      "prg_sps>1 requires prg_multistage=true");
    }
    _K = 1;
    stages_alloc(ks, ts, 1, _KK);
  }

  m_resize(_mdl_us, _KK+1, _mdl_nu);
  m_resize(_mdl_xs, _KK+1, _mdl_nx);
  m_resize(_mdl_ys, _KK+1, _mdl_ny);

  // setup _mdl_xs and mdl_us with start values from model
  set_mdl_x0(Omu_Model::_mdl_x_start);
  set_mdl_u0(Omu_Model::_mdl_u_start);

  // setup scaling of time
  v_copy(Omu_Program::ts(), _taus);
  if (!_mdl_is_fmu && _t_scale_idx >= 0) {
    // skip initialization for FMU because its start value is used
    for (kk = 0; kk < _KK; kk++)
      _mdl_us[kk][_t_scale_idx] = 1.0;
  }
}

//--------------------------------------------------------------------------
void Prg_DynamicOpt::setup(int k,
			   Omu_VariableVec &x, Omu_VariableVec &u,
			   Omu_VariableVec &c)
{
  int i, idx, isc, iscf, j;

  // complete general setup
  if (k == 0) {
    // take over possibly modified _mdl_p
    write_mx_args(_mdl_p);

    // obtain numbers of optimization variables (active model variables)
    _nu = 0;
    _nsu = 0;
    _nsuc = 0;
    for (idx = 0; idx < _mdl_nu; idx++) {
      if (idx == _t_scale_idx)
        _t_scale_i = _nu;
      if (_mdl_u.active[idx])
	_nu++;
      // count soft constraints for u, i.e. rate-of-change
      if ((_mdl_der_u_soft.min[idx] > -Inf || _mdl_der_u_soft.max[idx] < Inf) 
	  && (_mdl_der_u_soft.weight1[idx] != 0.0
              || _mdl_der_u_soft.weight2[idx] != 0.0)
          && _mdl_u.active[idx]) {
	_mdl_der_u_soft.active[idx] = 1;
	_nsu++;
        if (_mdl_der_u_soft.min[idx] > -Inf) {
          _nsuc++;
        }
        if (_mdl_der_u_soft.max[idx] < Inf) {
          _nsuc++;
        }
      }
      else
        _mdl_der_u_soft.active[idx] = 0;
    }
    // initialize time scaling
    if (_t_scale_idx >= 0) {
      // take over specification of model variable
      _t_active = _mdl_u.active[_t_scale_idx];
      _t_scale_nominal = _mdl_u_nominal[_t_scale_idx];
      if (_mdl_u_order[_t_scale_idx] != 0)
	  m_error(E_FORMAT, "Prg_DynamicOpt::setup: "
		  "mdl_u_order[mdl_t_scale_idx] must be 0 (zero order hold)");
    }
    else
      _t_active = 0;

    _nx = _nu + _mdl_nx;
    _nc = 0;
    _ns = 0;
    _nsc = 0;
    _nc0 = 0;
    _ncf = 0;
    _nsf = 0;
    _nscf = 0;
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
      else
        _mdl_y.active[idx] = 0;
      if ((_mdl_y_soft.min[idx] > -Inf || _mdl_y_soft.max[idx] < Inf) 
	  && (_mdl_y_soft.weight1[idx] != 0.0 || _mdl_y_soft.weight2[idx] != 0.0)) {
	_mdl_y_soft.active[idx] = 1;
	_ns++;
        if (_mdl_y_soft.min[idx] > -Inf) {
          _nsc++;
        }
        if (_mdl_y_soft.max[idx] < Inf) {
          _nsc++;
        }
      }
      else
        _mdl_y_soft.active[idx] = 0;
      if (_mdl_y0.min[idx] > -Inf || _mdl_y0.max[idx] < Inf
	  || _mdl_y0.weight1[idx] != 0.0 || _mdl_y0.weight2[idx] != 0.0) {
	if (!_multistage)
	  m_error(E_FORMAT, "Prg_DynamicOpt::setup: "
		  "use of mdl_y0 requires prg_multistage=true");
	_mdl_y0.active[idx] = 1;
	_nc0++;
      }
      else
        _mdl_y0.active[idx] = 0;
      // Note: constraints at final time require _K > 0
      //       as not both: initial and final constraints are treated for _K == 0
      if (_K > 0 && (_mdl_yf.min[idx] > -Inf || _mdl_yf.max[idx] < Inf
                     || _mdl_yf.weight1[idx] != 0.0 || _mdl_yf.weight2[idx] != 0.0)) {
	_mdl_yf.active[idx] = 1;
	_ncf++;
      }
      else
        _mdl_yf.active[idx] = 0;
      if (_K > 0 && (_mdl_yf_soft.min[idx] > -Inf || _mdl_yf_soft.max[idx] < Inf) 
	  && (_mdl_yf_soft.weight1[idx] != 0.0 || _mdl_yf_soft.weight2[idx] != 0.0)) {
	_mdl_yf_soft.active[idx] = 1;
	_nsf++;
        if (_mdl_yf_soft.min[idx] > -Inf) {
          _nscf++;
        }
        if (_mdl_yf_soft.max[idx] < Inf) {
          _nscf++;
        }
      }
      else
        _mdl_yf_soft.active[idx] = 0;
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
      if (k == _K-1)
        // Note: allocate additional us for slacks at final time
        //       as the final states must have state equations
        u.alloc(_nu + _nsu + spsk*_ns + _nsf);
      else
        u.alloc(_nu + _nsu + spsk*_ns);
      if (k == 0)
	c.alloc(spsk*(_nc+_nsc) + _nsuc + _nc0);
      else
	c.alloc(spsk*(_nc+_nsc) + _nsuc);
    }
    else {
      spsk = _KK;
      upsk = _KK;
      x.alloc(_nu, _nx);
      u.alloc(upsk*(_nu+_nsu) + spsk*_ns);
      c.alloc(spsk*(_nc+_nsc) + upsk*(_nsuc+_nu));
    }
  }
  else {
    spsk = 1;
    upsk = 1;
    // at final time allocate additional final states for slack variables
    // as no control parameters allowed here
    // (propagate us and xs to final stage as they are accessed in update)
    x.alloc(_nx + spsk*_ns + _nsf);
    if (_K > 0)
      c.alloc(spsk*(_nc+_nsc) + _ncf + _nscf);
    else
      c.alloc(spsk*(_nc+_nsc) + _nc0 + _ncf + _nscf);
  }

  // setup states
  for (i = _mdl_nd, idx = 0; idx < _mdl_nu; idx++) {
    if (k == 0 && _mdl_u_order[idx] < 0 || _mdl_u_order[idx] > 1)
      m_error(E_FORMAT, "Prg_DynamicOpt::setup: "
              "mdl_u_order must be 0 or 1");
    if (_mdl_u.active[idx]) {
      x.initial[i] = _mdl_us[ks(k)][idx] / _mdl_u_nominal[idx];
      if (k == 0 && _mdl_u0_nfixed[idx] > 1 - _mdl_u_order[idx])
        // fixed initial value
        x.min[i] = x.max[i] = x.initial[i]; 
      else if (k == 0 && _mdl_u0_nfixed[idx] > 0 && _mdl_u_order[idx] == 0) {
        // apply rate of change bounds to initial step for zero order hold
        if (_mdl_der_u.min[idx] > -Inf)
          x.min[i] = x.initial[i]
            + _mdl_der_u.min[idx] / _mdl_u_nominal[idx]*_t_nominal;
        if (_mdl_der_u.max[idx] < Inf)
          x.max[i] = x.initial[i]
            + _mdl_der_u.max[idx] / _mdl_u_nominal[idx]*_t_nominal;
      }
      else if (_multistage && k >= _mdl_u0_nfixed[idx] + _mdl_u_order[idx] - 1
               || !_multistage && k == 0 && _mdl_u0_nfixed[idx] == 0) {
	// control bounds
        if (k==0 && _mdl_u_periodic[idx])
          // indicate periodic state to Hqp_Docp
          x.min[i] = Inf;
        else if (k == 0 && _mdl_u0.min[idx] > _mdl_u.min[idx])
          // at k==0, use more restrictive bound of u0_min/max and u_min/max
          x.min[i] = _mdl_u0.min[idx] / _mdl_u_nominal[idx];
        else if (k == _K && _mdl_uf.min[idx] > _mdl_u.min[idx])
          // at k==_K, use more restrictive bound of uf_min/max and u_min/max
          x.min[i] = _mdl_uf.min[idx] / _mdl_u_nominal[idx];
	else if (_mdl_u.min[idx] > -Inf)
	  x.min[i] = _mdl_u.min[idx] / _mdl_u_nominal[idx];
        if (k==0 && _mdl_u_periodic[idx])
          x.max[i] = Inf;
        else if (k == 0 && _mdl_u0.max[idx] < _mdl_u.max[idx])
          x.max[i] = _mdl_u0.max[idx] / _mdl_u_nominal[idx];
        else if (k == _K && _mdl_uf.max[idx] < _mdl_u.max[idx])
          x.max[i] = _mdl_uf.max[idx] / _mdl_u_nominal[idx];
	else if (_mdl_u.max[idx] < Inf)
	  x.max[i] = _mdl_u.max[idx] / _mdl_u_nominal[idx];
        // integer variables
        x.integer[i] = _mdl_u.integer[idx];
      }
      i++;
    }
  }
  for (idx = 0, i = 0; i < _nx; idx++, i++) {
    if (i == _mdl_nd && _nu > 0) {
      // skip control states
      i += _nu - 1;
      idx--;
      continue;
    }
    x.initial[i] = _mdl_xs[ks(k)][idx] / _mdl_x_nominal[idx];
    if (k == 0) {
      if (_mdl_x_periodic[idx])
        // indicate periodic state to Hqp_Docp
        x.min[i] = x.max[i] = Inf;
      else if (_mdl_x0_active[idx]) {
	if (!_multistage)
	  m_error(E_FORMAT, "Prg_DynamicOpt::setup: "
		  "mdl_x0_active=1 requires prg_multistage=true");
        // use the more restrictive bound of x0_min/max and x_min/max
        if (_mdl_x0.min[idx] > _mdl_x.min[idx])
          x.min[i] = _mdl_x0.min[idx] / _mdl_x_nominal[idx];
        else
          x.min[i] = _mdl_x.min[idx] / _mdl_x_nominal[idx];
        if (_mdl_x0.max[idx] < _mdl_x.max[idx])
          x.max[i] = _mdl_x0.max[idx] / _mdl_x_nominal[idx];
        else
          x.max[i] = _mdl_x.max[idx] / _mdl_x_nominal[idx];
        // integer variables
        x.integer[i] = _mdl_x.integer[idx];
      }
      else
	x.min[i] = x.max[i] = x.initial[i];
    }
    else {
      x.min[i] = _mdl_x.min[idx] / _mdl_x_nominal[idx];
      x.max[i] = _mdl_x.max[idx] / _mdl_x_nominal[idx];
      // integer variables
      x.integer[i] = _mdl_x.integer[idx];
    }
  }

  // setup control inputs
  for (i = 0, isc = 0, idx = 0; idx < _mdl_nu; idx++) {
    if (_mdl_u.active[idx]) {
      if (!_multistage && k < _K) {
	// treat control bounds via general constraints
	if (_mdl_u.min[idx] > -Inf) {
	  for (j = max(0, _mdl_u0_nfixed[idx] - _mdl_u_order[idx]); j < upsk; j++)
	    c.min[spsk*(_nc+_nsc) + upsk*_nsuc + i + j] =
	      _mdl_u.min[idx] / _mdl_u_nominal[idx];
	}
	if (_mdl_u.max[idx] < Inf) {
	  for (j = max(0, _mdl_u0_nfixed[idx] - _mdl_u_order[idx]); j < upsk; j++)
	    c.max[spsk*(_nc+_nsc) + upsk*_nsuc + i + j] =
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
        // rate of change soft bounds
        if (_mdl_der_u_soft.active[idx]) {
          if (_mdl_der_u_soft.min[idx] > -Inf) {
            for (j = 0; j < upsk; j++)
              c.min[spsk*(_nc+_nsc) + isc + j] =
                _mdl_der_u_soft.min[idx] / _mdl_u_nominal[idx]*_t_nominal;
            isc += upsk;
          }
          if (_mdl_der_u_soft.max[idx] < Inf) {
            for (j = 0; j < upsk; j++)
              c.max[spsk*(_nc+_nsc) + isc + j] =
                _mdl_der_u_soft.max[idx] / _mdl_u_nominal[idx]*_t_nominal;
            isc += upsk;
          }
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
	  if (_mdl_u_decimation[idx] == 0 ||
              (k + 1) % _mdl_u_decimation[idx] != 0) {
	    u.min[i] = u.max[i] = 0.0;
	  }
	}
	else {
	  for (j = 0; j < upsk; j++)
	    if (_mdl_u_decimation[idx] == 0 ||
                (j + 1) % _mdl_u_decimation[idx] != 0)
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
    if (_mdl_y_soft.active[idx]) {
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
  }

  // setup slack variables for soft constraints
  if (k < _K)
    for (i = upsk*_nu; i < upsk*(_nu+_nsu) + spsk*_ns; i++) {
      u.min[i] = 0.0;
      u.initial[i] = 0.0;
    }
  else
    for (i = _nx; i < _nx + spsk*_ns + _nsf; i++) {
      x.min[i] = 0.0;
      x.initial[i] = 0.0;
    }

  // setup state and output constraints at initial time
  if (k == 0) {
    for (i = spsk*(_nc+_nsc) + upsk*_nsuc, idx = 0; idx < _mdl_nx; idx++) {
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
    for (i = spsk*(_nc+_nsc), iscf = 0, idx = 0; idx < _mdl_ny; idx++) {
      if (_mdl_yf.active[idx]) {
        if (_mdl_yf.min[idx] > -Inf)
          c.min[i] = _mdl_yf.min[idx] / _mdl_y_nominal[idx];
        if (_mdl_yf.max[idx] < Inf)
          c.max[i] = _mdl_yf.max[idx] / _mdl_y_nominal[idx];
        i++;
      }
      // soft constraints at final time
      if (_mdl_yf_soft.active[idx]) {
        if (_mdl_yf_soft.min[idx] > -Inf) {
          c.min[spsk*(_nc+_nsc) + _ncf + iscf] =
            _mdl_yf_soft.min[idx] / _mdl_y_nominal[idx];
          iscf++;
        }
        if (_mdl_yf_soft.max[idx] < Inf) {
          c.max[spsk*(_nc+_nsc) + _ncf + iscf] =
            _mdl_yf_soft.max[idx] / _mdl_y_nominal[idx];
          iscf++;
        }
      }
    }
  }
}

//--------------------------------------------------------------------------
void Prg_DynamicOpt::setup_struct(int k,
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

  // setup struct of Jacobians for wrt. states and active inputs
  if (k < _K) {
    m_zero(F.Jdx);
    m_zero(F.Ju);

    for (i = 0; i < _mdl_nd; i++) {
      // no dependency of discrete states on states 
      for (j = 0; j < _nx; j++)
        F.Jx[i][j] = 0.0;
    }

    for (i = _mdl_nd, idx = 0; idx < _mdl_nu; idx++) {
      if (_mdl_u.active[idx]) {
	// no dependency of active inputs on states
	for (j = 0; j < _nx; j++)
	  F.Jx[i][j] = 0.0;
	if (_mdl_u_order[idx] == 0) {
	  // generally no dependecy on time derivative for zero order hold
          // it is still integrated to get contiguous vector though
	  F.Jdx[i][i] = -1.0;
	  // F.Ju is zero for zero order hold
	}
        else {
          F.Jdx[i][i] = -1.0; // continuous integration for first order hold
	  if (_multistage) {
	    // F.Ju is constant for multistage and not _t_active
	    F.Ju[i][i-_mdl_nd] = 1.0/_t_nominal;
            // an additional element needs to be allocated for _t_active
            if (_t_active)
              F.Ju[i][_t_scale_i] = 1.0;
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
    if (_multistage && !_t_active)
      F.set_linear(Omu_Dependent::WRT_u);

    // explicit ODE for continuous-time equations
    for (i = _mdl_nd + _nu; i < _nx; i++)
      F.Jdx[i][i] = -1.0;
    F.set_linear(Omu_Dependent::WRT_dx);
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

  // slack variables appear linear in all optimization constraints
  if (_multistage) {
    if (k < _K) {
      for (j = _nu; j < u->dim; j++) {
        F.set_linear_variable(Omu_Dependent::WRT_u, j);
        f.set_linear_variable(Omu_Dependent::WRT_u, j);
        c.set_linear_variable(Omu_Dependent::WRT_u, j);
        for (i = 0; i < F->dim; i++)
          F.Ju[i][j] = 0.0;
        for (i = 0; i < f->dim; i++)
          f.Ju[i][j] = 0.0;
      }
    }
    else {
      for (j = _nx; j < x->dim; j++) {
        c.set_linear_variable(Omu_Dependent::WRT_x, j);
      }
    }
  }
}

//--------------------------------------------------------------------------
void Prg_DynamicOpt::init_simulation(int k,
				     Omu_VariableVec &x, Omu_VariableVec &u)
{
  // initial states
  if (k == 0 || _multistage > 1)
    v_copy(x.initial, x);
}

//--------------------------------------------------------------------------
void Prg_DynamicOpt::update(int kk,
			    const Omu_StateVec &x, const Omu_Vec &u,
			    const Omu_StateVec &xf,
			    Omu_DependentVec &f, Omu_Dependent &f0,
			    Omu_DependentVec &c)
{
  int i, idx, is, isc, isf, iscf;
  int spsk = kk < _KK? (_multistage? _sps: _KK): 1; // one sample at final time
  int upsk = _multistage? 1: _KK;
  double tscale_1 = _t_active? x[_mdl_nd + _t_scale_i]* _t_scale_nominal: 1.0;
  double tscale = _t_active && kk < _KK?
    (x[_mdl_nd + _t_scale_i] + u[_t_scale_i*upsk + kk%upsk])* _t_scale_nominal: tscale_1;
  double help;

  // initialize model inputs
  real_T *mdl_u = NULL;
  if (ssGetNumInputPorts(_SS) > 0) {
    if (ssGetInputPortRequiredContiguous(_SS, 0))
      mdl_u = (real_T *)ssGetInputPortRealSignal(_SS, 0);
    else
      mdl_u = (real_T *)*ssGetInputPortRealSignalPtrs(_SS, 0);
  }
  for (i = _mdl_nd, idx = 0; idx < _mdl_nu; idx++) {
    if (_mdl_u.active[idx]) {
      mdl_u[idx] = x[i++] * _mdl_u_nominal[idx];
      _mdl_us[kk][idx] = mdl_u[idx];
    }
    else if (kk == _KK && _KK > 0 && _mdl_u_order[idx] == 0)
      // hold last but one value in case of zero order hold
      mdl_u[idx] = _mdl_us[kk-1][idx];
    else
      mdl_u[idx] = _mdl_us[kk][idx];
  }

  // set simulation time
  ssSetT(_SS, _taus[kk]);

  // pass current states to model
  real_T *mdl_xd = ssGetDiscStates(_SS);
  real_T *mdl_xc = ssGetContStates(_SS);
  for (i = 0; i < _mdl_nd; i++)
    mdl_xd[i] = x[i] * _mdl_x_nominal[i];
  for (i = _mdl_nd; i < _mdl_nx; i++)
    mdl_xc[i - _mdl_nd] = x[_nu + i] * _mdl_x_nominal[i];

  // call mdlOutputs
  SMETHOD_CALL2(mdlOutputs, _SS, 0); 

  // obtain model outputs
  real_T *mdl_y = NULL;
  if (ssGetNumOutputPorts(_SS) > 0)
    mdl_y = ssGetOutputPortRealSignal(_SS, 0);

  // correct model outputs with bias
  for (idx = 0; idx < _mdl_ny; idx++)
    mdl_y[idx] += _mdl_y_bias[idx];

  // utilize model outputs for constraints and optimization objective
  f0 = 0.0;

  double dt; 	// factor for linear interpol. (trapezoidal integration rule)
  double dt0; 	// factor for zero order hold
  double dtu; 	// factor used for controlled inputs in objective function
  double dty; 	// factor used for outputs in objective function
  if (kk == 0) {
    if (_KK > 0) {
      dt = 0.5 * (ts(kk+1) - ts(kk)) * tscale;
      dt0 = (ts(kk+1) - ts(kk)) * tscale;
    } else {
      // steady-state problem
      dt = dt0 = 1.0;
    }
  } else if (kk == _KK) {
    dt = 0.5 * (ts(kk) - ts(kk-1)) * tscale_1;
    dt0 = 0.0;
  } else {
    dt = 0.5 * ((ts(kk) - ts(kk-1)) * tscale_1 + (ts(kk+1) - ts(kk)) * tscale);
    dt0 = (ts(kk+1) - ts(kk)) * tscale;
  }

  // contribution of controlled inputs
  for (i = 0, is = 0, isc = 0, idx = 0; idx < _mdl_nu; idx++) {
    if (_mdl_u.active[idx]) {
      // control inputs
      dtu = (_mdl_u_order[idx] == 0)? dt0: dt;
      help = (_mdl_us[kk][idx] - _mdl_u.ref[idx]) / _mdl_u_nominal[idx];
      f0 += dtu * _mdl_u.weight1[idx] * help;
      f0 += dtu * _mdl_u.weight2[idx] * help*help;
      // rates of change
      if (kk < _KK) {
        _mdl_der_u[idx] = u[i*upsk + kk%upsk]*_mdl_u_nominal[idx]/_t_nominal;
	help = u[i*upsk + kk%upsk]/_t_nominal - _mdl_der_u.ref[idx]
	  / _mdl_u_nominal[idx];
	f0 += dt0 * _mdl_der_u.weight1[idx] * help;
	f0 += dt0 * _mdl_der_u.weight2[idx] * help*help;
        // soft constraints on rates of change
        if (_mdl_der_u_soft.active[idx]) {
          help = u[upsk*_nu + is + kk%upsk]/_t_nominal;
          f0 += dt0 * (_mdl_der_u_soft.weight1[idx]*help
                       + _mdl_der_u_soft.weight2[idx]*help*help);
          if (_mdl_der_u_soft.min[idx] > -Inf) {
            c[spsk*(_nc+_nsc) + isc + kk%upsk] = u[i*upsk + kk%upsk]/_t_nominal + help;
            isc += upsk;
          }
          if (_mdl_der_u_soft.max[idx] < Inf) {
            c[spsk*(_nc+_nsc) + isc + kk%upsk] = u[i*upsk + kk%upsk]/_t_nominal - help;
            isc += upsk;
          }
          is += upsk;
        }
      }
      i++;
    }
  }
  // if not multistage, treat control bounds via general constraints
  // (note: do not access xf to avoid non-linear dependency in update)
  if (!_multistage && kk < _KK) {
    for (i = 0; i < _nu; i++)
      c[spsk*(_nc+_nsc) + upsk*(_nsuc + i) + kk%upsk] = x[i] +
	(ts(kk+1)-ts(kk)) * u[i*upsk + kk%upsk]/_t_nominal; // == xf[i]
  }
  // contribution of active outputs
  for (i = 0, is = 0, isc = 0, idx = 0; idx < _mdl_ny; idx++) {
    dty = (_mdl_y_order[idx] == 0)? dt0: dt;
    // assign used outputs to constraints
    if (_mdl_y.active[idx]) {
      c[i + kk%spsk] = mdl_y[idx] / _mdl_y_nominal[idx];
      i += spsk;
    }
    // calculate objective term
    help = (mdl_y[idx] - _mdl_y.ref[idx]) / _mdl_y_nominal[idx];
    if (_mdl_y.weight1[idx] != 0.0)
      f0 += dty * _mdl_y.weight1[idx] * help;
    if (_mdl_y.weight2[idx] != 0.0)
      f0 += dty * _mdl_y.weight2[idx] * help*help;
    // consider soft constraints
    if (_mdl_y_soft.active[idx] == 1) {
      if (kk < _KK)
	help = u[upsk*(_nu+_nsu) + is + kk%spsk];
      else
	help = x[_nx + is + kk%spsk];

      f0 += dty * (_mdl_y_soft.weight1[idx]*help
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
    i = spsk*(_nc+_nsc) + upsk*_nsuc;
    // constraints on derivatives of initial states
    if (_mdl_nx > _mdl_nd) {
      SMETHOD_CALL(mdlDerivatives, _SS);
      real_T *mdl_dx = ssGetdX(_SS);
      for (idx = 0; idx < _mdl_nx; idx++) {
        if (_mdl_x0_active[idx]
            && (_mdl_der_x0_min[idx] > -Inf || _mdl_der_x0_max[idx] < Inf)) {
          c[i] = tscale * mdl_dx[idx] / _mdl_x_nominal[idx];
          i++;
        }
      }
    }
    for (idx = 0; idx < _mdl_ny; idx++) {
      // assign used outputs to constraints
      if (_mdl_y0.active[idx]) {
	c[i] = mdl_y[idx] / _mdl_y_nominal[idx];
	i++;
      }
      // initial objective terms
      help = (mdl_y[idx] - _mdl_y.ref[idx]) / _mdl_y_nominal[idx];
      if (_mdl_y0.weight1[idx] != 0.0)
	f0 += _mdl_y0.weight1[idx] * help;
      if (_mdl_y0.weight2[idx] != 0.0)
	f0 += _mdl_y0.weight2[idx] * help*help;
    }
  }

  // additional terms at final time
  if (kk == _KK) {
    for (i = _nc+_nsc, isf = 0, iscf = 0, idx = 0; idx < _mdl_ny; idx++) {
      // assign used outputs to constraints
      if (_mdl_yf.active[idx]) {
	c[i] = mdl_y[idx] / _mdl_y_nominal[idx];
	i++;
      }
      // final objective terms
      help = (mdl_y[idx] - _mdl_y.ref[idx]) / _mdl_y_nominal[idx];
      if (_mdl_yf.weight1[idx] != 0.0)
	f0 += _mdl_yf.weight1[idx] * help;
      if (_mdl_yf.weight2[idx] != 0.0)
	f0 += _mdl_yf.weight2[idx] * help*help;
      // soft constraints at final time
      if (_mdl_yf_soft.active[idx] == 1) {
        help = x[_nx + spsk*_ns + isf];

        f0 += _mdl_yf_soft.weight1[idx]*help
                + _mdl_yf_soft.weight2[idx]*help*help;

        if (_mdl_yf_soft.min[idx] > -Inf) {
          c[spsk*(_nc+_nsc) + _ncf + iscf] = mdl_y[idx] / _mdl_y_nominal[idx] + help;
          iscf++;
        }
        if (_mdl_yf_soft.max[idx] < Inf) {
          c[spsk*(_nc+_nsc) + _ncf + iscf] = mdl_y[idx] / _mdl_y_nominal[idx] - help;
          iscf++;
        }
        isf++;
      }
    }
  }

  // store values of model states
  for (i = 0; i < _mdl_nd; i++)
    _mdl_xs[kk][i] = value(mdl_xd[i]);
  for (; i < _mdl_nx; i++)
    _mdl_xs[kk][i] = value(mdl_xc[i - _mdl_nd]);

  // store values of model outputs
  for (i = 0; i < _mdl_ny; i++)
    _mdl_ys[kk][i] = value(mdl_y[i]);

  // store initial and final values
  if (kk == 0) {
    for (i = 0; i < _mdl_nd; i++)
      _mdl_x0[i] = value(mdl_xd[i]);
    for (; i < _mdl_nx; i++)
      _mdl_x0[i] = value(mdl_xc[i - _mdl_nd]);
    for (i = 0; i < _mdl_nu; i++)
      _mdl_u0[i] = value(mdl_u[i]);
    for (i = 0; i < _mdl_ny; i++)
      _mdl_y0[i] = value(mdl_y[i]);
  }
  else if (kk == _KK) {
    for (i = 0; i < _mdl_nu; i++)
      _mdl_uf[i] = value(mdl_u[i]);
    for (i = 0; i < _mdl_ny; i++)
      _mdl_yf[i] = value(mdl_y[i]);
  }

  // junction conditions for subsequent stage
  if (kk < _KK) {
    if (_mdl_nd > 0) {
      setContinuousTask(false);
      setSampleHit(true);
      // call mdlUpdate to get discrete events processed
      if (ssGetmdlUpdate(_SS) != NULL) {
        SMETHOD_CALL2(mdlUpdate, _SS, 0);
      }
      // read discrete states from model
      for (i = 0; i < _mdl_nd; i++) {
        f[i] = mdl_xd[i] / _mdl_x_nominal[i];
      }
    }
    // update controlled inputs from optimizer
    for (i = _mdl_nd, idx = 0; idx < _mdl_nu; idx++) {
      if (_mdl_u.active[idx]) {
	if (_mdl_u_order[idx] == 0 && (!_multistage || (kk+1)%spsk == 0)) {
	  // zero order hold at end of stage:
	  // apply step in u for subsequent stage
          help = ((ts(kk+1) - ts(kk-spsk/upsk+1))/_t_nominal
                  * u[(i-_mdl_nd)*upsk + kk%upsk]);
          if (_t_active && i-_mdl_nd != _t_scale_i)
            f[i] = x[i] + tscale * help;
          else
            f[i] = x[i] + help;
        }
	else
	  // piecewise linear interpolation or within stage with zoh
	  f[i] = xf[i];
	i++;
      }
    }
    assert(i == _mdl_nd + _nu); // problem structure must not have changed
    // read continuous states from model
    for (; i < _nx; i++)
      f[i] = xf[i]; //mdl_xc[i - _mdl_nd - _nu];
    // re-use slack variables of last but one sample period at final time
    // note: this is only because no control parameters allowed at final time
    // and because states need to be defined with state equations
    if (kk == _KK-1) {
      for (; i < _nx + _ns; i++)
	f[i] = u[upsk*(_nu+_nsu) + (i-_nx+1)*spsk-1];
      for (; i < _nx + _ns + _nsf; i++)
	f[i] = u[upsk*(_nu+_nsu) + spsk*_ns + i - _nx - _ns];
    }
  }

  // obtain Jacobians if required
  if (f.is_required_J() || f0.is_required_g() || c.is_required_J()) {
    update_grds(kk, x, u, xf, f, f0, c);
  }
}

//--------------------------------------------------------------------------
void Prg_DynamicOpt::update_grds(int kk, 
				 const Omu_StateVec &x, const Omu_Vec &u,
				 const Omu_StateVec &xf,
				 Omu_DependentVec &f, Omu_Dependent &f0,
				 Omu_DependentVec  &c)
{
  int i, idx, ii, iidx0, iidx, isc, iscdx, is, j, jdx, rdx;
  int isf, iscf, iscfdx, ixf;
  int spsk = kk < _KK? (_multistage? _sps: _KK): 1; // one sample at final time
  int upsk = _multistage? 1: _KK;
  double tscale_1 = _t_active? x[_mdl_nd + _t_scale_i]* _t_scale_nominal: 1.0;
  double tscale = _t_active && kk < _KK?
    (x[_mdl_nd + _t_scale_i] + u[_t_scale_i*upsk + kk%upsk])* _t_scale_nominal: tscale_1;
  int t_scale_iu = _t_scale_i*upsk + kk%upsk;
  int t_scale_ix = _mdl_nd + _t_scale_i;
  double help;
  int analyticJac = 0;

  if (ssGetmdlJacobian(_SS) == NULL || _t_active) {
    // todo: use analytic Jacobian if _t_active
    // call predefined update_grds for numerical differentiation
    Omu_Program::update_grds(kk, x, u, xf, f, f0, c);
  }
  else {
    analyticJac = 1;
    // exploit mdlJacobian to obtain c.Jx
    int mdl_u_idx, mdl_x_idx, mdl_dx_idx, mdl_y_idx;
    int mdl_nc = _mdl_nx - _mdl_nd; // number of continuous states
    real_T *pr = ssGetJacobianPr(_SS);
    int_T *ir = ssGetJacobianIr(_SS);
    int_T *jc = ssGetJacobianJc(_SS);

    // restore states after mdlUpdate has been processed
    // call mdlOutputs
    real_T *mdl_xd = ssGetDiscStates(_SS);
    real_T *mdl_xc = ssGetContStates(_SS);
    for (i = 0; i < _mdl_nd; i++)
      mdl_xd[i] = x[i] * _mdl_x_nominal[i];
    for (i = _mdl_nd; i < _mdl_nx; i++)
      mdl_xc[i - _mdl_nd] = x[_nu + i] * _mdl_x_nominal[i];
    SMETHOD_CALL2(mdlOutputs, _SS, 0); 

    mdlJacobian(_SS);

    m_zero(c.Jx);
    m_zero(c.Ju);
    m_zero(f.Jx);
    // Jacobian wrt S-function states (ddxdy/dxcxd)
    for (jdx = 0; jdx < _mdl_nx; jdx++) {
      j = jdx < mdl_nc? _mdl_nd + _nu + jdx: jdx - mdl_nc;
      mdl_x_idx = jdx < mdl_nc? _mdl_nd + jdx: jdx - mdl_nc;
      for (i = 0, idx = _mdl_nx, isc = spsk*_nc, iscdx = _mdl_nx,
             iscf = spsk*(_nc+_nsc) + _ncf, iscfdx = _mdl_nx,
	     ii = _nc+_nsc, iidx0 = 0, iidx = _mdl_nx,
	     rdx = jc[jdx]; rdx < jc[jdx+1]; rdx++) {
        if (ir[rdx] < mdl_nc) {
          mdl_dx_idx = _mdl_nd + ir[rdx];
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
        else if (ir[rdx] < _mdl_nx && kk < _KK) {
          // junction conditions for discrete states
          ixf = ir[rdx] - mdl_nc;
          f.Jx[ixf][j] = pr[rdx] /
            _mdl_x_nominal[ixf] * _mdl_x_nominal[mdl_x_idx];
        }
        else if (ir[rdx] >= _mdl_nx) {
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
	  // soft constraints at final time
	  if (kk == _KK && _mdl_yf_soft.active[mdl_y_idx] == 1) {
	    // need to loop to obtain iscf considering active outputs
	    for (; iscfdx < ir[rdx]; iscfdx++) {
	      if (_mdl_yf_soft.min[iscfdx - _mdl_nx] > -Inf)
		iscf++;
	      if (_mdl_yf_soft.max[iscfdx - _mdl_nx] < Inf)
		iscf++;
	    }
	    if (_mdl_yf_soft.min[mdl_y_idx] > -Inf) {
	      c.Jx[iscf][j] = pr[rdx] /
		_mdl_y_nominal[mdl_y_idx] * _mdl_x_nominal[mdl_x_idx];
	      iscf++;
	    }
	    if (_mdl_yf_soft.max[mdl_y_idx] < Inf) {
	      c.Jx[iscf][j] = pr[rdx] /
		_mdl_y_nominal[mdl_y_idx] * _mdl_x_nominal[mdl_x_idx];
	    }
	    if (_mdl_yf_soft.min[mdl_y_idx] > -Inf) {
	      iscf--;
	    }
	  }
	}
      }
    }
    // Jacobian wrt S-function inputs (ddxdy/du)
    // (note: these are states for the optimizer)
    for (j = _mdl_nd, jdx = _mdl_nx; jdx < _mdl_nx + _mdl_nu; jdx++) {
      mdl_u_idx = jdx - _mdl_nx;
      if (_mdl_u.active[mdl_u_idx]) {
	for (i = 0, idx = _mdl_nx, isc = spsk*_nc, iscdx = _mdl_nx,
               iscf = spsk*(_nc+_nsc) + _ncf, iscfdx = _mdl_nx,
	       ii = _nc+_nsc, iidx0 = 0, iidx = _mdl_nx,
	       rdx = jc[jdx]; rdx < jc[jdx+1]; rdx++) {
          if (ir[rdx] < mdl_nc) {
            mdl_dx_idx = _mdl_nd + ir[rdx];
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
	  else if (ir[rdx] < _mdl_nx && kk < _KK) {
            // junction conditions for discrete states
            ixf = ir[rdx] - mdl_nc;
            f.Jx[ixf][j] = pr[rdx] /
              _mdl_x_nominal[ixf] * _mdl_u_nominal[mdl_u_idx];
          }
	  else if (ir[rdx] >= _mdl_nx) {
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
	    // soft constraints at final time
	    if (kk == _KK && _mdl_yf_soft.active[mdl_y_idx] == 1) {
	      // need to loop to obtain iscf considering active outputs
	      for (; iscfdx < ir[rdx]; iscfdx++) {
		if (_mdl_yf_soft.min[iscfdx - _mdl_nx] > -Inf)
		  iscf++;
		if (_mdl_yf_soft.max[iscfdx - _mdl_nx] < Inf)
		  iscf++;
	      }
	      if (_mdl_yf_soft.min[mdl_y_idx] > -Inf) {
		c.Jx[iscf][j] = pr[rdx] /
		  _mdl_y_nominal[mdl_y_idx] * _mdl_u_nominal[mdl_u_idx];
		iscf++;
	      }
	      if (_mdl_yf_soft.max[mdl_y_idx] < Inf) {
		c.Jx[iscf][j] = pr[rdx] /
		  _mdl_y_nominal[mdl_y_idx] * _mdl_u_nominal[mdl_u_idx];
	      }
	      if (_mdl_yf_soft.min[mdl_y_idx] > -Inf) {
		iscf--;
	      }
	    }
	  }
	}
	j++;
      }
    }
    // contribution of slack variables for soft constraints
    for (is = 0, isc = spsk*_nc, isf = 0, iscf = 0,
           idx = 0; idx < _mdl_ny; idx++) {
      if (_mdl_y_soft.active[idx] == 1) {
	if (_mdl_y_soft.min[idx] > -Inf) {
	  if (kk < _KK)
	    c.Ju[isc + kk%spsk][upsk*(_nu+_nsu) + is + kk%spsk] += 1.0;
	  else
	    c.Jx[isc + kk%spsk][_nx + is + kk%spsk] += 1.0;
	  isc += spsk;
	}
	if (_mdl_y_soft.max[idx] < Inf) {
	  if (kk < _KK)
	    c.Ju[isc + kk%spsk][upsk*(_nu+_nsu) + is + kk%spsk] -= 1.0;
	  else
	    c.Jx[isc + kk%spsk][_nx + is + kk%spsk] -= 1.0;
	  isc += spsk;
	}
	is += spsk;
      }
      // soft constraints at final time
      if (kk == _KK && _mdl_yf_soft.active[idx] == 1) {
	if (_mdl_yf_soft.min[idx] > -Inf) {
          c.Jx[spsk*(_nc+_nsc) + _ncf + iscf][_nx + spsk*_ns + isf] += 1.0;
	  iscf++;
	}
	if (_mdl_yf_soft.max[idx] < Inf) {
          c.Jx[spsk*(_nc+_nsc) + _ncf + iscf][_nx + spsk*_ns + isf] -= 1.0;
	  iscf++;
	}
	isf++;
      }
    }
    // control bounds if not multistage
    if (!_multistage && kk < _KK) {
      for (i = 0; i < _nu; i++) {
	c.Jx[spsk*(_nc+_nsc) + upsk*(_nsuc + i) + kk%upsk][i] = 1.0;
	c.Ju[spsk*(_nc+_nsc) + upsk*(_nsuc + i) + kk%upsk][i*upsk + kk%upsk] =
	  (ts(kk+1)-ts(kk))/_t_nominal;
      }
    }
  }

  // junction conditions for continuous-time equations
  // Note: f.Jx for discrete states is filled above
  if (kk < _KK) {
    // default values
    m_zero(f.Ju);
    m_zero(f.Jxf);
    // modifications for controlled inputs with zero order hold
    for (i = _mdl_nd, idx = 0; idx < _mdl_nu; idx++) {
      if (_mdl_u.active[idx]) {
	if (_mdl_u_order[idx] == 0 && (!_multistage || (kk+1)%spsk == 0)) {
	  // zero order hold at end of stage:
	  // apply step in u for subsequent stage
	  f.Jx[i][i] = 1.0;
          help = (ts(kk+1) - ts(kk-spsk/upsk+1))/_t_nominal;
          if (_t_active && i-_mdl_nd != _t_scale_i) {
            f.Ju[i][(i-_mdl_nd)*upsk + kk%upsk] = tscale*help;
            f.Jx[i][t_scale_ix] = _t_scale_nominal*help*u[(i-_mdl_nd)*upsk + kk%upsk];
            f.Ju[i][t_scale_iu] = f.Jx[i][t_scale_ix];
          }
          else
            f.Ju[i][(i-_mdl_nd)*upsk + kk%upsk] = help;
	}
        else {
          f.Jxf[i][i] = 1.0;
        }
	i++;
      }
    }
    assert(i == _mdl_nd + _nu); // problem structure must not have changed
    // continuous states
    for (; i < _nx; i++)
      f.Jxf[i][i] = 1.0;

    if (kk == _KK-1) {
      // slack variables at final time
      for (i = _nx; i < _nx + _ns; i++)
	// note: f.Jxf has only _nx columns, but _nx+_ns rows
	// that is why don't set f.Jxf[i][i] = 0.0;
	f.Ju[i][upsk*(_nu+_nsu) + (i-_nx+1)*spsk-1] = 1.0;
      for (; i < _nx + _ns; i++)
	f.Ju[i][upsk*(_nu+_nsu) + spsk*_ns + i - _nx - _ns] = 1.0;
    }
  }

  // apply chain rule to calculate gradient of f0
  // (and Jacobian for soft constraints on rates of change)
  // (do this as well if Omu_Program::update_grds was used, as complete
  //  numeric differentiation gives bad results for quadratic terms)
  double dt; 	// factor for linear interpol. (trapezoidal integration rule)
  double dt0; 	// factor for zero order hold
  double dtu; 	// factor used for controlled inputs in objective function
  double dty; 	// factor used for outputs in objective function
  double ddtx; 	// help variable for partial derivative of a dt
  double ddtu; 	// help variable for partial derivative of a dt
  double ddt0; 	// help variable for partial derivative of a dt
  if (kk == 0) {
    if (_KK > 0) {
      dt = 0.5 * (ts(kk+1) - ts(kk)) * tscale;
      ddtx = 0.5 * (ts(kk+1) - ts(kk)) * _t_scale_nominal;
      ddtu = 0.5 * (ts(kk+1) - ts(kk)) * _t_scale_nominal;
      dt0 = (ts(kk+1) - ts(kk)) * tscale;
      ddt0 = (ts(kk+1) - ts(kk)) * _t_scale_nominal;
    } else {
      // steady-state problem
      dt = dt0 = 1.0;
      ddtx = 0.0;
      ddtu = 0.0;
      ddt0 = 0.0;
    }
  } else if (kk == _KK) {
    dt = 0.5 * (ts(kk) - ts(kk-1)) * tscale_1;
    ddtx = 0.5 * (ts(kk) - ts(kk-1)) * _t_scale_nominal;
    ddtu = 0.0;
    dt0 = 0.0;
    ddt0 = 0.0;
  } else {
    dt = 0.5 * ((ts(kk) - ts(kk-1)) * tscale_1 + (ts(kk+1) - ts(kk)) * tscale);
    ddtx = 0.5 * (ts(kk+1) - ts(kk-1)) * _t_scale_nominal;
    ddtu = 0.5 * (ts(kk+1) - ts(kk)) * _t_scale_nominal;
    dt0 = (ts(kk+1) - ts(kk)) * tscale;
    ddt0 = (ts(kk+1) - ts(kk)) * _t_scale_nominal;
  }

  v_zero(f0.gx);
  v_zero(f0.gu);
  // controlled inputs
  for (i = _mdl_nd, is = 0, isc = 0, idx = 0; idx < _mdl_nu; idx++) {
    // contribution of objective term
    if (_mdl_u.active[idx]) {
      // controlled inputs
      dtu = (_mdl_u_order[idx] == 0)? dt0: dt;
      help = (x[i] - _mdl_u.ref[idx]/_mdl_u_nominal[idx]);
      f0.gx[i] += dtu * _mdl_u.weight1[idx];
      f0.gx[i] += dtu * _mdl_u.weight2[idx] * 2.0 * help;
      if (_t_active) {
        if (_mdl_u_order[idx] == 0) {
          f0.gx[t_scale_ix] += ddt0 * _mdl_u.weight1[idx] * help;
          f0.gx[t_scale_ix] += ddt0 * _mdl_u.weight2[idx] * help*help;
          if (kk < _KK) {
            f0.gu[t_scale_iu] += ddt0 * _mdl_u.weight1[idx] * help;
            f0.gu[t_scale_iu] += ddt0 * _mdl_u.weight2[idx] * help*help;
          }
        }
        else if (_mdl_u_order[idx] > 0) {
          f0.gx[t_scale_ix] += ddtx * _mdl_u.weight1[idx] * help;
          f0.gx[t_scale_ix] += ddtx * _mdl_u.weight2[idx] * help*help;
          if (kk < _KK) {
            f0.gu[t_scale_iu] += ddtu * _mdl_u.weight1[idx] * help;
            f0.gu[t_scale_iu] += ddtu * _mdl_u.weight2[idx] * help*help;
          }
        }
      }
      // rates of change
      if (kk < _KK) {
	help = u[(i-_mdl_nd)*upsk + kk%upsk]/_t_nominal - _mdl_der_u.ref[idx]
	  / _mdl_u_nominal[idx];
	f0.gu[(i-_mdl_nd)*upsk + kk%upsk] += dt0 * _mdl_der_u.weight1[idx]/_t_nominal;
	f0.gu[(i-_mdl_nd)*upsk + kk%upsk] += dt0 * _mdl_der_u.weight2[idx]/_t_nominal * 2.0 * help;
        if (_t_active) {
          f0.gx[t_scale_ix] += ddt0 * _mdl_der_u.weight1[idx] * help;
          f0.gx[t_scale_ix] += ddt0 * _mdl_der_u.weight2[idx] * help*help;
          f0.gu[t_scale_iu] += ddt0 * _mdl_der_u.weight1[idx] * help;
          f0.gu[t_scale_iu] += ddt0 * _mdl_der_u.weight2[idx] * help*help;
        }
        // soft constraints on rates of change
        if (_mdl_der_u_soft.active[idx]) {
          help = u[upsk*_nu + is + kk%upsk]/_t_nominal;
          f0.gu[upsk*_nu + is + kk%upsk] += dt0 * _mdl_der_u_soft.weight1[idx]/_t_nominal;
          f0.gu[upsk*_nu + is + kk%upsk] += dt0 * _mdl_der_u_soft.weight2[idx]/_t_nominal * 2.0 * help;
          if (_t_active) {
            f0.gx[t_scale_ix] += ddt0 * _mdl_der_u_soft.weight1[idx] * help;
            f0.gx[t_scale_ix] += ddt0 * _mdl_der_u_soft.weight2[idx] * help*help;
            f0.gu[t_scale_iu] += ddt0 * _mdl_der_u_soft.weight1[idx] * help;
            f0.gu[t_scale_iu] += ddt0 * _mdl_der_u_soft.weight2[idx] * help*help;
          }
          if (analyticJac && _mdl_der_u_soft.min[idx] > -Inf) {
            c.Ju[spsk*(_nc+_nsc) + isc + kk%upsk][i*upsk + kk%upsk]
              += 1.0/_t_nominal;
            c.Ju[spsk*(_nc+_nsc) + isc + kk%upsk][upsk*_nu + is + kk%upsk]
              += 1.0/_t_nominal;
            isc += upsk;
          }
          if (analyticJac && _mdl_der_u_soft.max[idx] < Inf) {
            c.Ju[spsk*(_nc+_nsc) + isc + kk%upsk][i*upsk + kk%upsk]
              += 1.0/_t_nominal;
            c.Ju[spsk*(_nc+_nsc) + isc + kk%upsk][upsk*_nu + is + kk%upsk]
              -= 1.0/_t_nominal;
            isc += upsk;
          }
          is += upsk;
        }
      }
      i++;
    }
  }
  // active outputs
  for (i = 0, is = 0, idx = 0; idx < _mdl_ny; idx++) {
    dty = (_mdl_y_order[idx] == 0)? dt0: dt;
    // contribution of objective term
    if (_mdl_y.active[idx] == 1) {
      help = (c[i + kk%spsk] - _mdl_y.ref[idx]/_mdl_y_nominal[idx]);
      if (_mdl_y.weight1[idx] != 0.0) {
        for (j = 0; j < _nx; j++)
          f0.gx[j] += dty * _mdl_y.weight1[idx] * c.Jx[i + kk%spsk][j];
        if (_t_active) {
          f0.gx[t_scale_ix] += ddtx * _mdl_y.weight1[idx] * help;
          if (kk < _KK)
            f0.gu[t_scale_iu] += ddtu * _mdl_y.weight1[idx] * help;
        }
      }
      if (_mdl_y.weight2[idx] != 0.0) {
        for (j = 0; j < _nx; j++)
          f0.gx[j] += dty * _mdl_y.weight2[idx] *
            2.0 * help * c.Jx[i + kk%spsk][j];
        if (_t_active) {
          f0.gx[t_scale_ix] += ddtx * _mdl_y.weight2[idx] * help*help;
          if (kk < _KK)
            f0.gu[t_scale_iu] += ddtu * _mdl_y.weight2[idx] * help*help;
        }
      }
      i += spsk;
    }

    // contribution of soft constraints
    if (_mdl_y_soft.active[idx] == 1) {
      if (kk < _KK) {
        help = u[upsk*(_nu+_nsu) + is + kk%spsk];
	f0.gu[upsk*(_nu+_nsu) + is + kk%spsk]
	  += dty * (_mdl_y_soft.weight1[idx]
		    + 2.0*_mdl_y_soft.weight2[idx]*help);
      }
      else {
        help = x[_nx + is + kk%spsk];
	f0.gx[_nx + is + kk%spsk]
	  += dty * (_mdl_y_soft.weight1[idx]
		    + 2.0*_mdl_y_soft.weight2[idx]*help);
        
      }
      if (_t_active) {
        f0.gx[t_scale_ix] += ddtx * (_mdl_y_soft.weight1[idx]*help
                                     + _mdl_y_soft.weight2[idx]*help*help);
        if (kk < _KK)
          f0.gu[t_scale_iu] += ddtu * (_mdl_y_soft.weight1[idx]*help
                                       + _mdl_y_soft.weight2[idx]*help*help);
      }
      is += spsk;
    }
  }

  // additional terms at initial time
  if (kk == 0) {
    for (i = spsk*(_nc+_nsc) + upsk*_nsuc, idx = 0; idx < _mdl_ny; idx++) {
      // contributions of final objective terms
      if (_mdl_y0.active[idx]) {
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
	i++;
      }
    }
  }
  // additional terms at final time
  if (kk == _KK) {
    for (i = _nc+_nsc, isf = 0, idx = 0; idx < _mdl_ny; idx++) {
      // contributions of final objective terms
      if (_mdl_yf.active[idx]) {
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
	i++;
      }
      // contribution of soft constraints
      if (_mdl_yf_soft.active[idx] == 1) {
        f0.gx[_nx + spsk*_ns + isf] += 
          (_mdl_yf_soft.weight1[idx]
           + 2.0*_mdl_yf_soft.weight2[idx]*x[_nx + spsk*_ns + isf]);
        isf++;
      }
    }
  }
}

//--------------------------------------------------------------------------
void Prg_DynamicOpt::consistic(int kk, double t,
			       const Omu_StateVec &x, const Omu_Vec &u,
			       Omu_DependentVec &xt)
{
  int i, idx;

  // initialize model inputs
  real_T *mdl_u = NULL;
  if (ssGetNumInputPorts(_SS) > 0) {
    if (ssGetInputPortRequiredContiguous(_SS, 0))
      mdl_u = (real_T *)ssGetInputPortRealSignal(_SS, 0);
    else
      mdl_u = (real_T *)*ssGetInputPortRealSignalPtrs(_SS, 0);
  }
  for (i = _mdl_nd, idx = 0; idx < _mdl_nu; idx++) {
    if (_mdl_u.active[idx])
      mdl_u[idx] = x[i++] * _mdl_u_nominal[idx];
    else if (kk == _KK && _KK > 0 && _mdl_u_order[idx] == 0)
      // hold last but one value in case of zero order hold
      mdl_u[idx] = _mdl_us[kk-1][idx];
    else
      mdl_u[idx] = _mdl_us[kk][idx];
  }

  // initialize scaled time 
  if (kk == 0 || _t_scale_idx < 0)
    _taus[kk] = ts(kk);
  else
    _taus[kk] = _taus[kk-1] +  _mdl_us[kk-1][_t_scale_idx] * (ts(kk) - ts(kk-1));

  // set simulation time
  ssSetT(_SS, _taus[kk]);

  // initialize model in first stage
  // This way model initial conditions get applied as functions of parameters.
  // Don't initialize if time changed, e.g. for subsequent simulation calls
  if ((_mdl_needs_init || kk == 0 && t == _t0_setup_model)
      && ssGetmdlInitializeConditions(_SS) != NULL) {
    // initialize model
    SMETHOD_CALL(mdlInitializeConditions, _SS);
    _mdl_needs_init = false;
  }

  // disable discrete sample times
  setSampleHit(false);
  // mark continuous sample time and process mdlUpdate in case of success
  if (setContinuousTask(true)) {
    // pass states from optimizer to model
    real_T *mdl_xd = ssGetDiscStates(_SS);
    real_T *mdl_xc = ssGetContStates(_SS);
    for (i = 0; i < _mdl_nd; i++)
      mdl_xd[i] = x[i] * _mdl_x_nominal[i];
    for (i = _mdl_nd; i < _mdl_nx; i++)
      mdl_xc[i - _mdl_nd] = x[_nu + i] * _mdl_x_nominal[i];

    // call mdlUpdate to get discrete events processed
    // Note: this is done once at the beginning of a sample interval;
    // no event processing takes place during the integration.
    if (ssGetmdlUpdate(_SS) != NULL) {
      // also call mdlOutputs as done by Simulink before each mdlUpdate
      SMETHOD_CALL2(mdlOutputs, _SS, 0); 
      SMETHOD_CALL2(mdlUpdate, _SS, 0);
    }

    // take over optimized control inputs from optimizer
    for (i = 0; i < _nu; i++)
      xt[_mdl_nd + i] = x[_mdl_nd + i];

    // read back states from model
    // Note: take optimized initial states from the optimizer
    // to prevent their overriding by the model, e.g. from parameters;
    // they are read once in Omu_Model::setup_model() instead.
    for (i = 0; i < _mdl_nd; i++) {
      if (kk == 0 && _mdl_x0_active[i])
        xt[i] = x[i];
      else
        xt[i] = mdl_xd[i] / _mdl_x_nominal[i];
    }
    for (i = _mdl_nd; i < _mdl_nx; i++) {
      if (kk == 0 && _mdl_x0_active[i])
        xt[_nu + i] = x[_nu + i];
      else
        xt[_nu + i] = mdl_xc[i - _mdl_nd] / _mdl_x_nominal[i];
    }
  }
  else {
    // take over states from optimizer
    for (i = 0; i < _nx; i++)
      xt[i] = x[i];
  }

  // take over slack variables from optimizer at final time
  if (kk == _KK) {
    for (i = 0; i < _ns+_nsf; i++)
      xt[_nx + i] = x[_nx + i];
  }
}

//--------------------------------------------------------------------------
void Prg_DynamicOpt::continuous(int kk, double t,
				const Omu_StateVec &x, const Omu_Vec &u,
				const Omu_StateVec &dx, Omu_DependentVec &F)
{
  int i, idx;
  int upsk = _multistage? 1: _KK;
  double rt = (t - ts(kk)) / (ts(kk+1) - ts(kk));
  double tscale = _t_active?
    (x[_mdl_nd + _t_scale_i] + u[_t_scale_i*upsk + kk%upsk])* _t_scale_nominal: 1.0;

  // initialize model inputs
  real_T *mdl_u = NULL;
  if (ssGetNumInputPorts(_SS) > 0) {
    if (ssGetInputPortRequiredContiguous(_SS, 0))
      mdl_u = (real_T *)ssGetInputPortRealSignal(_SS, 0);
    else
      mdl_u = (real_T *)*ssGetInputPortRealSignalPtrs(_SS, 0);
  }
  for (i = _mdl_nd, idx = 0; idx < _mdl_nu; idx++) {
    if (_mdl_u.active[idx])
      mdl_u[idx] = x[i++] * _mdl_u_nominal[idx];
    else if (_mdl_u_order[idx] == 0)
      mdl_u[idx] = _mdl_us[kk][idx];
    else
      mdl_u[idx] = _mdl_us[kk][idx] * (1 - rt) + _mdl_us[kk+1][idx] * rt;
  }

  // obtain scaled simulation time
  if (kk < _KK)
    _taus[kk+1] = _taus[kk] + tscale * (ts(kk+1) - ts(kk));
  double tau = _taus[kk] * (1 - rt) + _taus[kk+1] * rt;

  // set simulation time
  ssSetT(_SS, tau);

  // pass current states to model
  real_T *mdl_xd = ssGetDiscStates(_SS);
  real_T *mdl_xc = ssGetContStates(_SS);
  for (i = 0; i < _mdl_nd; i++)
    mdl_xd[i] = x[i] * _mdl_x_nominal[i];
  for (i = _mdl_nd; i < _mdl_nx; i++)
    mdl_xc[i - _mdl_nd] = x[_nu + i] * _mdl_x_nominal[i];

  // set model outputs before calculating derivatives
  // (note: this is required for sub-blocks providing inputs other blocks)
  // (furthermore a model may check for discrete changes in mdlOutputs)
  SMETHOD_CALL2(mdlOutputs, _SS, 0); 

  // evaluate continuous model equations
  if (_mdl_nx > _mdl_nd)
    SMETHOD_CALL(mdlDerivatives, _SS);

  // get model derivatives and change to residual form
  real_T *mdl_dx = ssGetdX(_SS);
  for (i = _mdl_nd; i < _mdl_nx; i++)
    F[_nu + i] = tscale*mdl_dx[i - _mdl_nd]/_mdl_x_nominal[i] - dx[_nu + i];

  // model equations for controlled inputs
  for (i = _mdl_nd, idx = 0; idx < _mdl_nu; idx++) {
    if (_mdl_u.active[idx]) {
      if (_mdl_u_order[idx] == 0)
	// zero order hold
	F[i] = 0.0 - dx[i];
      else
	// piecewise linear interpolation
	F[i] = tscale*u[(i-_mdl_nd)*upsk + kk%upsk]/_t_nominal - dx[i];
      i++;
    }
  }
  assert(i == _mdl_nd + _nu); // problem structure must not have changed

  // obtain Jacobians if required
  if (F.is_required_J())
    continuous_grds(kk, t, x, u, dx, F);
}

//--------------------------------------------------------------------------
void Prg_DynamicOpt::continuous_grds(int kk, double t,
				     const Omu_StateVec &x,
				     const Omu_Vec &u,
				     const Omu_StateVec &dx,
				     Omu_DependentVec &F)
{
  if (ssGetmdlJacobian(_SS) == NULL || _t_active) {
    // todo: use analytic Jacobian if _t_active
    // call predefined continuous_grds for numerical differentiation
    Omu_Program::continuous_grds(kk, t, x, u, dx, F, Omu_Dependent::WRT_x);
  }
  else {
    // exploit mdlJacobian
    int i, idx, j, jdx, rdx;
    int mdl_u_idx, mdl_x_idx;
    int mdl_nc = _mdl_nx - _mdl_nd; // number of continuous states
    real_T *pr = ssGetJacobianPr(_SS);
    int_T *ir = ssGetJacobianIr(_SS);
    int_T *jc = ssGetJacobianJc(_SS);

    mdlJacobian(_SS);

    // obtain Jx
    m_zero(F.Jx);
    // Jacobian wrt S-function states (ddxc/dxcxd)
    for (jdx = 0; jdx < _mdl_nx; jdx++) {
      mdl_x_idx = jdx < mdl_nc? _mdl_nd + jdx: jdx - mdl_nc;
      for (rdx = jc[jdx]; rdx < jc[jdx+1]; rdx++) {
	idx = ir[rdx];
	if (idx >= mdl_nc)
	  break;
	i = _mdl_nd + _nu + idx;
	j = jdx < mdl_nc? _mdl_nd + _nu + jdx: jdx - mdl_nc;
	F.Jx[i][j] = pr[rdx] /
	  _mdl_x_nominal[_mdl_nd + idx] * _mdl_x_nominal[mdl_x_idx];
      }
    }
    // Jacobian wrt S-function inputs (ddxc/du)
    // (note: these are states for the optimizer)
    for (j = _mdl_nd, jdx = _mdl_nx; jdx < _mdl_nx + _mdl_nu; jdx++) {
      mdl_u_idx = jdx - _mdl_nx;
      if (_mdl_u.active[mdl_u_idx]) {
	for (rdx = jc[jdx]; rdx < jc[jdx+1]; rdx++) {
	  idx = ir[rdx];
	  if (idx >= mdl_nc)
	    break;
	  i = _mdl_nd + _nu + idx;
	  F.Jx[i][j] = pr[rdx] /
	    _mdl_x_nominal[_mdl_nd + idx] * _mdl_u_nominal[mdl_u_idx];
	}
	j++;
      }
    }
  }

  // set F.Ju analytically
  // note: this can't generally be done in setup_struct as
  //   - if _t_active: F.Ju depends on tscale
  //   - if not _multistage: structure of F.Ju changes over sample periods
  if (_t_active || !_multistage) {
    int i, idx;
    int upsk = _multistage? 1: _KK;
    double tscale = _t_active?
      (x[_mdl_nd + _t_scale_i] + u[_t_scale_i*upsk + kk%upsk])* _t_scale_nominal: 1.0;
    m_zero(F.Ju);
    for (i = _mdl_nd, idx = 0; idx < _mdl_nu; idx++) {
      if (_mdl_u.active[idx]) {
        if (_mdl_u_order[idx] > 0) {
          // piecewise linear interpolation
          F.Ju[i][(i-_mdl_nd)*upsk + kk%upsk] = tscale/_t_nominal;
          if (_t_active)
            F.Ju[i][_t_scale_i*upsk + kk%upsk] =
              _t_scale_nominal * u[(i-_mdl_nd)*upsk + kk%upsk]/_t_nominal;
        }
        i++;
      }
    }
  }
}


//==========================================================================

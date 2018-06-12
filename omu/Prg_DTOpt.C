/*
 * Prg_DTOpt.C -- class definition
 *
 */

/*
    Copyright (C) 1997--2018  Ruediger Franke

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

#include "Prg_DTOpt.h"
#include "Hqp_omp.h"

#include <stdlib.h>

#include <If_Real.h>
#include <If_RealVec.h>
#include <If_RealMat.h>
#include <If_IntVec.h>
#include <If_Int.h>
#include <If_IntVec.h>
#include <If_Bool.h>
#include <If_Method.h>

// redefine assert to throw an error instead of aborting
#undef assert
#define assert(expr) if (!(expr)) m_error(E_INTERN, "assert(" #expr ")");

#define GET_SET_CB(vartype, prefix, name) \
  GET_CB(vartype, prefix, name), \
  IF_SET_CB(vartype, Prg_DTOpt, set_##name)

#define GET_CB(vartype, prefix, name) \
  prefix#name, \
  IF_GET_CB(vartype, Prg_DTOpt, name)

// Throw E_CONV as errors occuring during solution are generally
// due to bad values and need to be treated (esp. during stepsize check).
#undef SMETHOD_ERROR
#define SMETHOD_ERROR E_CONV

typedef If_Method<Prg_DTOpt> If_Cmd;

IF_CLASS_DEFINE("DTOpt", Prg_DTOpt, Hqp_SqpProgram);

inline double absmax(double a, double b)
{
  return fabs(a) > fabs(b)? a: b;
}

//--------------------------------------------------------------------------
Prg_DTOpt::Prg_DTOpt()
  : Hqp_Docp(omp_get_max_threads())
  , Omu_Model(ncpu())
{
  _t_scale_idx = -1;
  _t_active = 0;
  _t_scale_i = -1;
  _t_scale_nominal = 1.0;
  _stages_ok = false;
  _ad = true;
  _fscale = 1.0;

  _mdl_x0_active = iv_get(_mdl_nx);
  _mdl_u_order = iv_get(_mdl_nu);
  _mdl_u0_nfixed = iv_get(_mdl_nu);
  _mdl_u_decimation = iv_get(_mdl_nu);
  _mdl_u_periodic = iv_get(_mdl_nu);
  _mdl_x_periodic = iv_get(_mdl_nx);
  _t_nominal = 1.0;
  _mdl_y_order = iv_get(_mdl_ny);
  _mdl_y_bias = v_get(_mdl_ny);
  _mdl_y_lambda = v_get(_mdl_ny);
  _c_lambda = VNULL;
  v_set(_mdl_x0, 0.0);
  iv_set(_mdl_x0_active, 0);
  iv_set(_mdl_u_order, 1);
  iv_set(_mdl_u0_nfixed, 0);
  iv_set(_mdl_u_decimation, 1);
  iv_set(_mdl_u_periodic, 0);
  iv_set(_mdl_x_periodic, 0);
  iv_set(_mdl_y_order, 1);
  v_set(_mdl_y_bias, 0.0);
  v_set(_mdl_y_lambda, 0.0);

  // numbers of optimization variables
  _nu = 0;
  _ndu = 0;
  _nx = _ndu + _mdl_nx;
  _nc = 0;
  _ns = 0;
  _nsc = 0;
  _nsu = 0;
  _nsuc = 0;
  _nc0 = 0;
  _ncf = 0;
  _nsf = 0;
  _nscf = 0;

  _K = 0;
  _mdl_us = m_get(_K+1, _mdl_nu);
  _mdl_xs = m_get(_K+1, _mdl_nx);
  _mdl_ys = m_get(_K+1, _mdl_ny);
  _t0 = 0.0;
  _tf = 1.0;
  _ts = v_get(_K+1);
  _taus = v_get(_K+1);

  _ifList.append(new If_Cmd("prg_setup_stages",
                            &Prg_DTOpt::setup_stages, this));
  _ifList.append(new If_Int(GET_SET_CB(int, "prg_", K)));
  _ifList.append(new If_Real(GET_SET_CB(double, "prg_", t0)));
  _ifList.append(new If_Real(GET_SET_CB(double, "prg_", tf)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "prg_", ts)));
  _ifList.append(new If_Int(GET_SET_CB(int, "", mdl_t_scale_idx)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "prg_", taus)));
  _ifList.append(new If_Bool(GET_SET_CB(bool, "prg_", ad)));
  _ifList.append(new If_Real(GET_SET_CB(double, "prg_", fscale)));

  // redefine mdl_x0 in order to also consider mdl_xs
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_x0)));
  _ifList.append(new If_IntVec(GET_SET_CB(const IVECP, "", mdl_x0_active)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_x0_min)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, "", mdl_x0_max)));
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
  _ifList.append(new If_RealVec(GET_CB(const VECP, "", mdl_y_lambda)));
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
Prg_DTOpt::~Prg_DTOpt()
{
  v_free(_taus);
  v_free(_ts);
  m_free(_mdl_ys);
  m_free(_mdl_xs);
  m_free(_mdl_us);
  v_free(_mdl_y_bias);
  v_free(_mdl_y_lambda);
  iv_free(_mdl_y_order);
  iv_free(_mdl_x_periodic);
  iv_free(_mdl_u_periodic);
  iv_free(_mdl_u_decimation);
  iv_free(_mdl_u0_nfixed);
  iv_free(_mdl_u_order);
  iv_free(_mdl_x0_active);
}

//--------------------------------------------------------------------------
void Prg_DTOpt::setup_model()
{
  if (_mdl_logging >= If_LogInfo)
    If_Log("Info", "Prg_DTOpt::setup_model");

  // limit number of CPUs to number of stages
  if (ncpu() > _K + 1)
    set_ncpu(_K + 1);

  // load FMU or S-function
  Omu_Model::setup_model(_t0, ncpu());

  // no continuous states
  assert(_mdl_nd == _mdl_nx);

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
  iv_resize(_mdl_u_order, _mdl_nu);
  iv_resize(_mdl_u0_nfixed, _mdl_nu);
  iv_resize(_mdl_u_decimation, _mdl_nu);
  iv_resize(_mdl_u_periodic, _mdl_nu);
  iv_resize(_mdl_x_periodic, _mdl_nx);
  iv_resize(_mdl_y_order, _mdl_ny);
  v_resize(_mdl_y_bias, _mdl_ny);
  v_resize(_mdl_y_lambda, _mdl_ny);
  iv_set(_mdl_x0_active, 0);
  iv_set(_mdl_u_order, 1);
  iv_set(_mdl_u0_nfixed, 0);
  iv_set(_mdl_u_decimation, 1);
  iv_set(_mdl_u_periodic, 0);
  iv_set(_mdl_x_periodic, 0);
  iv_set(_mdl_y_order, 1);
  v_set(_mdl_y_bias, 0.0);
  v_set(_mdl_y_lambda, 0.0);
}

//--------------------------------------------------------------------------
void Prg_DTOpt::setup_horizon(int &k0, int &kf)
{
  // call setup_stages if this wasn't already done through the interface
  if (!_stages_ok)
    setup_stages();

  _stages_ok = false;	// for subsequent modifications

  k0 = 0;
  kf = _K;
}

//--------------------------------------------------------------------------
void Prg_DTOpt::setup_stages()
{
  int k, j;

  if (_mdl_logging >= If_LogInfo)
    If_Log("Info", "Prg_DTOpt::setup_stages");

  // setup model
  if (_mdl_needs_setup)
    setup_model();

  // get model structure
  if (_mdl_nx > _mdl_nd)
    m_error(E_FORMAT, "Prg_DTOpt::setup_stages: "
            "prg_name DTOpt requires model without continuous states");

  // setup optimization program
  v_resize(_ts, _K + 1);
  _ts[0] = _t0; // explicitly assign first value to avoid div. by zero for _K=0
  for (k = 1; k <= _K; k++)
    _ts[k] = _t0 + (double)k * (_tf - _t0) / (double)_K;

  m_resize(_mdl_us, _K+1, _mdl_nu);
  m_resize(_mdl_xs, _K+1, _mdl_nx);
  m_resize(_mdl_ys, _K+1, _mdl_ny);

  // setup _mdl_xs and mdl_us with start values from model
  set_mdl_x0(Omu_Model::_mdl_x_start);
  set_mdl_u0(Omu_Model::_mdl_u_start);

  // setup scaling of time
  v_copy(_ts, _taus);
  if (!_mdl_is_fmu && _t_scale_idx >= 0) {
    // skip initialization for FMU because its start value is used
    for (k = 0; k < _K; k++)
      _mdl_us[k][_t_scale_idx] = 1.0;
  }

  _stages_ok = true;
}

//--------------------------------------------------------------------------
void Prg_DTOpt::setup_vars(int k,
			   VECP x, VECP x_min, VECP x_max, IVECP x_int,
			   VECP u, VECP u_min, VECP u_max, IVECP u_int,
			   VECP c, VECP c_min, VECP c_max)
{
  int i, idx, isc, iscf, j;

  // complete general setup
  if (k == 0) {
    // take over possibly modified _mdl_p
    write_mx_args(_mdl_p);

    // obtain numbers of optimization variables (active model variables)
    _nu = 0;
    _ndu = 0;
    _nsu = 0;
    _nsuc = 0;
    for (idx = 0; idx < _mdl_nu; idx++) {
      if (idx == _t_scale_idx)
        _t_scale_i = _nu;
      if (_mdl_u.active[idx]) {
	_nu++;
        if (true || _mdl_der_u.min[idx] > -Inf || _mdl_der_u.max[idx] < Inf ||
            _mdl_der_u.weight1[idx] != 0.0 || _mdl_der_u.weight2[idx] != 0.0 ||
            (_mdl_der_u_soft.min[idx] > -Inf || _mdl_der_u_soft.max[idx] < Inf) 
            && (_mdl_der_u_soft.weight1[idx] != 0.0
                || _mdl_der_u_soft.weight2[idx] != 0.0) ||
            _mdl_u_periodic[idx]) {
          // introduce additional state for input as derivative needed
          // Note: always introduce state for now as this appears more reliable
          _mdl_der_u.active[idx] = 1;
          _ndu++;
        }
      }
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
	  m_error(E_FORMAT, "Prg_DTOpt::setup: "
		  "mdl_u_order[mdl_t_scale_idx] must be 0 (zero order hold)");
    }
    else
      _t_active = 0;

    _nx = _ndu + _mdl_nx;
    _nc = 0;
    _ns = 0;
    _nsc = 0;
    _nc0 = 0;
    _ncf = 0;
    _nsf = 0;
    _nscf = 0;
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
    _t_nominal = (ts(_K) - ts(0)) / _K;
  }

  // allocate states and controls for optimization
  int spsk, upsk;
  int nxk = k < _K? _nx: _nu - _ndu + _nx;
  spsk = 1;
  upsk = 1;
  if (k < _K) {
    alloc_vars(x, x_min, x_max, x_int, _nx);
    if (k == _K-1)
      // Note: allocate additional us for slacks at final time
      //       as the final states must have state equations
      alloc_vars(u, u_min, u_max, u_int, _nu + _nsu + _ns + _nsf);
    else
      alloc_vars(u, u_min, u_max, u_int, _nu + _nsu + _ns);
    if (k == 0)
      alloc_vars(c, c_min, c_max, IVNULL, _nc + _nsc + _nsuc + _nc0);
    else
      alloc_vars(c, c_min, c_max, IVNULL, _nc + _nsc + _nsuc);
  }
  else {
    // at final time allocate additional final states for slack variables
    // as no control parameters allowed here
    // (propagate us and xs to final stage as they are accessed in update)
    alloc_vars(x, x_min, x_max, x_int, nxk + _ns + _nsf);
    if (_K > 0)
      alloc_vars(c, c_min, c_max, IVNULL, _nc + _nsc + _ncf + _nscf);
    else
      alloc_vars(c, c_min, c_max, IVNULL, _nc + _nsc + _nc0 + _ncf + _nscf);
  }

  // setup states and control inputs
  for (i = _mdl_nd, j = 0, idx = 0; idx < _mdl_nu; idx++) {
    if (k == 0 && _mdl_u_order[idx] < 0 || _mdl_u_order[idx] > 1)
      m_error(E_FORMAT, "Prg_DTOpt::setup: "
              "mdl_u_order must be 0 or 1");
    if (_mdl_u.active[idx]) {
      if (_mdl_der_u.active[idx] || _K == 0) {
        x[i] = _mdl_us[k][idx] / _mdl_u_nominal[idx];
        if (k == 0 && _mdl_u0_nfixed[idx] > 1 - _mdl_u_order[idx])
          // fixed initial value
          x_min[i] = x_max[i] = x[i]; 
        else if (k == 0 && _mdl_u0_nfixed[idx] > 0 && _mdl_u_order[idx] == 0) {
          // apply rate of change bounds to initial step for zero order hold
          if (_mdl_der_u.min[idx] > -Inf)
            x_min[i] = x[i]
              + _mdl_der_u.min[idx] / _mdl_u_nominal[idx]*_t_nominal;
          if (_mdl_der_u.max[idx] < Inf)
            x_max[i] = x[i]
              + _mdl_der_u.max[idx] / _mdl_u_nominal[idx]*_t_nominal;
        }
        else if (k >= _mdl_u0_nfixed[idx] + _mdl_u_order[idx] - 1) {
          // control bounds
          if (k==0 && _mdl_u_periodic[idx])
            // indicate periodic state to Hqp_Docp
            x_min[i] = Inf;
          else if (k == 0 && _mdl_u0.min[idx] > _mdl_u.min[idx])
            // at k==0, use more restrictive bound of u0_min/max and u_min/max
            x_min[i] = _mdl_u0.min[idx] / _mdl_u_nominal[idx];
          else if (k == _K && _mdl_uf.min[idx] > _mdl_u.min[idx])
            // at k==_K, use more restrictive bound of uf_min/max and u_min/max
            x_min[i] = _mdl_uf.min[idx] / _mdl_u_nominal[idx];
          else if (_mdl_u.min[idx] > -Inf)
            x_min[i] = _mdl_u.min[idx] / _mdl_u_nominal[idx];
          if (k==0 && _mdl_u_periodic[idx])
            x_max[i] = Inf;
          else if (k == 0 && _mdl_u0.max[idx] < _mdl_u.max[idx])
            x_max[i] = _mdl_u0.max[idx] / _mdl_u_nominal[idx];
          else if (k == _K && _mdl_uf.max[idx] < _mdl_u.max[idx])
            x_max[i] = _mdl_uf.max[idx] / _mdl_u_nominal[idx];
          else if (_mdl_u.max[idx] < Inf)
            x_max[i] = _mdl_u.max[idx] / _mdl_u_nominal[idx];
          // integer variables
          x_int[i] = _mdl_u.integer[idx];
        }
        i++;
      }
      else if (k < _K) {
        u[j] = _mdl_us[k][idx] / _mdl_u_nominal[idx];
        if (k == 0 && _mdl_u0_nfixed[idx] > 1 - _mdl_u_order[idx])
          // fixed initial value
          u_min[j] = u_max[idx] = u[idx]; 
        else if (k >= _mdl_u0_nfixed[idx] + _mdl_u_order[idx] - 1) {
          // control bounds
          if (k == 0 && _mdl_u0.min[idx] > _mdl_u.min[idx])
            // at k==0, use more restrictive bound of u0_min/max and u_min/max
            u_min[j] = _mdl_u0.min[idx] / _mdl_u_nominal[idx];
          else if (k == _K && _mdl_uf.min[idx] > _mdl_u.min[idx])
            // at k==_K, use more restrictive bound of uf_min/max and u_min/max
            u_min[j] = _mdl_uf.min[idx] / _mdl_u_nominal[idx];
          else if (_mdl_u.min[idx] > -Inf)
            u_min[j] = _mdl_u.min[idx] / _mdl_u_nominal[idx];
          else if (k == 0 && _mdl_u0.max[idx] < _mdl_u.max[idx])
            u_max[j] = _mdl_u0.max[idx] / _mdl_u_nominal[idx];
          else if (k == _K && _mdl_uf.max[idx] < _mdl_u.max[idx])
            u_max[j] = _mdl_uf.max[idx] / _mdl_u_nominal[idx];
          else if (_mdl_u.max[idx] < Inf)
            u_max[j] = _mdl_u.max[idx] / _mdl_u_nominal[idx];
          // integer variables
          u_int[j] = _mdl_u.integer[idx];
        }
        j++;
      }
    }
  }
  for (idx = 0, i = 0; i < _nx; idx++, i++) {
    if (i == _mdl_nd && _ndu > 0) {
      // skip control states
      i += _ndu - 1;
      idx--;
      continue;
    }
    x[i] = _mdl_xs[k][idx] / _mdl_x_nominal[idx];
    if (k == 0) {
      if (_mdl_x_periodic[idx])
        // indicate periodic state to Hqp_Docp
        x_min[i] = x_max[i] = Inf;
      else if (_mdl_x0_active[idx]) {
        // use the more restrictive bound of x0_min/max and x_min/max
        if (_mdl_x0.min[idx] > _mdl_x.min[idx])
          x_min[i] = _mdl_x0.min[idx] / _mdl_x_nominal[idx];
        else
          x_min[i] = _mdl_x.min[idx] / _mdl_x_nominal[idx];
        if (_mdl_x0.max[idx] < _mdl_x.max[idx])
          x_max[i] = _mdl_x0.max[idx] / _mdl_x_nominal[idx];
        else
          x_max[i] = _mdl_x.max[idx] / _mdl_x_nominal[idx];
        // integer variables
        x_int[i] = _mdl_x.integer[idx];
      }
      else
	x_min[i] = x_max[i] = x[i];
    }
    else {
      x_min[i] = _mdl_x.min[idx] / _mdl_x_nominal[idx];
      x_max[i] = _mdl_x.max[idx] / _mdl_x_nominal[idx];
      // integer variables
      x_int[i] = _mdl_x.integer[idx];
    }
  }

  // setup rate of change bounds for control inputs
  for (i = 0, isc = 0, idx = 0; idx < _mdl_nu; idx++) {
    if (_mdl_der_u.active[idx]) {
      if (k < _K) {
	// rate of change bounds
	if (_mdl_der_u.min[idx] > -Inf) {
	  for (j = 0; j < upsk; j++)
	    u_min[i+j] =
	      _mdl_der_u.min[idx] / _mdl_u_nominal[idx]*_t_nominal;
	}
	if (_mdl_der_u.max[idx] < Inf) {
	  for (j = 0; j < upsk; j++)
	    u_max[i+j] =
	      _mdl_der_u.max[idx] / _mdl_u_nominal[idx]*_t_nominal;
	}
        // rate of change soft bounds
        if (_mdl_der_u_soft.active[idx]) {
          if (_mdl_der_u_soft.min[idx] > -Inf) {
            for (j = 0; j < upsk; j++)
              c_min[spsk*(_nc+_nsc) + isc + j] =
                _mdl_der_u_soft.min[idx] / _mdl_u_nominal[idx]*_t_nominal;
            isc += upsk;
          }
          if (_mdl_der_u_soft.max[idx] < Inf) {
            for (j = 0; j < upsk; j++)
              c_max[spsk*(_nc+_nsc) + isc + j] =
                _mdl_der_u_soft.max[idx] / _mdl_u_nominal[idx]*_t_nominal;
            isc += upsk;
          }
        }
	// initial guess
        u[i] =
          (_mdl_us[k+1][idx] - _mdl_us[k][idx])
          / (ts(k+1) - ts(k))
          / _mdl_u_nominal[idx]*_t_nominal;
	// override rate of change bounds for fixed controls
	if (k < _mdl_u0_nfixed[idx] + _mdl_u_order[idx] - 2) {
	  u_min[i] = u_max[i] = u[i];
	}
	// override rate of change bounds for joined controls
	if (_mdl_u_decimation[idx] == 0 ||
            (k + 1) % _mdl_u_decimation[idx] != 0) {
          u_min[i] = u_max[i] = 0.0;
        }
      }
      // for zero order hold: constraint u parameter of last interval to zero
      // as control must not change at final time, compared to last interval
      if (k == _K-1 && _mdl_u_order[idx] == 0) {
        u_min[i] = u_max[i] = u[i] = 0.0;
      }
      i += upsk;
    }
    else if (_mdl_u.active[idx])
      i += upsk;
  }

  // setup constraints on model outputs
  for (i = 0, isc = 0, idx = 0; idx < _mdl_ny; idx++) {
    if (_mdl_y.min[idx] > -Inf)
      for (j = 0; j < spsk; j++)
	c_min[i+j] = _mdl_y.min[idx] / _mdl_y_nominal[idx];
    if (_mdl_y.max[idx] < Inf)
      for (j = 0; j < spsk; j++)
	c_max[i+j] = _mdl_y.max[idx] / _mdl_y_nominal[idx];
    if (_mdl_y.active[idx])
      i += spsk;
    if (_mdl_y_soft.active[idx]) {
      if (_mdl_y_soft.min[idx] > -Inf) {
        for (j = 0; j < spsk; j++)
          c_min[spsk*_nc + isc + j] = _mdl_y_soft.min[idx] / _mdl_y_nominal[idx];
        isc += spsk;
      }
      if (_mdl_y_soft.max[idx] < Inf) {
        for (j = 0; j < spsk; j++)
          c_max[spsk*_nc + isc + j] = _mdl_y_soft.max[idx] / _mdl_y_nominal[idx];
        isc += spsk;
      }
    }
  }

  // setup slack variables for soft constraints
  if (k < _K)
    for (i = upsk*_nu; i < upsk*(_nu+_nsu) + spsk*_ns; i++) {
      u_min[i] = 0.0;
      u[i] = 0.0;
    }
  else
    for (i = nxk; i < nxk + spsk*_ns + _nsf; i++) {
      x_min[i] = 0.0;
      x[i] = 0.0;
    }

  // setup state and output constraints at initial time
  if (k == 0) {
    for (idx = 0; idx < _mdl_ny; idx++) {
      if (_mdl_y0.active[idx]) {
        if (_mdl_y0.min[idx] > -Inf)
          c_min[i] = _mdl_y0.min[idx] / _mdl_y_nominal[idx];
        if (_mdl_y0.max[idx] < Inf)
          c_max[i] = _mdl_y0.max[idx] / _mdl_y_nominal[idx];
	i++;
      }
    }
  }

  // setup output constraints at final time
  if (k == _K) {
    for (i = spsk*(_nc+_nsc), iscf = 0, idx = 0; idx < _mdl_ny; idx++) {
      if (_mdl_yf.active[idx]) {
        if (_mdl_yf.min[idx] > -Inf)
          c_min[i] = _mdl_yf.min[idx] / _mdl_y_nominal[idx];
        if (_mdl_yf.max[idx] < Inf)
          c_max[i] = _mdl_yf.max[idx] / _mdl_y_nominal[idx];
        i++;
      }
      // soft constraints at final time
      if (_mdl_yf_soft.active[idx]) {
        if (_mdl_yf_soft.min[idx] > -Inf) {
          c_min[spsk*(_nc+_nsc) + _ncf + iscf] =
            _mdl_yf_soft.min[idx] / _mdl_y_nominal[idx];
          iscf++;
        }
        if (_mdl_yf_soft.max[idx] < Inf) {
          c_max[spsk*(_nc+_nsc) + _ncf + iscf] =
            _mdl_yf_soft.max[idx] / _mdl_y_nominal[idx];
          iscf++;
        }
      }
    }
  }
}

//--------------------------------------------------------------------------
void Prg_DTOpt::setup_struct(int k, const VECP x, const VECP u,
			     MATP fx, MATP fu, IVECP f_lin,
			     VECP f0x, VECP f0u, int &f0_lin,
			     MATP cx, MATP cu, IVECP c_lin,
			     MATP Lxx, MATP Luu, MATP Lxu)
{
  int i, j, idx, jdx, offs;
  SimStruct *S = _SS[0];

  if (_mdl_is_fmu && _mdl_jac && ssGetmdlJacobian(S) != NULL) {
    if (k == 0) {
      for (idx = 0; idx < _mdl_ny; idx++) {
        _mdl_jac_y_active[idx] = _mdl_y.active[idx] || _mdl_y_soft.active[idx]
          || _mdl_y0.active[idx] || _mdl_yf.active[idx] || _mdl_yf_soft.active[idx];
      }
      iv_copy(_mdl_u.active, _mdl_jac_u_active);
      setup_jac();
    }

    real_T *pr = ssGetJacobianPr(S);
    int i_end = ssGetJacobianJc(S)[_mdl_nx + _mdl_nu];
    for (i = 0; i < i_end; i++)
      *pr++ = 1.0;
    fetch_jac(S, k, 1.0, x, u, fx, fu, cx, cu);
  }
}

//--------------------------------------------------------------------------
void Prg_DTOpt::init_simulation(int k,
				VECP x, VECP u)
{
  int i, idx;

  if (_mdl_logging >= If_LogInfo)
    If_Log("Info", "Prg_DTOpt::init_simulation at k = %d", k);

  // initial states
  if (k == 0) {
    for (i = _mdl_nd, idx = 0; idx < _mdl_nu; idx++) {
      if (_mdl_der_u.active[idx] || k == _K && _mdl_u.active[idx])
        x[i++] = _mdl_us[k][idx] / _mdl_u_nominal[idx];
    }
    for (idx = 0, i = 0; i < _nx; idx++, i++) {
      if (i == _mdl_nd && _ndu > 0) {
        // skip control states
        i += _ndu - 1;
        idx--;
        continue;
      }
      x[i] = _mdl_xs[k][idx] / _mdl_x_nominal[idx];
    }
  }
}

//--------------------------------------------------------------------------
void Prg_DTOpt::update_vals(int k, const VECP x, const VECP u,
			    VECP f, Real &f0, VECP c)
{
  int i, j, idx, is, isc, isf, iscf;
  int nxk = k < _K? _nx: _nu - _ndu + _nx;
  int spsk = 1;
  int upsk = 1;
  double tscale_1 = _t_active? x[_mdl_nd + _t_scale_i]* _t_scale_nominal: 1.0;
  double tscale = _t_active && k < _K?
    (x[_mdl_nd + _t_scale_i] + u[_t_scale_i*upsk + k%upsk])* _t_scale_nominal: tscale_1;
  double help;
  int tn = omp_get_thread_num();
  SimStruct *S = _SS[tn];
  bool with_lambda = (k == 0 && _c_lambda != VNULL);

  if (_mdl_logging >= If_LogInfo)
    If_Log("Info", "Prg_DTOpt::update_vals at k = %d, tn = %d", k, tn);

  if (with_lambda)
    v_zero(_mdl_y_lambda);

  // initialize model inputs
  real_T *mdl_u = NULL;
  if (ssGetNumInputPorts(S) > 0) {
    if (ssGetInputPortRequiredContiguous(S, 0))
      mdl_u = (real_T *)ssGetInputPortRealSignal(S, 0);
    else
      mdl_u = (real_T *)*ssGetInputPortRealSignalPtrs(S, 0);
  }
  for (i = _mdl_nd, j = 0, idx = 0; idx < _mdl_nu; idx++) {
    if (_mdl_u.active[idx]) {
      if (_mdl_der_u.active[idx] || k == _K)
        mdl_u[idx] = x[i++] * _mdl_u_nominal[idx];
      else
        mdl_u[idx] = u[j++] * _mdl_u_nominal[idx];
      _mdl_us[k][idx] = mdl_u[idx];
    }
    else if (k == _K && _K > 0 && _mdl_u_order[idx] == 0)
      // hold last but one value in case of zero order hold
      mdl_u[idx] = _mdl_us[k-1][idx];
    else
      mdl_u[idx] = _mdl_us[k][idx];
  }

  // initialize scaled time 
  if (k == 0 || _t_scale_idx < 0)
    _taus[k] = ts(k);
  else
    _taus[k] = _taus[k-1] +  _mdl_us[k-1][_t_scale_idx] * (ts(k) - ts(k-1));

  // set simulation time
  ssSetT(S, _taus[k]);
  if (k < _K) {
    if (_t_scale_idx < 0)
      setSampleTime(S, ts(k+1) - ts(k));
    else
      setSampleTime(S, mdl_u[_t_scale_idx] * (ts(k+1) - ts(k)));
  }

  // pass current states to model
  real_T *mdl_xd = ssGetDiscStates(S);
  real_T *mdl_xc = ssGetContStates(S);
  bool needs_update = true;
  for (i = 0; i < _mdl_nd; i++) {
    mdl_xd[i] = x[i] * _mdl_x_nominal[i];
  }
  for (i = _mdl_nd; i < _mdl_nx; i++) {
    mdl_xc[i - _mdl_nd] = x[_nu + i] * _mdl_x_nominal[i];
  }

  // initialize model in first stage
  // Don't initialize if time changed, e.g. for subsequent simulation calls
  if ((_mdl_needs_init[tn] || k == 0 && ts(k) == _t0_setup_model)
      && ssGetmdlInitializeConditions(S) != NULL) {
    // initialize model
    SMETHOD_CALL(mdlInitializeConditions, S);
    _mdl_needs_init[tn] = 0;
    // call mdlOutputs because implementation may be delayed
    // call mdlUpdate and disable continuous task to trigger initial clock
    if (_mdl_is_fmu) {
      setSampleHit(S, true);
      setContinuousTask(S, false);
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
    // take active initial states from solver
    for (i = 0; i < _mdl_nd; i++) {
      if (_mdl_x0_active[i] || k > 0) {
        mdl_xd[i] = x[i] * _mdl_x_nominal[i];
        needs_update = true;
      }
    }
    for (i = _mdl_nd; i < _mdl_nx; i++) {
      if (_mdl_x0_active[i] || k > 0) {
        mdl_xc[i - _mdl_nd] = x[_nu + i] * _mdl_x_nominal[i];
        needs_update = true;
      }
    }
  }

  // call mdlOutputs/mdlUpdate to get outputs for current states
  // enable continuous task to not update states here
  if (needs_update) {
    if (_mdl_is_fmu) {
      setSampleHit(S, true);
      setContinuousTask(S, true);
    }
    // also call mdlOutputs as done by Simulink before each mdlUpdate
    SMETHOD_CALL2(mdlOutputs, S, 0);
    if (ssGetmdlUpdate(S) != NULL) {
      SMETHOD_CALL2(mdlUpdate, S, 0);
    }
    if (_mdl_is_fmu) {
      setSampleHit(S, false);
    }
  }

  // obtain model outputs
  real_T *mdl_y = NULL;
  if (ssGetNumOutputPorts(S) > 0)
    mdl_y = ssGetOutputPortRealSignal(S, 0);

  // correct model outputs with bias
  for (idx = 0; idx < _mdl_ny; idx++)
    mdl_y[idx] += _mdl_y_bias[idx];

  // utilize model outputs for constraints and optimization objective
  f0 = 0.0;

  double dt; 	// factor for linear interpol. (trapezoidal integration rule)
  double dt0; 	// factor for zero order hold
  double dtu; 	// factor used for controlled inputs in objective function
  double dty; 	// factor used for outputs in objective function
  if (k == 0) {
    if (_K > 0) {
      dt = 0.5 * (ts(k+1) - ts(k)) * tscale;
      dt0 = (ts(k+1) - ts(k)) * tscale;
    } else {
      // steady-state problem
      dt = dt0 = 1.0;
    }
  } else if (k == _K) {
    dt = 0.5 * (ts(k) - ts(k-1)) * tscale_1;
    dt0 = 0.0;
  } else {
    dt = 0.5 * ((ts(k) - ts(k-1)) * tscale_1 + (ts(k+1) - ts(k)) * tscale);
    dt0 = (ts(k+1) - ts(k)) * tscale;
  }

  // contribution of controlled inputs
  for (i = 0, is = 0, isc = 0, idx = 0; idx < _mdl_nu; idx++) {
    if (_mdl_u.active[idx]) {
      // control inputs
      dtu = (_mdl_u_order[idx] == 0)? dt0: dt;
      help = (_mdl_us[k][idx] - _mdl_u.ref[idx]) / _mdl_u_nominal[idx];
      f0 += dtu * _mdl_u.weight1[idx] * help;
      f0 += dtu * _mdl_u.weight2[idx] * help*help;
      // rates of change
      if (k < _K) {
        _mdl_der_u[idx] = u[i*upsk + k%upsk]*_mdl_u_nominal[idx]/_t_nominal;
	help = u[i*upsk + k%upsk]/_t_nominal - _mdl_der_u.ref[idx]
	  / _mdl_u_nominal[idx];
	f0 += dt0 * _mdl_der_u.weight1[idx] * help;
	f0 += dt0 * _mdl_der_u.weight2[idx] * help*help;
        // soft constraints on rates of change
        if (_mdl_der_u_soft.active[idx]) {
          help = u[upsk*_nu + is + k%upsk]/_t_nominal;
          f0 += dt0 * (_mdl_der_u_soft.weight1[idx]*help
                       + _mdl_der_u_soft.weight2[idx]*help*help);
          if (_mdl_der_u_soft.min[idx] > -Inf) {
            c[spsk*(_nc+_nsc) + isc + k%upsk] = u[i*upsk + k%upsk]/_t_nominal + help;
            isc += upsk;
          }
          if (_mdl_der_u_soft.max[idx] < Inf) {
            c[spsk*(_nc+_nsc) + isc + k%upsk] = u[i*upsk + k%upsk]/_t_nominal - help;
            isc += upsk;
          }
          is += upsk;
        }
      }
      i++;
    }
  }
  // contribution of active outputs
  for (i = 0, is = 0, isc = 0, idx = 0; idx < _mdl_ny; idx++) {
    dty = (_mdl_y_order[idx] == 0)? dt0: dt;
    // assign used outputs to constraints
    if (_mdl_y.active[idx]) {
      c[i + k%spsk] = mdl_y[idx] / _mdl_y_nominal[idx];
      if (with_lambda)
        _mdl_y_lambda[idx] = absmax(_c_lambda[i + k%spsk], _mdl_y_lambda[idx]);
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
      if (k < _K)
	help = u[upsk*(_nu+_nsu) + is + k%spsk];
      else
	help = x[nxk + is + k%spsk];

      f0 += dty * (_mdl_y_soft.weight1[idx]*help
		   + _mdl_y_soft.weight2[idx]*help*help);

      if (_mdl_y_soft.min[idx] > -Inf) {
	c[spsk*_nc + isc + k%spsk] = mdl_y[idx] / _mdl_y_nominal[idx] + help;
	if (with_lambda)
	  _mdl_y_lambda[idx] = absmax(_c_lambda[spsk*_nc + isc + k%spsk], _mdl_y_lambda[idx]);
	isc += spsk;
      }
      if (_mdl_y_soft.max[idx] < Inf) {
	c[spsk*_nc + isc + k%spsk] = mdl_y[idx] / _mdl_y_nominal[idx] - help;
	if (with_lambda)
	  _mdl_y_lambda[idx] = absmax(_c_lambda[spsk*_nc + isc + k%spsk], _mdl_y_lambda[idx]);
	isc += spsk;
      }
      is += spsk;
    }
  }

  // additional terms at initial time
  if (k == 0) {
    i = spsk*(_nc+_nsc) + upsk*_nsuc;
    // constraints on initial outputs
    for (idx = 0; idx < _mdl_ny; idx++) {
      // assign used outputs to constraints
      if (_mdl_y0.active[idx]) {
	c[i] = mdl_y[idx] / _mdl_y_nominal[idx];
	if (with_lambda)
	  _mdl_y_lambda[idx] = absmax(_c_lambda[i], _mdl_y_lambda[idx]);
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
  if (k == _K) {
    for (i = _nc+_nsc, isf = 0, iscf = 0, idx = 0; idx < _mdl_ny; idx++) {
      // assign used outputs to constraints
      if (_mdl_yf.active[idx]) {
	c[i] = mdl_y[idx] / _mdl_y_nominal[idx];
	if (with_lambda)
	  _mdl_y_lambda[idx] = absmax(_c_lambda[i], _mdl_y_lambda[idx]);
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
        help = x[nxk + spsk*_ns + isf];

        f0 += _mdl_yf_soft.weight1[idx]*help
                + _mdl_yf_soft.weight2[idx]*help*help;

        if (_mdl_yf_soft.min[idx] > -Inf) {
          c[spsk*(_nc+_nsc) + _ncf + iscf] = mdl_y[idx] / _mdl_y_nominal[idx] + help;
          if (with_lambda)
            _mdl_y_lambda[idx] = absmax(_c_lambda[spsk*(_nc+_nsc) + _ncf + iscf], _mdl_y_lambda[idx]);
          iscf++;
        }
        if (_mdl_yf_soft.max[idx] < Inf) {
          c[spsk*(_nc+_nsc) + _ncf + iscf] = mdl_y[idx] / _mdl_y_nominal[idx] - help;
          if (with_lambda)
            _mdl_y_lambda[idx] = absmax(_c_lambda[spsk*(_nc+_nsc) + _ncf + iscf], _mdl_y_lambda[idx]);
          iscf++;
        }
        isf++;
      }
    }
  }

  // store values of model states
  for (i = 0; i < _mdl_nd; i++)
    _mdl_xs[k][i] = value(mdl_xd[i]);
  for (; i < _mdl_nx; i++)
    _mdl_xs[k][i] = value(mdl_xc[i - _mdl_nd]);

  // store values of model outputs
  for (i = 0; i < _mdl_ny; i++)
    _mdl_ys[k][i] = value(mdl_y[i]);

  // store initial and final values
  if (k == 0) {
    for (i = 0; i < _mdl_nd; i++)
      _mdl_x0[i] = value(mdl_xd[i]);
    for (; i < _mdl_nx; i++)
      _mdl_x0[i] = value(mdl_xc[i - _mdl_nd]);
    for (i = 0; i < _mdl_nu; i++)
      _mdl_u0[i] = value(mdl_u[i]);
    for (i = 0; i < _mdl_ny; i++)
      _mdl_y0[i] = value(mdl_y[i]);
  }
  else if (k == _K) {
    for (i = 0; i < _mdl_nu; i++)
      _mdl_uf[i] = value(mdl_u[i]);
    for (i = 0; i < _mdl_ny; i++)
      _mdl_yf[i] = value(mdl_y[i]);
  }

  // junction conditions for subsequent stage
  if (k < _K) {
    if (_mdl_nd > 0) {
      // call mdlUpdate to get discrete events processed
      if (ssGetmdlUpdate(S) != NULL) {
        setContinuousTask(S, false);
        setSampleHit(S, true);
        if (_mdl_is_fmu) {
          // obtain discrete states at end of sample interval
          ssSetT(S, _taus[k + 1]);
          SMETHOD_CALL2(mdlOutputs, S, 0);
        }
        SMETHOD_CALL2(mdlUpdate, S, 0);
        setSampleHit(S, false);
        setContinuousTask(S, true);
      }
      // read discrete states from model
      for (i = 0; i < _mdl_nd; i++) {
        f[i] = mdl_xd[i] / _mdl_x_nominal[i];
      }
    }
    // update controlled inputs from optimizer
    for (i = _mdl_nd, idx = 0; idx < _mdl_nu; idx++) {
      if (_mdl_der_u.active[idx]) {
	if (true || _mdl_u_order[idx] == 0 && (k+1)%spsk == 0) {
	  // zero order hold at end of stage:
	  // apply step in u for subsequent stage
          help = ((ts(k+1) - ts(k-spsk/upsk+1))/_t_nominal
                  * u[(i-_mdl_nd)*upsk + k%upsk]);
          if (_t_active && i-_mdl_nd != _t_scale_i)
            f[i] = x[i] + tscale * help;
          else
            f[i] = x[i] + help;
        }
	//else
	  // piecewise linear interpolation or within stage with zoh
	  //f[i] = xf[i];
	i++;
      }
      else if (_mdl_u.active[idx] && k == _K - 1) {
        // propagate control inputs to final time
        f[i] = mdl_u[idx];
        i++;
      }
    }
    // problem structure must not have changed
    assert(i == _mdl_nd + _ndu || k == _K - 1 && i == _mdl_nd + _nu);
    // re-use slack variables of last but one sample period at final time
    // note: this is only because no control parameters allowed at final time
    // and because states need to be defined with state equations
    if (k == _K-1) {
      int nxK = _nx + _nu - _ndu; // nxk of final stage
      for (; i < nxK + _ns; i++)
	f[i] = u[upsk*(_nu+_nsu) + (i-nxK+1)*spsk-1];
      for (; i < nxK + _ns + _nsf; i++)
	f[i] = u[upsk*(_nu+_nsu) + spsk*_ns + i - nxK - _ns];
    }
  }
  f0 *= _fscale;
}

//--------------------------------------------------------------------------
void Prg_DTOpt::update_stage(int k, const VECP x, const VECP u,
			     VECP f, Real &f0, VECP c,
			     MATP fx, MATP fu, VECP f0x, VECP f0u,
			     MATP cx, MATP cu,
			     const VECP rf, const VECP rc,
			     MATP Lxx, MATP Luu, MATP Lxu)
{
  int i, iu, idx, ii, iidx0, iidx, isc, iscdx, is, j, ju, jdx, rdx;
  int isf, iscf, iscfdx, ixf;
  int nxk = k < _K? _nx: _nu - _ndu + _nx;
  int spsk = 1;
  int upsk = 1;
  double tscale_1 = _t_active? x[_mdl_nd + _t_scale_i]* _t_scale_nominal: 1.0;
  double tscale = _t_active && k < _K?
    (x[_mdl_nd + _t_scale_i] + u[_t_scale_i*upsk + k%upsk])* _t_scale_nominal: tscale_1;
  int t_scale_iu = _t_scale_i*upsk + k%upsk;
  int t_scale_ix = _mdl_nd + _t_scale_i;
  double help;
  int analyticJac = 0;
  int tn = omp_get_thread_num();
  SimStruct *S = _SS[tn];

  if (_mdl_logging >= If_LogInfo)
    If_Log("Info", "Prg_DTOpt::update_stage at k = %d, tn = %d", k, tn);

  if (!_ad || !_mdl_jac || ssGetmdlJacobian(S) == NULL || _t_active) {
    // todo: use analytic Jacobian if _t_active
    // call predefined update_stage for numerical differentiation
    Hqp_Docp::update_stage(k, x, u, f, f0, c, fx, fu, f0x, f0u, cx, cu,
                           rf, rc, Lxx, Luu, Lxu);
    if (!_ad)
      return;
  }
  else {
    if (k == 0)
      _c_lambda = rc;
    update_vals(k, x, u, f, f0, c);
    if (k == 0)
      _c_lambda = VNULL;
    analyticJac = 1;

    // restore states after mdlUpdate has been processed
    // call mdlOutputs
    real_T *mdl_xd = ssGetDiscStates(S);
    real_T *mdl_xc = ssGetContStates(S);
    for (i = 0; i < _mdl_nd; i++)
      mdl_xd[i] = x[i] * _mdl_x_nominal[i];
    for (i = _mdl_nd; i < _mdl_nx; i++)
      mdl_xc[i - _mdl_nd] = x[_nu + i] * _mdl_x_nominal[i];
    SMETHOD_CALL2(mdlOutputs, S, 0); 

    SMETHOD_CALL(mdlJacobian, S);

    fetch_jac(S, k, tscale, x, u, fx, fu, cx, cu);
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
  if (k == 0) {
    if (_K > 0) {
      dt = 0.5 * (ts(k+1) - ts(k)) * tscale;
      ddtx = 0.5 * (ts(k+1) - ts(k)) * _t_scale_nominal;
      ddtu = 0.5 * (ts(k+1) - ts(k)) * _t_scale_nominal;
      dt0 = (ts(k+1) - ts(k)) * tscale;
      ddt0 = (ts(k+1) - ts(k)) * _t_scale_nominal;
    } else {
      // steady-state problem
      dt = dt0 = 1.0;
      ddtx = 0.0;
      ddtu = 0.0;
      ddt0 = 0.0;
    }
  } else if (k == _K) {
    dt = 0.5 * (ts(k) - ts(k-1)) * tscale_1;
    ddtx = 0.5 * (ts(k) - ts(k-1)) * _t_scale_nominal;
    ddtu = 0.0;
    dt0 = 0.0;
    ddt0 = 0.0;
  } else {
    dt = 0.5 * ((ts(k) - ts(k-1)) * tscale_1 + (ts(k+1) - ts(k)) * tscale);
    ddtx = 0.5 * (ts(k+1) - ts(k-1)) * _t_scale_nominal;
    ddtu = 0.5 * (ts(k+1) - ts(k)) * _t_scale_nominal;
    dt0 = (ts(k+1) - ts(k)) * tscale;
    ddt0 = (ts(k+1) - ts(k)) * _t_scale_nominal;
  }

  v_zero(f0x);
  v_zero(f0u);
  // controlled inputs
  for (i = _mdl_nd, iu = 0, is = 0, isc = 0, idx = 0; idx < _mdl_nu; idx++) {
    // contribution of objective term
    if (_mdl_der_u.active[idx] || k == _K && _mdl_u.active[idx]) {
      // controlled inputs
      dtu = (_mdl_u_order[idx] == 0)? dt0: dt;
      help = (x[i] - _mdl_u.ref[idx]/_mdl_u_nominal[idx]);
      f0x[i] += dtu * _mdl_u.weight1[idx];
      f0x[i] += dtu * _mdl_u.weight2[idx] * 2.0 * help;
      if (_t_active) {
        if (_mdl_u_order[idx] == 0) {
          f0x[t_scale_ix] += ddt0 * _mdl_u.weight1[idx] * help;
          f0x[t_scale_ix] += ddt0 * _mdl_u.weight2[idx] * help*help;
          if (k < _K) {
            f0u[t_scale_iu] += ddt0 * _mdl_u.weight1[idx] * help;
            f0u[t_scale_iu] += ddt0 * _mdl_u.weight2[idx] * help*help;
          }
        }
        else if (_mdl_u_order[idx] > 0) {
          f0x[t_scale_ix] += ddtx * _mdl_u.weight1[idx] * help;
          f0x[t_scale_ix] += ddtx * _mdl_u.weight2[idx] * help*help;
          if (k < _K) {
            f0u[t_scale_iu] += ddtu * _mdl_u.weight1[idx] * help;
            f0u[t_scale_iu] += ddtu * _mdl_u.weight2[idx] * help*help;
          }
        }
      }
      // rates of change
      if (k < _K) {
	help = u[(i-_mdl_nd)*upsk + k%upsk]/_t_nominal - _mdl_der_u.ref[idx]
	  / _mdl_u_nominal[idx];
	f0u[(i-_mdl_nd)*upsk + k%upsk] += dt0 * _mdl_der_u.weight1[idx]/_t_nominal;
	f0u[(i-_mdl_nd)*upsk + k%upsk] += dt0 * _mdl_der_u.weight2[idx]/_t_nominal * 2.0 * help;
        if (_t_active) {
          f0x[t_scale_ix] += ddt0 * _mdl_der_u.weight1[idx] * help;
          f0x[t_scale_ix] += ddt0 * _mdl_der_u.weight2[idx] * help*help;
          f0u[t_scale_iu] += ddt0 * _mdl_der_u.weight1[idx] * help;
          f0u[t_scale_iu] += ddt0 * _mdl_der_u.weight2[idx] * help*help;
        }
        // soft constraints on rates of change
        if (_mdl_der_u_soft.active[idx]) {
          help = u[upsk*_nu + is + k%upsk]/_t_nominal;
          f0u[upsk*_nu + is + k%upsk] += dt0 * _mdl_der_u_soft.weight1[idx]/_t_nominal;
          f0u[upsk*_nu + is + k%upsk] += dt0 * _mdl_der_u_soft.weight2[idx]/_t_nominal * 2.0 * help;
          if (_t_active) {
            f0x[t_scale_ix] += ddt0 * _mdl_der_u_soft.weight1[idx] * help;
            f0x[t_scale_ix] += ddt0 * _mdl_der_u_soft.weight2[idx] * help*help;
            f0u[t_scale_iu] += ddt0 * _mdl_der_u_soft.weight1[idx] * help;
            f0u[t_scale_iu] += ddt0 * _mdl_der_u_soft.weight2[idx] * help*help;
          }
          if (analyticJac && _mdl_der_u_soft.min[idx] > -Inf) {
            cu[spsk*(_nc+_nsc) + isc + k%upsk][i*upsk + k%upsk]
              += 1.0/_t_nominal;
            cu[spsk*(_nc+_nsc) + isc + k%upsk][upsk*_nu + is + k%upsk]
              += 1.0/_t_nominal;
            isc += upsk;
          }
          if (analyticJac && _mdl_der_u_soft.max[idx] < Inf) {
            cu[spsk*(_nc+_nsc) + isc + k%upsk][i*upsk + k%upsk]
              += 1.0/_t_nominal;
            cu[spsk*(_nc+_nsc) + isc + k%upsk][upsk*_nu + is + k%upsk]
              -= 1.0/_t_nominal;
            isc += upsk;
          }
          is += upsk;
        }
      }
      i++;
    }
    else if (_mdl_u.active[idx]) {
      // controlled inputs
      dtu = (_mdl_u_order[idx] == 0)? dt0: dt;
      help = (u[iu] - _mdl_u.ref[idx]/_mdl_u_nominal[idx]);
      f0u[iu] += dtu * _mdl_u.weight1[idx];
      f0u[iu] += dtu * _mdl_u.weight2[idx] * 2.0 * help;
      iu++;
    }
  }
  // active outputs
  for (i = 0, is = 0, idx = 0; idx < _mdl_ny; idx++) {
    dty = (_mdl_y_order[idx] == 0)? dt0: dt;
    // contribution of objective term
    if (_mdl_y.active[idx] == 1) {
      help = (c[i + k%spsk] - _mdl_y.ref[idx]/_mdl_y_nominal[idx]);
      if (_mdl_y.weight1[idx] != 0.0) {
        for (j = 0; j < nxk; j++)
          f0x[j] += dty * _mdl_y.weight1[idx] * cx[i + k%spsk][j];
        if (k < _K && _nu > _ndu)
          for (j = 0; j < _nu; j++)
            f0u[j] += dty * _mdl_y.weight1[idx] * cu[i + k%spsk][j];
        if (_t_active) {
          f0x[t_scale_ix] += ddtx * _mdl_y.weight1[idx] * help;
          if (k < _K)
            f0u[t_scale_iu] += ddtu * _mdl_y.weight1[idx] * help;
        }
      }
      if (_mdl_y.weight2[idx] != 0.0) {
        for (j = 0; j < nxk; j++)
          f0x[j] += dty * _mdl_y.weight2[idx] *
            2.0 * help * cx[i + k%spsk][j];
        if (k < _K && _nu > _ndu)
          for (j = 0; j < _nu; j++)
            f0u[j] += dty * _mdl_y.weight2[idx] *
              2.0 * help * cu[i + k%spsk][j];
        if (_t_active) {
          f0x[t_scale_ix] += ddtx * _mdl_y.weight2[idx] * help*help;
          if (k < _K)
            f0u[t_scale_iu] += ddtu * _mdl_y.weight2[idx] * help*help;
        }
      }
      i += spsk;
    }

    // contribution of soft constraints
    if (_mdl_y_soft.active[idx] == 1) {
      if (k < _K) {
        help = u[upsk*(_nu+_nsu) + is + k%spsk];
	f0u[upsk*(_nu+_nsu) + is + k%spsk]
	  += dty * (_mdl_y_soft.weight1[idx]
		    + 2.0*_mdl_y_soft.weight2[idx]*help);
      }
      else {
        help = x[nxk + is + k%spsk];
	f0x[nxk + is + k%spsk]
	  += dty * (_mdl_y_soft.weight1[idx]
		    + 2.0*_mdl_y_soft.weight2[idx]*help);
        
      }
      if (_t_active) {
        f0x[t_scale_ix] += ddtx * (_mdl_y_soft.weight1[idx]*help
                                     + _mdl_y_soft.weight2[idx]*help*help);
        if (k < _K)
          f0u[t_scale_iu] += ddtu * (_mdl_y_soft.weight1[idx]*help
                                       + _mdl_y_soft.weight2[idx]*help*help);
      }
      is += spsk;
    }
  }

  // additional terms at initial time
  if (k == 0) {
    for (i = spsk*(_nc+_nsc) + upsk*_nsuc, idx = 0; idx < _mdl_ny; idx++) {
      // initial objective terms
      if (_mdl_y0.active[idx]) {
        if (_mdl_y0.weight1[idx] != 0.0) {
          for (j = 0; j < nxk; j++)
            f0x[j] += _mdl_y0.weight1[idx] * cx[i][j];
          if (k < _K && _nu > _ndu)
            for (j = 0; j < _nu; j++)
              f0u[j] += _mdl_y0.weight1[idx] * cu[i][j];
        }
        if (_mdl_y0.weight2[idx] != 0.0) {
          for (j = 0; j < nxk; j++)
            f0x[j] += _mdl_y0.weight2[idx] *
              2.0 * (c[i] - _mdl_y.ref[idx]/_mdl_y_nominal[idx]) *
              cx[i][j];
          if (k < _K && _nu > _ndu)
            for (j = 0; j < _nu; j++)
              f0u[j] += _mdl_y0.weight2[idx] *
                2.0 * (c[i] - _mdl_y.ref[idx]/_mdl_y_nominal[idx]) *
                cu[i][j];
        }
	i++;
      }
    }
  }
  // additional terms at final time
  if (k == _K) {
    for (i = _nc+_nsc, isf = 0, idx = 0; idx < _mdl_ny; idx++) {
      // contributions of final objective terms
      if (_mdl_yf.active[idx]) {
        if (_mdl_yf.weight1[idx] != 0.0) {
          for (j = 0; j < nxk; j++)
            f0x[j] += _mdl_yf.weight1[idx] * cx[i][j];
        }
        if (_mdl_yf.weight2[idx] != 0.0) {
          for (j = 0; j < nxk; j++)
            f0x[j] += _mdl_yf.weight2[idx] *
              2.0 * (c[i] - _mdl_y.ref[idx]/_mdl_y_nominal[idx]) *
              cx[i][j];
        }
	i++;
      }
      // contribution of soft constraints
      if (_mdl_yf_soft.active[idx] == 1) {
        f0x[nxk + spsk*_ns + isf] += 
          (_mdl_yf_soft.weight1[idx]
           + 2.0*_mdl_yf_soft.weight2[idx]*x[nxk + spsk*_ns + isf]);
        isf++;
      }
    }
  }
  sv_mlt(_fscale, f0x, f0x);
  sv_mlt(_fscale, f0u, f0u);
}

//--------------------------------------------------------------------------
void Prg_DTOpt::fetch_jac(SimStruct *S,
                          int k, double tscale, const VECP x, const VECP u,
                          MATP fx, MATP fu, MATP cx, MATP cu)
{
  int i, iu, idx, ii, iidx0, iidx, isc, iscdx, is, j, ju, jdx, rdx;
  int isf, iscf, iscfdx, ixf;
  int nxk = k < _K? _nx: _nu - _ndu + _nx;
  int spsk = 1;
  int upsk = 1;
  int t_scale_iu = _t_scale_i*upsk + k%upsk;
  int t_scale_ix = _mdl_nd + _t_scale_i;
  double help;
  int mdl_u_idx, mdl_x_idx, mdl_y_idx;
  int mdl_nc = _mdl_nx - _mdl_nd; // number of continuous states
  real_T *pr = ssGetJacobianPr(S);
  int_T *ir = ssGetJacobianIr(S);
  int_T *jc = ssGetJacobianJc(S);

  // default values
  m_zero(fx);
  m_zero(fu);
  m_zero(cx);
  m_zero(cu);
  // Jacobian wrt S-function states (ddxdy/dxcxd)
  for (jdx = 0; jdx < _mdl_nx; jdx++) {
    j = jdx < mdl_nc? _mdl_nd + _nu + jdx: jdx - mdl_nc;
    mdl_x_idx = jdx < mdl_nc? _mdl_nd + jdx: jdx - mdl_nc;
    for (i = 0, idx = _mdl_nx, isc = spsk*_nc, iscdx = _mdl_nx,
           iscf = spsk*(_nc+_nsc) + _ncf, iscfdx = _mdl_nx,
           ii = _nc+_nsc, iidx0 = 0, iidx = _mdl_nx,
           rdx = jc[jdx]; rdx < jc[jdx+1]; rdx++) {
      if (ir[rdx] < _mdl_nx && k < _K) {
        // junction conditions for discrete states
        ixf = ir[rdx] - mdl_nc;
        fx[ixf][j] = pr[rdx] /
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
          cx[i + k%spsk][j] = pr[rdx] /
            _mdl_y_nominal[mdl_y_idx] * _mdl_x_nominal[mdl_x_idx];
        }
        // soft constraints
        if (_mdl_y_soft.active[mdl_y_idx] == 1) {
          // need to loop to obtain isc considering active outputs
          for (; iscdx < ir[rdx]; iscdx++) {
            if (_mdl_y_soft.active[iscdx - _mdl_nx] == 1) {
              if (_mdl_y_soft.min[iscdx - _mdl_nx] > -Inf)
                isc += spsk;
              if (_mdl_y_soft.max[iscdx - _mdl_nx] < Inf)
                isc += spsk;
            }
          }
          if (_mdl_y_soft.min[mdl_y_idx] > -Inf) {
            cx[isc + k%spsk][j] = pr[rdx] /
              _mdl_y_nominal[mdl_y_idx] * _mdl_x_nominal[mdl_x_idx];
            isc += spsk;
          }
          if (_mdl_y_soft.max[mdl_y_idx] < Inf) {
            cx[isc + k%spsk][j] = pr[rdx] /
              _mdl_y_nominal[mdl_y_idx] * _mdl_x_nominal[mdl_x_idx];
          }
          if (_mdl_y_soft.min[mdl_y_idx] > -Inf) {
            isc -= spsk;
	  }
        }
        // constraints (used model outputs) at initial time
        if (k == 0 && _mdl_y0.active[mdl_y_idx]) {
          // need to loop to obtain ii considering active outputs
          for (; iidx < ir[rdx]; iidx++) {
            if (_mdl_y0.active[iidx - _mdl_nx])
              ii++;
          }
          cx[ii][j] = pr[rdx] /
            _mdl_y_nominal[mdl_y_idx] * _mdl_x_nominal[mdl_x_idx];
        }
        // constraints (used model outputs) at final time
        if (k == _K && _mdl_yf.active[mdl_y_idx]) {
          // need to loop to obtain ii considering active outputs
          for (; iidx < ir[rdx]; iidx++) {
            if (_mdl_yf.active[iidx - _mdl_nx])
              ii++;
          }
          cx[ii][j] = pr[rdx] /
            _mdl_y_nominal[mdl_y_idx] * _mdl_x_nominal[mdl_x_idx];
        }
        // soft constraints at final time
        if (k == _K && _mdl_yf_soft.active[mdl_y_idx] == 1) {
          // need to loop to obtain iscf considering active outputs
          for (; iscfdx < ir[rdx]; iscfdx++) {
            if (_mdl_yf_soft.active[iscfdx - _mdl_nx] == 1) {
              if (_mdl_yf_soft.min[iscfdx - _mdl_nx] > -Inf)
                iscf++;
              if (_mdl_yf_soft.max[iscfdx - _mdl_nx] < Inf)
                iscf++;
            }
          }
          if (_mdl_yf_soft.min[mdl_y_idx] > -Inf) {
            cx[iscf][j] = pr[rdx] /
              _mdl_y_nominal[mdl_y_idx] * _mdl_x_nominal[mdl_x_idx];
            iscf++;
          }
          if (_mdl_yf_soft.max[mdl_y_idx] < Inf) {
            cx[iscf][j] = pr[rdx] /
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
  for (j = _mdl_nd, ju = 0, jdx = _mdl_nx; jdx < _mdl_nx + _mdl_nu; jdx++) {
    mdl_u_idx = jdx - _mdl_nx;
    if (_mdl_der_u.active[mdl_u_idx] || k == _K && _mdl_u.active[mdl_u_idx]) {
      for (i = 0, idx = _mdl_nx, isc = spsk*_nc, iscdx = _mdl_nx,
             iscf = spsk*(_nc+_nsc) + _ncf, iscfdx = _mdl_nx,
             ii = _nc+_nsc, iidx0 = 0, iidx = _mdl_nx,
             rdx = jc[jdx]; rdx < jc[jdx+1]; rdx++) {
        if (ir[rdx] < _mdl_nx && k < _K) {
          // junction conditions for discrete states
          ixf = ir[rdx] - mdl_nc;
          fx[ixf][j] = pr[rdx] /
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
            cx[i + k%spsk][j] = pr[rdx] /
              _mdl_y_nominal[mdl_y_idx] * _mdl_u_nominal[mdl_u_idx];
          }
          // soft constraints
          if (_mdl_y_soft.active[mdl_y_idx] == 1) {
            // need to loop to obtain isc considering active outputs
            for (; iscdx < ir[rdx]; iscdx++) {
              if (_mdl_y_soft.active[iscdx - _mdl_nx] == 1) {
                if (_mdl_y_soft.min[iscdx - _mdl_nx] > -Inf)
                  isc += spsk;
                if (_mdl_y_soft.max[iscdx - _mdl_nx] < Inf)
                  isc += spsk;
              }
            }
            if (_mdl_y_soft.min[mdl_y_idx] > -Inf) {
              cx[isc + k%spsk][j] = pr[rdx] /
                _mdl_y_nominal[mdl_y_idx] * _mdl_u_nominal[mdl_u_idx];
              isc += spsk;
            }
            if (_mdl_y_soft.max[mdl_y_idx] < Inf) {
              cx[isc + k%spsk][j] = pr[rdx] /
                _mdl_y_nominal[mdl_y_idx] * _mdl_u_nominal[mdl_u_idx];
            }
            if (_mdl_y_soft.min[mdl_y_idx] > -Inf) {
              isc -= spsk;
            }
          }
          // constraints at initial time
          if (k == 0 && _mdl_y0.active[mdl_y_idx]) {
            // need to loop to obtain ii considering active outputs
            for (; iidx < ir[rdx]; iidx++) {
              if (_mdl_y0.active[iidx - _mdl_nx])
                ii++;
            }
            cx[ii][j] = pr[rdx] /
              _mdl_y_nominal[mdl_y_idx] * _mdl_u_nominal[mdl_u_idx];
          }
          // constraints at final time
          if (k == _K && _mdl_yf.active[mdl_y_idx]) {
            // need to loop to obtain ii considering active outputs
            for (; iidx < ir[rdx]; iidx++) {
              if (_mdl_yf.active[iidx - _mdl_nx])
                ii++;
            }
            cx[ii][j] = pr[rdx] /
              _mdl_y_nominal[mdl_y_idx] * _mdl_u_nominal[mdl_u_idx];
          }
          // soft constraints at final time
          if (k == _K && _mdl_yf_soft.active[mdl_y_idx] == 1) {
            // need to loop to obtain iscf considering active outputs
            for (; iscfdx < ir[rdx]; iscfdx++) {
              if (_mdl_yf_soft.active[iscfdx - _mdl_nx] == 1) {
                if (_mdl_yf_soft.min[iscfdx - _mdl_nx] > -Inf)
                  iscf++;
                if (_mdl_yf_soft.max[iscfdx - _mdl_nx] < Inf)
                  iscf++;
              }
            }
            if (_mdl_yf_soft.min[mdl_y_idx] > -Inf) {
              cx[iscf][j] = pr[rdx] /
                _mdl_y_nominal[mdl_y_idx] * _mdl_u_nominal[mdl_u_idx];
              iscf++;
            }
            if (_mdl_yf_soft.max[mdl_y_idx] < Inf) {
              cx[iscf][j] = pr[rdx] /
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
    else if (_mdl_u.active[mdl_u_idx]) {
      for (i = 0, idx = _mdl_nx, isc = spsk*_nc, iscdx = _mdl_nx,
             iscf = spsk*(_nc+_nsc) + _ncf, iscfdx = _mdl_nx,
             ii = _nc+_nsc, iidx0 = 0, iidx = _mdl_nx,
             rdx = jc[jdx]; rdx < jc[jdx+1]; rdx++) {
        if (ir[rdx] < _mdl_nx && k < _K) {
          // junction conditions for discrete states
          ixf = ir[rdx] - mdl_nc;
          fu[ixf][ju] = pr[rdx] /
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
            cu[i + k%spsk][ju] = pr[rdx] /
              _mdl_y_nominal[mdl_y_idx] * _mdl_u_nominal[mdl_u_idx];
          }
          // soft constraints
          if (_mdl_y_soft.active[mdl_y_idx] == 1) {
            // need to loop to obtain isc considering active outputs
            for (; iscdx < ir[rdx]; iscdx++) {
              if (_mdl_y_soft.active[iscdx - _mdl_nx] == 1) {
                if (_mdl_y_soft.min[iscdx - _mdl_nx] > -Inf)
                  isc += spsk;
                if (_mdl_y_soft.max[iscdx - _mdl_nx] < Inf)
                  isc += spsk;
              }
            }
            if (_mdl_y_soft.min[mdl_y_idx] > -Inf) {
              cu[isc + k%spsk][ju] = pr[rdx] /
                _mdl_y_nominal[mdl_y_idx] * _mdl_u_nominal[mdl_u_idx];
              isc += spsk;
            }
            if (_mdl_y_soft.max[mdl_y_idx] < Inf) {
              cu[isc + k%spsk][ju] = pr[rdx] /
                _mdl_y_nominal[mdl_y_idx] * _mdl_u_nominal[mdl_u_idx];
            }
            if (_mdl_y_soft.min[mdl_y_idx] > -Inf) {
              isc -= spsk;
            }
          }
          // constraints at initial time
          if (k == 0 && _mdl_y0.active[mdl_y_idx]) {
            // need to loop to obtain ii considering active outputs
            for (; iidx < ir[rdx]; iidx++) {
              if (_mdl_y0.active[iidx - _mdl_nx])
                ii++;
            }
            cu[ii][ju] = pr[rdx] /
              _mdl_y_nominal[mdl_y_idx] * _mdl_u_nominal[mdl_u_idx];
          }
          // constraints at final time
          if (k == _K && _mdl_yf.active[mdl_y_idx]) {
            // need to loop to obtain ii considering active outputs
            for (; iidx < ir[rdx]; iidx++) {
              if (_mdl_yf.active[iidx - _mdl_nx])
                ii++;
            }
            cu[ii][ju] = pr[rdx] /
              _mdl_y_nominal[mdl_y_idx] * _mdl_u_nominal[mdl_u_idx];
          }
          // soft constraints at final time
          if (k == _K && _mdl_yf_soft.active[mdl_y_idx] == 1) {
            // need to loop to obtain iscf considering active outputs
            for (; iscfdx < ir[rdx]; iscfdx++) {
              if (_mdl_yf_soft.active[iscfdx - _mdl_nx] == 1) {
                if (_mdl_yf_soft.min[iscfdx - _mdl_nx] > -Inf)
                  iscf++;
                if (_mdl_yf_soft.max[iscfdx - _mdl_nx] < Inf)
                  iscf++;
              }
            }
            if (_mdl_yf_soft.min[mdl_y_idx] > -Inf) {
              cu[iscf][ju] = pr[rdx] /
                _mdl_y_nominal[mdl_y_idx] * _mdl_u_nominal[mdl_u_idx];
              iscf++;
            }
            if (_mdl_yf_soft.max[mdl_y_idx] < Inf) {
              cu[iscf][ju] = pr[rdx] /
                _mdl_y_nominal[mdl_y_idx] * _mdl_u_nominal[mdl_u_idx];
            }
            if (_mdl_yf_soft.min[mdl_y_idx] > -Inf) {
              iscf--;
            }
          }
        }
      }
      ju++;
    }
  }
  // contribution of slack variables for soft constraints
  for (is = 0, isc = spsk*_nc, isf = 0, iscf = 0,
         idx = 0; idx < _mdl_ny; idx++) {
    if (_mdl_y_soft.active[idx] == 1) {
      if (_mdl_y_soft.min[idx] > -Inf) {
        if (k < _K)
          cu[isc + k%spsk][upsk*(_nu+_nsu) + is + k%spsk] += 1.0;
        else
          cx[isc + k%spsk][nxk + is + k%spsk] += 1.0;
        isc += spsk;
      }
      if (_mdl_y_soft.max[idx] < Inf) {
        if (k < _K)
          cu[isc + k%spsk][upsk*(_nu+_nsu) + is + k%spsk] -= 1.0;
        else
          cx[isc + k%spsk][nxk + is + k%spsk] -= 1.0;
        isc += spsk;
      }
      is += spsk;
    }
    // soft constraints at final time
    if (k == _K && _mdl_yf_soft.active[idx] == 1) {
      if (_mdl_yf_soft.min[idx] > -Inf) {
        cx[spsk*(_nc+_nsc) + _ncf + iscf][nxk + spsk*_ns + isf] += 1.0;
        iscf++;
      }
      if (_mdl_yf_soft.max[idx] < Inf) {
        cx[spsk*(_nc+_nsc) + _ncf + iscf][nxk + spsk*_ns + isf] -= 1.0;
        iscf++;
      }
      isf++;
    }
  }

  // junction conditions for continuous-time equations
  // Note: fx for discrete states is filled above
  if (k < _K) {
    // modifications for controlled inputs with zero order hold
    for (i = _mdl_nd, idx = 0; idx < _mdl_nu; idx++) {
      if (_mdl_der_u.active[idx]) {
        if (true || _mdl_u_order[idx] == 0 && (k+1)%spsk == 0) {
          // zero order hold at end of stage:
          // apply step in u for subsequent stage
          fx[i][i] = 1.0;
          help = (ts(k+1) - ts(k-spsk/upsk+1))/_t_nominal;
          if (_t_active && i-_mdl_nd != _t_scale_i) {
            fu[i][(i-_mdl_nd)*upsk + k%upsk] = tscale*help;
            fx[i][t_scale_ix] = _t_scale_nominal*help*u[(i-_mdl_nd)*upsk + k%upsk];
            fu[i][t_scale_iu] = fx[i][t_scale_ix];
          }
          else
            fu[i][(i-_mdl_nd)*upsk + k%upsk] = help;
        }
        //else {
        //  fxf[i][i] = 1.0;
        //}
        i++;
      }
      else if (_mdl_u.active[idx] && k == _K - 1) {
        // propagate control inputs to final time
        fu[i][i-_mdl_nd] = 1.0;
        i++;
      }
    }
    // problem structure must not have changed
    assert(i == _mdl_nd + _ndu || k == _K - 1 && i == _mdl_nd + _nu);
    // slack variables at final time
    if (k == _K-1) {
      int nxK = _nx + _nu - _ndu; // nxk of final stage
      for (; i < nxK + _ns; i++)
        fu[i][upsk*(_nu+_nsu) + (i-nxK+1)*spsk-1] = 1.0;
      for (; i < nxK + _ns + _nsf; i++)
        fu[i][upsk*(_nu+_nsu) + spsk*_ns + i - nxK - _ns] = 1.0;
    }
  }
}

//==========================================================================

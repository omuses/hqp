/*
 * Prg_SFunctionEst.C -- class definition
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
#define assert(expr) if (!(expr)) error(E_INTERN, "assert(" #expr ")");

IF_CLASS_DEFINE("SFunctionEst", Prg_SFunctionEst, Omu_Program);

// define own Inf so that iftcl works
static const double INF = 8e88;

//#define PRG_WITH_U   1

//--------------------------------------------------------------------------
Prg_SFunctionEst::Prg_SFunctionEst()
{
  _K = 1;
  _KK = 1;
  _multistage = true;

  _mx_p = NULL;
  _mdl_args_p_idx = -1;

  _mdl_np = 0;

  _mdl_p_active = iv_get(_mdl_np);
  _mdl_x0_active = iv_get(_mdl_nx);
  _mdl_x0_steady = iv_get(_mdl_nx);
  _mdl_y_active = iv_get(_mdl_ny);
  iv_zero(_mdl_p_active);
  iv_zero(_mdl_x0_active);
  iv_zero(_mdl_x0_steady);
  iv_zero(_mdl_y_active);
  _np = 0;
  _nx0 = 0;
  _ny = 0;
  _nx = _mdl_nx;

  _mdl_us = m_get(_KK, _mdl_nu);
  _mdl_ys = m_get(_KK, _mdl_ny);
  _prg_ys_ref = m_get(_KK, 1);
  _prg_y_nominal = v_get(1);
  v_ones(_prg_y_nominal);

  _mdl_p = v_get(_mdl_np);
  _mdl_p_nominal = v_get(_mdl_np);
  _mdl_p_min = v_get(_mdl_np);
  _mdl_p_max = v_get(_mdl_np);
  v_set(_mdl_p_nominal, 1.0);
  v_set(_mdl_p_min, 0.0);
  v_set(_mdl_p_max, INF);

  _mdl_x0 = v_get(_mdl_nx);
  _mdl_x0_min = v_get(_mdl_nx);
  _mdl_x0_nominal = v_get(_mdl_nx);
  _mdl_x0_max = v_get(_mdl_nx);
  v_set(_mdl_x0_nominal, 1.0);
  v_set(_mdl_x0_min, -INF);
  v_set(_mdl_x0_max, INF);

  _M = m_resize(m_get(1, 1), _K*_ny, _np+_nx0);
  _dxdpx0 = m_resize(m_get(1, 1), _nx, _np+_nx0);
  _dxfdpx0 = m_resize(m_get(1, 1), _nx, _np+_nx0);

  _ifList.append(new If_Int("mdl_args_p_idx", &_mdl_args_p_idx));

  _ifList.append(new If_IntVec("mdl_p_active", &_mdl_p_active));
  _ifList.append(new If_IntVec("mdl_x0_active", &_mdl_x0_active));
  _ifList.append(new If_IntVec("mdl_x0_steady", &_mdl_x0_steady));
  _ifList.append(new If_IntVec("mdl_y_active", &_mdl_y_active));

  _ifList.append(new If_RealVec("mdl_p", &_mdl_p));
  _ifList.append(new If_RealVec("mdl_p_nominal", &_mdl_p_nominal));
  _ifList.append(new If_RealVec("mdl_p_min", &_mdl_p_min));
  _ifList.append(new If_RealVec("mdl_p_max", &_mdl_p_max));
  _ifList.append(new If_RealVec("mdl_x0", &_mdl_x0));
  _ifList.append(new If_RealVec("mdl_x0_nominal", &_mdl_x0_nominal));
  _ifList.append(new If_RealVec("mdl_x0_min", &_mdl_x0_min));
  _ifList.append(new If_RealVec("mdl_x0_max", &_mdl_x0_max));

  _ifList.append(new If_RealMat("mdl_us", &_mdl_us));
  _ifList.append(new If_RealMat("mdl_ys", &_mdl_ys));
  _ifList.append(new If_Bool("prg_multistage", &_multistage));

  _ifList.append(new If_RealMat("prg_ys_ref", &_prg_ys_ref));
  _ifList.append(new If_RealVec("prg_y_nominal", &_prg_y_nominal));

  _ifList.append(new If_RealMat("prg_M", &_M));
}

//--------------------------------------------------------------------------
Prg_SFunctionEst::~Prg_SFunctionEst()
{
  m_free(_dxfdpx0);
  m_free(_dxdpx0);
  m_free(_M);
  iv_free(_mdl_p_active);
  iv_free(_mdl_x0_steady);
  iv_free(_mdl_x0_active);
  iv_free(_mdl_y_active);
  v_free(_prg_y_nominal);
  m_free(_prg_ys_ref);
  m_free(_mdl_ys);
  m_free(_mdl_us);
  v_free(_mdl_p_min);
  v_free(_mdl_p_max);
  v_free(_mdl_p_nominal);
  v_free(_mdl_p);
  v_free(_mdl_x0_min);
  v_free(_mdl_x0_max);
  v_free(_mdl_x0_nominal);
  v_free(_mdl_x0);
}

//--------------------------------------------------------------------------
void Prg_SFunctionEst::setup_stages(IVECP ks, VECP ts)
{
  int i;

  // setup optimization problem
  if (_multistage) {
    _K = _KK;
    stages_alloc(ks, ts, _K, 1);
  }
  else {
    _K = 2;
    v_resize(ts, _KK + 1);
    iv_resize(ks, 3);
    ks[0] = 0;
    ks[1] = _KK - 1;
    ks[2] = _KK;
  }

  // setup S-function
  setup_sfun();

  // check for optional S-function methods that are required
  assert(ssGetmdlDerivatives(_S) != NULL);

  if (_mdl_args_p_idx >= 0) {
    _mx_p = mxGetCell(_mx_args, _mdl_args_p_idx);
    _mdl_np = mxGetNumberOfElements(_mx_p);
  }
  else {
    _mx_p = NULL;
    _mdl_np = 0;
  }
  
  // adapt sizes of model vectors
  iv_resize(_mdl_p_active, _mdl_np);
  iv_resize(_mdl_x0_active, _mdl_nx);
  iv_resize(_mdl_x0_steady, _mdl_nx);
  iv_resize(_mdl_y_active, _mdl_ny);
  iv_zero(_mdl_p_active);
  iv_zero(_mdl_x0_active);
  iv_zero(_mdl_x0_steady);
  iv_zero(_mdl_y_active);

  v_resize(_mdl_p, _mdl_np);
  v_resize(_mdl_p_nominal, _mdl_np);
  v_resize(_mdl_p_min, _mdl_np);
  v_resize(_mdl_p_max, _mdl_np);
  v_set(_mdl_p_nominal, 1.0);
  v_set(_mdl_p_min, 0.0);
  v_set(_mdl_p_max, INF);

  v_resize(_mdl_x0, _mdl_nx);
  v_resize(_mdl_x0_nominal, _mdl_nx);
  v_resize(_mdl_x0_min, _mdl_nx);
  v_resize(_mdl_x0_max, _mdl_nx);
  v_set(_mdl_x0_nominal, 1.0);
  v_set(_mdl_x0_min, -INF);
  v_set(_mdl_x0_max, INF);

  m_resize(_mdl_us, _KK, _mdl_nu);
  m_resize(_mdl_ys, _KK, _mdl_ny);

  // store parameters in _mdl_p
  for (i = 0; i < _mdl_np; i++)
    _mdl_p[i] = mxGetPr(_mx_p)[i];
}

//--------------------------------------------------------------------------
void Prg_SFunctionEst::setup(int k,
			     Omu_Vector &x, Omu_Vector &u, Omu_Vector &c)
{
  int i;

  // take over possibly modified _mdl_p
  for (i = 0; i < _mdl_np; i++)
    mxGetPr(_mx_p)[i] = _mdl_p[i];

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
    m_resize(_prg_ys_ref, _KK, _ny);
    v_resize(_prg_y_nominal, _ny);
    v_ones(_prg_y_nominal);

    m_resize(_M, _KK*_ny, _np+_nx0);
    m_resize(_dxdpx0, _nx, _np+_nx0);
    m_resize(_dxfdpx0, _nx, _np+_nx0);
  }

  if (k < _K)
    x.alloc(_nx);
  else
    x.alloc(_np);

  if (k == 0) {
    // setup parameters
    int ip = 0;
    for (i = 0; i < _mdl_np; i++) {
      if (_mdl_p_active[i]) {
	x.initial[ip] = _mdl_p[i] / _mdl_p_nominal[i];
	if (_mdl_p_min[i] > -INF)
	  x.min[ip] = _mdl_p_min[i] / _mdl_p_nominal[i];
	if (_mdl_p_max[i] < INF)
	  x.max[ip] = _mdl_p_max[i] / _mdl_p_nominal[i];
	++ip;
      }
    }
    assert(ip == _np);	// _np must not have changed since setup_stages

    // setup initial states
    for (i = _np; i < _nx; i++) {
      x.initial[i] = _mdl_x0[i-_np] / _mdl_x0_nominal[i-_np];
      if (!_mdl_x0_active[i-_np])
	x.min[i] = x.max[i] = x.initial[i];
      else {
	if (_mdl_x0_min[i-_np] > -INF)
	  x.min[i] = _mdl_x0_min[i-_np] / _mdl_x0_nominal[i-_np];
	if (_mdl_x0_max[i-_np] < INF)
	  x.max[i] = _mdl_x0_max[i-_np] / _mdl_x0_nominal[i-_np];
      }
    }
  }
  // constraints for assigning relevant model outputs to opt vars
  // and for steady initial states
  if (k < _K) {
    int nc = _ny;
    if (k == 0) {
      for (i = 0; i < _mdl_nx; i++) {
	if (_mdl_x0_steady[i])
	  ++nc;
      }
    }
    c.alloc(nc);
#if PRG_WITH_U
    u.alloc(_ny);
    for (i = 0; i < _ny; i++) {
      c.min[i] = c.max[i] = 0.0;
    }
#endif
    for (i = _ny; i < nc; i++) {
      c.min[i] = c.max[i] = 0.0;
    }
  }
}

//--------------------------------------------------------------------------
void Prg_SFunctionEst::setup_struct(int k,
				    const Omu_Vector &x, const Omu_Vector &u,
				    Omu_DependentVec &xt, Omu_DependentVec &F,
				    Omu_DependentVec &f,
				    Omu_Dependent &f0, Omu_DependentVec &c)
{
  int i, j;

  // consistic just takes states from optimizer
  m_ident(xt.Jx);
  m_zero(xt.Ju);
  xt.set_linear();

  if (k < _K-1) {
    // explicit ODE for continuous-time equations
    for (i = 0; i < _np; i++) {
      for (j = 0; j < _nx; j++)
	F.Jx[i][j] = 0.0;
#if PRG_WITH_U
      for (j = 0; j < _ny; j++)
	F.Ju[i][j] = 0.0;
#endif
    }
    m_zero(F.Jxp);
    for (i = _np; i < _nx; i++)
      F.Jxp[i][i] = -1.0;
    F.set_linear(Omu_Dependent::WRT_xp);
  }
  else {
    // no continuous equations in last stage
    m_zero(F.Jx);
    m_zero(F.Ju);
    m_zero(F.Jxp);
    F.set_linear();
  }
}

//--------------------------------------------------------------------------
void Prg_SFunctionEst::init_simulation(int k,
				       Omu_Vector &x, Omu_Vector &u)
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
  int i, j;
  //double dt = kk < _KK? ts(kk+1)-ts(kk): 0.0;

  if (kk < _KK) {
    // set simulation time
    ssSetT(_S, ts(kk));

    // pass estimated parameters to model
    int ip = 0;
    for (i = 0; i < _mdl_np; i++) {
      if (_mdl_p_active[i]) {
	mxGetPr(_mx_p)[i] = x[ip++] * _mdl_p_nominal[i];
      }
    }
    assert(ip == _np);	// _np must not have changed since setup_stages

    // initialize model inputs
    real_T **mdl_u = (real_T **)ssGetInputPortRealSignalPtrs(_S, 0);
    for (i = 0; i < _mdl_nu; i++)
      *mdl_u[i] = _mdl_us[kk][i];

    // initialize model (this is required as parameters change
    // and as time steps back, compared to previous continous call)
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
      mdl_x[i] = x[_np + i] * _mdl_x0_nominal[i];

    // obtain model outputs
    mdlOutputs(_S, 0);
    if (ssGetErrorStatus(_S)) {
      fprintf(stderr, "Error from mdlOutputs: %s\n", ssGetErrorStatus(_S));
      ssSetErrorStatus(_S, NULL);
      error(E_RANGE, "mdlOutputs");
    }

    // store outputs in constraints
    real_T *mdl_y = ssGetOutputPortRealSignal(_S, 0);
    int iy = 0;
    for (i = 0; i < _mdl_ny; i++) {
      if (_mdl_y_active[i]) {
#if PRG_WITH_U
	c[iy] = mdl_y[i] - u[iy];
#else
	c[iy] = mdl_y[i];
#endif
	++iy;
      }
    }

    // calculate objective term
    f0 = 0.0;
    double help;
    for (i = 0; i < _ny; i++) {
#if PRG_WITH_U
      help = (u[i] - _prg_ys_ref[kk][i]) / _prg_y_nominal[i];
#else
      help = (c[i] - _prg_ys_ref[kk][i]) / _prg_y_nominal[i];
#endif
      f0 += help*help;
    }
    //f0 *= dt;

    // steady initial state constraints
    if (kk == 0) {
      real_T *mdl_xp = ssGetdX(_S);
      mdlDerivatives(_S);
      if (ssGetErrorStatus(_S)) {
	fprintf(stderr, "Error from mdlDerivatives: %s\n",
		ssGetErrorStatus(_S));
	ssSetErrorStatus(_S, NULL);
	error(E_RANGE, "mdlDerivatives");
      }
      for (i = 0; i < _mdl_nx; i++) {
	if (_mdl_x0_steady[i]) {
	  c[iy] = mdl_xp[i]/_mdl_x0_nominal[i];
	  ++iy;
	}
      }
    }

    // discrete-time state equations for parameters
    for (i = 0; i < _np; i++)
      f[i] = x[i];

    // junction conditions for continuous-time equations
    if (kk < _KK-1) {
      for (i = _np; i < _nx; i++)
	f[i] = xf[i];
    }

    // store values of model outputs
    for (i = 0; i < _mdl_ny; i++)
      _mdl_ys[kk][i] = value(mdl_y[i]);
  }

  // obtain Jacobians if required
  if (f.is_required_J() || f0.is_required_g() || c.is_required_J()) {

    // store current model states and parameters
    if (kk == 0) {
      for (i = 0; i < _mdl_np; i++)
	_mdl_p[i] = mxGetPr(_mx_p)[i];
      real_T *mdl_x = ssGetContStates(_S);
      for (i = 0; i < _mdl_nx; i++)
	_mdl_x0[i] = mdl_x[i];
    }

    // call predefined update for numerical differentiation
    _S->mdlInfo->reservedForFutureInt[0] = 1;
    Omu_Program::update_grds(kk, x, u, xf, f, f0, c);
    _S->mdlInfo->reservedForFutureInt[0] = 0;

    // correct gradient of f0 as complete nd gives bad results
    if (kk < _KK) {
#if PRG_WITH_U
      for (i = 0; i < _ny; i++) {
	f0.gu[i] = /* dt * */
	  2.0 * (u[i] - _prg_ys_ref[kk][i])
	  / _prg_y_nominal[i]/_prg_y_nominal[i];
      }
#else
      v_zero(f0.gx);
      for (i = 0; i < _ny; i++) {
	for (j = 0; j < _nx; j++)
	  f0.gx[j] += /* dt * */
	    2.0 * (c[i] - _prg_ys_ref[kk][i]) * c.Jx[i][j]
	    / _prg_y_nominal[i]/_prg_y_nominal[i];
      }
#endif
    }

    // compute M
    if (kk < _KK) {
      if (kk == 0) {
	m_zero(_dxdpx0);
	for (i = 0; i < _np; i++)
	  _dxdpx0[i][i] = 1.0;
	for (i = 0, j = _np; i < _mdl_nx; i++)
	  if (_mdl_x0_active[i])
	    _dxdpx0[_np + i][j++] = 1.0;
      }
      else {
	m_copy(_dxfdpx0, _dxdpx0);
      }
      m_mlt(c.Jx, _dxdpx0, _dxfdpx0);
      for (i = 0; i < _ny; i++)
	for (j = 0; j < _np + _nx0; j++)
	  _M[kk*_ny + i][j] = _dxfdpx0[i][j];
      if (kk < _KK - 1)
	m_mlt(xf.Sx, _dxdpx0, _dxfdpx0);
    }
  }
}

//--------------------------------------------------------------------------
void Prg_SFunctionEst::consistic(int kk, double t,
				 const Omu_StateVec &x, const Omu_Vec &u,
				 Omu_DependentVec &xt)
{
  int i, idx;

  // initialize model in first stage
  if (kk == 0 && ssGetmdlInitializeConditions(_S) != NULL) {

    // set simulation time
    ssSetT(_S, t);

    // pass estimated parameters to model
    int ip = 0;
    for (i = 0; i < _mdl_np; i++) {
      if (_mdl_p_active[i]) {
	mxGetPr(_mx_p)[i] = x[ip++] * _mdl_p_nominal[i];
      }
    }
    assert(ip == _np);	// _np must not have changed since setup_stages

    // initialize model inputs
    real_T **mdl_u = (real_T **)ssGetInputPortRealSignalPtrs(_S, 0);
    for (i = 0; i < _mdl_nu; i++)
      *mdl_u[i] = _mdl_us[kk][i];

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
void Prg_SFunctionEst::continuous(int kk, double t,
				  const Omu_StateVec &x, const Omu_Vec &u,
				  const Omu_StateVec &xp, Omu_DependentVec &F)
{
  int i;

  if (kk < _KK-1) {

    // set simulation time
    ssSetT(_S, t);

    // pass estimated parameters to model
    int ip = 0;
    for (i = 0; i < _mdl_np; i++) {
      if (_mdl_p_active[i]) {
	mxGetPr(_mx_p)[i] = x[ip++] * _mdl_p_nominal[i];
      }
    }
    assert(ip == _np);	// _np must not have changed since setup_stages

    // initialize model inputs using linear interpolation over time
    double rt = (t - ts(kk)) / (ts(kk+1) - ts(kk));
    real_T **mdl_u = (real_T **)ssGetInputPortRealSignalPtrs(_S, 0);
    for (i = 0; i < _mdl_nu; i++)
      *mdl_u[i] = _mdl_us[kk][i] * (1 - rt) + _mdl_us[kk+1][i] * rt;

    if (_mdl_np > 0 && ssGetmdlInitializeConditions(_S) != NULL) {
      // initialize model (this is required as parameters change)
      mdlInitializeConditions(_S);
      if (ssGetErrorStatus(_S)) {
	fprintf(stderr, "Error from mdlInitializeConditions: %s\n",
		ssGetErrorStatus(_S));
	ssSetErrorStatus(_S, NULL);
	error(E_RANGE, "mdlInitializeConditions");
      }
      // call mdlOutput as mdlInitializeConditions might not actually
      // perform initialization
      mdlOutputs(_S, 0);
      if (ssGetErrorStatus(_S)) {
	fprintf(stderr, "Error from mdlOutputs: %s\n", ssGetErrorStatus(_S));
	ssSetErrorStatus(_S, NULL);
	error(E_RANGE, "mdlOutputs");
      }
    }

    // pass current states to model
    real_T *mdl_x = ssGetContStates(_S);
    for (i = 0; i < _mdl_nx; i++)
      mdl_x[i] = x[_np + i] * _mdl_x0_nominal[i];

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
      F[_np + i] = mdl_xp[i]/_mdl_x0_nominal[i] - xp[_np + i];
  }

  // obtain Jacobians if required
  if (F.is_required_J()) {
    // call predefined continuous for numerical differentiation
    _S->mdlInfo->reservedForFutureInt[0] = 1;
    Omu_Program::continuous_grds(kk, t, x, u, xp, F);
    _S->mdlInfo->reservedForFutureInt[0] = 0;
  }
}


//==========================================================================

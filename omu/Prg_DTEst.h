/**
 * @file Prg_DTEst.h
 *    Estimation of initial states and parameters for a model given
 *    as Functional Model Unit (FMU) or as S-function.
 *
 * rf, 7/25/00
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

#ifndef Prg_DTEst_H
#define Prg_DTEst_H

#include <Hqp_Docp.h>
#include "Omu_Model.h"
#include "Omu_Variables.h"

/**
   Estimation of parameters and initial states for a model given
   as S-function or Functional Model Unit (FMU).
   Compared to Prg_DynamicEst, Prg_DTEst supports multi-threading for
   parallel multiple shooting and exploits the sparse structure of an FMU.
   It does not support continuous-time equations.
 */
class Prg_DTEst: public Hqp_Docp, public Omu_Model {

 protected:
  Omu_VariableVec _mdl_p;	///< model parameters (note: default min is 0)
  Omu_VariableVec _mdl_x0;	///< initial states
  Omu_VariableVec _mdl_x;	///< state bounds

  IVECP		_mdl_p_active;	///< indicate estimated parameters
  VECP		_mdl_p_confidence;///< confidence intervals for est parameters
  IVECP		_mdl_x0_active; ///< indicate estimated states
  VECP		_mdl_x0_confidence;///< confidence intervals for est x0
  IVECP 	_mdl_u_order; 	///< interpolation order (default: 1 (linear))
  IVECP		_mdl_y_active; 	///< indicate measured outputs
  VECP 		_mdl_p_nominal;	///< nominal parameter values (for scaling)

  int 		_K;	///< number of stages
  double 	_t0; 	///< initial time of optimization time horizon (default: 0)
  double 	_tf; 	///< final time of optimization time horizon (default: 1)
  VECP		_ts;	///< time steps

  int		_nx;	///< number of states for optimizer
  int		_np;	///< number of estimated parameters
  int		_nx0;	///< number of estimated initial states
  int		_ny;	///< number of reference outputs for estimation

  int		_nex;	///< number of experiments used for estimation
  MATP 		_mdl_x0s;///< initial states for each experiment
  IVECP		_exs; 	///< index identifying experiment in estimation data

  // vectors for inputs, states, and outputs
  MATP	_mdl_us;	///< given model inputs
  MATP	_mdl_xs;	///< calculated model states
  MATP	_mdl_ys;	///< calculated model outputs
  MATP	_ys_ref;	///< reference values for active outputs

  // confidence things
  double _ssr;		///< sum of squared output residuals
  MATP 	_M2;		///< measurement matrix M=dy/d(p,x0)
  MATP 	_P2;		///< precision matrix P=(M^TM)^(-1)
  MATP 	_COV;		///< covariance matrix COV=
  MATP 	_dydx;		///< help variable dy/dx
  MATP 	_dxdpx0;	///< help variable dx/d(p,x0)
  MATP 	_dydpx0;	///< help variable dy/d(p,x0)
  MATP 	_dfdx;		///< help variable df/dx
  MATP 	_dfdpx0;	///< help variable df/d(p,x0)

  bool 	_stages_ok;	///< setup_stages() was called separately
  bool  _within_grds;	///< calculating derivatives
  bool 	_ad;		///< flag about use of algorithmic differentiation
  double _fscale;	///< scaling of optimization criterion

  void setup_stages();

  /**
   * @name Implementation of predefined methods.
   * @see Hqp_Docp, Omu_Model
   */
  //@{
  void setup_model();

  void setup_horizon(int &k0, int &kf);

  void setup_vars(int k,
		  VECP x, VECP x_min, VECP x_max, IVECP x_int,
		  VECP u, VECP u_min, VECP u_max, IVECP u_int,
		  VECP c, VECP c_min, VECP c_max);

  void setup_struct(int k, const VECP x, const VECP u,
		    MATP fx, MATP fu, IVECP f_lin,
		    VECP f0x, VECP f0u, int &f0_lin,
		    MATP cx, MATP cu, IVECP c_lin,
		    MATP Lxx, MATP Luu, MATP Lxu);

  void init_simulation(int k,
		       VECP x, VECP u);

  void update_vals(int k, const VECP x, const VECP u,
		   VECP f, Real &f0, VECP c);

  void update_stage(int k, const VECP x, const VECP u,
		    VECP f, Real &f0, VECP c,
		    MATP fx, MATP fu, VECP f0x, VECP f0u,
		    MATP cx, MATP cu,
		    const VECP rf, const VECP rc,
		    MATP Lxx, MATP Luu, MATP Lxu);
  //@}

  /// obtain analytic Jacobian wrt model states and inputs
  void fetch_jacxu(SimStruct *S,
                   int k, const VECP x, const VECP u,
                   MATP fx, MATP fu, MATP cx, MATP cu);

  /// write active parameters that are packed in p to _mx_args
  void write_active_mx_args(VECP p);

 public:

  Prg_DTEst();		///< constructor
  ~Prg_DTEst();		///< destructor

  const char *name() {return "DTEst";} ///< name DTEst

  /**
   * @name Access methods for program specific members (If prefix: prg_)
   */
  //@{

  int K() const {return _K;}	///< get number of stages
  void set_K(int K) {_K = K;}	///< set number of stages
  double t0() const {return _t0;}///< get initial time of horizon
  /// set initial time
  void set_t0(double t0) {
    _t0 = t0;
    if (_ts->dim > 0)
      _ts[0] = t0;
  }
  double tf() const {return _tf;}///< get final time of horizon
  /// set final time
  void set_tf(double tf) {
    _tf = tf;
    if (_ts->dim > 0)
      _ts[(int)_ts->dim - 1] = tf;
  }
  /// vector of start time points in each sample period
  const VECP ts() const {return _ts;}
  /// set vector of start times
  void set_ts(const VECP n_ts) {
    v_copy_elements(n_ts, _ts);
    if (_ts->dim > 0) {
      _t0 = _ts[0];
      _tf = _ts[(int)_ts->dim - 1];
    }
  }	
  /// get start time point of sample period k
  double ts(int k) const {return _ts[k];}
  /// scaling of optimization criterion

  /// flag about use of automatic differentiation
  bool ad() const {return _ad;}
  /// Set flag about use of automatic differentiation.
  void set_ad(bool val) {_ad = val;}

  double fscale() const {return _fscale;}
  /// Set scaling of optimization criterion.
  void set_fscale(double val) {_fscale = val;}

  /// number of experiments used for estimation
  int nex() const {return _nex;}
  /// set number of experiments
  void set_nex(int val) {_nex = val;}

  /// reference values for active model outputs (size: KK+1 . ny)
  const MATP ys_ref() const {return _ys_ref;}
  /// set reference outputs
  void set_ys_ref(const MATP val) {m_copy_elements(val, _ys_ref);}

  /// measurement matrix M=dy/d(p,x0) (size: (KK+1)*ny . np+nx0)
  const MATP M() const {return _M2;}

  /// precision matrix P=(M^TM)^(-1) (size: np+nx0 . np+nx0)
  const MATP P() const {return _P2;}

  /// covariance matrix (size: np+nx0 . np+nx0)
  const MATP COV() const {return _COV;}

  //@}

  /**
   * @name Read methods for model specific members (no If prefix).
   * Note that the prefix mdl_ is omitted in the detailed mathematical
   * problem description.
   */
  //@{
  /// model parameters 
  const VECP mdl_p() const {return _mdl_p;}

  /// lower bounds for model parameters (note: default is 0)
  const VECP mdl_p_min() const {return _mdl_p.min;}

  /// upper bounds for model parameters
  const VECP mdl_p_max() const {return _mdl_p.max;}

  /// indicate estimated parameters
  const IVECP mdl_p_active() const {return _mdl_p_active;}

  /// confidence intervals of estimated parameters
  const VECP mdl_p_confidence() const {return _mdl_p_confidence;}

  /// nominal parameter values (for scaling)
  const VECP mdl_p_nominal() const {return _mdl_p_nominal;}

  /// model initial states
  /// (read only, note: write initial states of experiments via mdl_x0s)
  const VECP mdl_x0() const {return _mdl_x0;}

  /// lower bounds for initial states
  const VECP mdl_x0_min() const {return _mdl_x0.min;}

  /// upper bounds for initial states
  const VECP mdl_x0_max() const {return _mdl_x0.max;}

  /// indicate estimated initial states
  const IVECP mdl_x0_active() const {return _mdl_x0_active;}

  /// confidence intervals of estimated initial states,
  /// averaged over multiple experiments
  const VECP mdl_x0_confidence() const {return _mdl_x0_confidence;}

  ///< interpolation order (0 (constant) or 1 (linear), default: 1)
  const IVECP mdl_u_order() const {return _mdl_u_order;}

  /// lower bounds on states at stage boundaries (default: -Inf)
  const VECP mdl_x_min() const {return _mdl_x.min;}

  /// upper bounds on states at stage boundaries (default: Inf)
  const VECP mdl_x_max() const {return _mdl_x.max;}

  /// indicate measured outputs
  const IVECP mdl_y_active() const {return _mdl_y_active;}

  /// initial states for each experiment (size: nex . mdl_nx)
  const MATP mdl_x0s() const {return _mdl_x0s;}

  /// model inputs (size: KK+1 . mdl_nu)
  const MATP mdl_us() const {return _mdl_us;}

  /// model states (size: KK+1 . mdl_nx)
  const MATP mdl_xs() const {return _mdl_xs;}

  /// model outputs (read only, size: KK+1 . mdl_ny)
  const MATP mdl_ys() const {return _mdl_ys;}

  //@}

  /**
   * @name Write methods for model specific members (no If prefix).
   */
  //@{
  /// set values of model parameter
  void set_mdl_p(const VECP v) {v_copy_elements(v, _mdl_p);}

  /// set lower bounds for model parameters
  void set_mdl_p_min(const VECP v) {v_copy_elements(v, _mdl_p.min);}

  /// set upper bounds for model parameters
  void set_mdl_p_max(const VECP v) {v_copy_elements(v, _mdl_p.max);}

  /// set estimated parameters
  void set_mdl_p_active(const IVECP v) {iv_copy_elements(v, _mdl_p_active);}

  // confidence intervals are read only

  /// set nominal parameters
  void set_mdl_p_nominal(const VECP v) {v_copy_elements(v, _mdl_p_nominal);}

  // model initial states are read only

  /// set lower bounds for initial states
  void set_mdl_x0_min(const VECP v) {v_copy_elements(v, _mdl_x0.min);}

  /// set upper bounds for initial states
  void set_mdl_x0_max(const VECP v) {v_copy_elements(v, _mdl_x0.max);}

  /// set estimated initial states
  void set_mdl_x0_active(const IVECP v) {iv_copy_elements(v, _mdl_x0_active);}

  // confidence intervals are read only

  /// set interpolation order for inputs
  void set_mdl_u_order(const IVECP v) {iv_copy_elements(v, _mdl_u_order);}

  /// set lower bounds on states
  void set_mdl_x_min(const VECP v)
  {v_copy_elements(v, _mdl_x.min);}

  /// set upper bounds on states
  void set_mdl_x_max(const VECP v)
  {v_copy_elements(v, _mdl_x.max);}

  /// set measured outputs
  void set_mdl_y_active(const IVECP v) {iv_copy_elements(v, _mdl_y_active);}

  /// set initial states
  void set_mdl_x0s(const MATP v) {
    m_copy_elements(v, _mdl_x0s);
    // copy initial states to all sample periods of all experiments
    // to achieve multiple shooting behavior when no simulation is performed
    int ex = 0;
    for (int k = 0; k <= _K; k++) {
      if (k > 0 && _ts[k] < _ts[k-1])
	ex++;
      for (int i = 0; i < _mdl_nx; i++)
	_mdl_xs[k][i] = v[ex][i];
    }
  }

  /// set model inputs
  void set_mdl_us(const MATP v) {m_copy_elements(v, _mdl_us);}

  /// set model states
  void set_mdl_xs(const MATP v) {
    m_copy_elements(v, _mdl_xs);
    // also set additionally treated mdl_x0
    for (int i = 0; i < _mdl_nx; i++)
      _mdl_x0[i] = _mdl_xs[0][i];
    // also set additionally treated mdl_x0s
    int ex = 0;
    for (int k = 0; k <= _K; k++) {
      if (k == 0 || _ts[k] < _ts[k-1]) {
	for (int i = 0; i < _mdl_nx; i++)
	  _mdl_x0s[ex][i] = _mdl_xs[k][i];
	ex++;
      }
    }
  }

  // model outputs are read only

  //@}
};  

#endif

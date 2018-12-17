/**
 * @file Prg_DTOpt.h
 *    Discrete-time Optimal Control Problem for a model given
 *    as Functional Model Unit (FMU) or as S-function.
 *
 * rf, 5/1/17
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

#ifndef Prg_DTOpt_H
#define Prg_DTOpt_H

#include <Hqp_Docp.h>
#include "Omu_Model.h"
#include "Omu_Variables.h"

/**
   Discrete-time optimal control program for a model given as S-function or 
   Functional Model Unit (FMU). The program is treated with
   multi-stage control vector parameterization.
   Compared to Prg_DynamicOpt, Prg_DTOpt supports multi-threading for
   parallel multiple shooting and exploits the sparse structure of an FMU.
   It does not support continuous-time equations and only covers one
   sample period per stage / shooting interval.
 */
class Prg_DTOpt: public Hqp_Docp, public Omu_Model {

 protected:
  Omu_VariableVec _mdl_x0;	///< initial states for optimization
  IVECP 	_mdl_x0_active;	///< free initial states (default: 0)
  Omu_VariableVec _mdl_u0;	///< initial inputs
  Omu_OptVarVec _mdl_y0; 	///< model outputs at initial time
  Omu_OptVarVec _mdl_y0_soft; 	///< relaxed model outputs at initial time
  Omu_OptVarVec _mdl_u; 	///< model inputs
  Omu_OptVarVec _mdl_der_u; 	///< rates of change of inputs
  Omu_OptVarVec _mdl_der_u_soft;///< relaxed rates of change of inputs
  Omu_VariableVec _mdl_x;	///< state bounds
  Omu_OptVarVec _mdl_y; 	///< model outputs
  Omu_OptVarVec _mdl_y_soft; 	///< attributes for relaxed output constraints
  Omu_OptVarVec _mdl_uf; 	///< model inputs at final time
  Omu_OptVarVec _mdl_yf; 	///< model outputs at final time
  Omu_OptVarVec _mdl_yf_soft; 	///< relaxed output constraints at final time

  IVECP 	_mdl_u_order; 	///< interpolation order (default: 1 (linear))
  IVECP 	_mdl_y_order; 	///< interpolation order (default: 1 (linear))

  int 		_K;	///< number of stages
  double 	_t0; 	///< initial time of optimization time horizon (default: 0)
  double 	_tf; 	///< final time of optimization time horizon (default: 1)
  VECP		_ts;	        ///< time steps
  VECP		_taus;		///< scaled time communicated to outside as prg_ts
  int 		_t_scale_idx; 	///< index into mdl_u vector for variable used for scaling of time
  int 		_t_active; 	///< time is being scaled
  int 		_t_scale_i; 	///< index of internal optimization variable for scaling of time
  double 	_t_scale_nominal;///< nominal value for time scaling (default: 1)
  double 	_t_nominal; 	///< nominal time (used internally for scaling)

  VECP 		_mdl_y_bias;	///< bias correction (offset) for outputs 
  VECP 		_mdl_y_lambda;	///< max Lagrange multiplier for constraints
  VECP 		_c_lambda;	///< Lagrange multipliers of constraints

  int		_nx;	///< number of states for optimizer
  int		_nu;	///< number of optimized control inputs
  int		_ndu;	///< number of optimized inputs with rate of change
  int		_nc;	///< number of constrained outputs
  int		_ns;	///< number of slack variables for soft constraints
  int		_nsc;	///< number of soft constraints
  int		_nsu;	///< number of slack variables for soft cns on inputs
  int		_nsuc;	///< number of soft constraints on inputs
  int		_nc0;	///< number of constrained/used outputs at initial time
  int		_ns0;	///< number of slacks for soft constraints at initial time
  int		_nsc0;	///< number of soft constraints at initial time
  int		_ncf;	///< number of constrained/used outputs at final time
  int		_nsf;	///< number of slacks for soft constraints at final time
  int		_nscf;	///< number of soft constraints at final time

  MATP	_mdl_us;	///< given model inputs (controls and disturbances)
  MATP	_mdl_xs;	///< given and calculated model states
  MATP	_mdl_ys;	///< calculated model outputs

  /**
   * Numbers of fixed control inputs at begin of time horizon (default: 0).
   * An active input is free for _mdl_u0_fixed=0, fixed at the initial time
   * point for _mdl_u0_nfixed=1, fixed up to the last time point of the
   * first stage for _mdl_u0_nfixed=2, and so on.
   */
  IVECP 	_mdl_u0_nfixed;

  /**
   * Decimation factor, i.e.~numbers of subsequent optimized inputs with 
   * equal values (default: 1). The controls of two subsequent stages are
   * equal for decimation=2, three for decimation=3, and so on.
   */
  IVECP 	_mdl_u_decimation;

  /**
   * Periodic controls, i.e.~first and last value are equal
   */
  IVECP 	_mdl_u_periodic;

  /**
   * Periodic states, i.e.~first and last value are equal
   */
  IVECP 	_mdl_x_periodic;

  bool 	_stages_ok;	///< setup_stages() was called separately
  bool 	_within_grds;	///< calculating derivatives
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

  void fetch_jac(SimStruct *S,
                 int k, double tscale, const VECP x, const VECP u,
                 MATP fx, MATP fu, MATP cx, MATP cu);

 public:

  Prg_DTOpt();		///< constructor
  ~Prg_DTOpt();		///< destructor

  const char *name() {return "DTOpt";} ///< name DTOpt

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
  //@}

  /**
   * @name Access methods for scaling of time
   */
  //@{
  /// optional index of a model input used for scaling of time (default: -1)
  int mdl_t_scale_idx() const {return _t_scale_idx;}
  /// set index of model input used for scaling of time
  void set_mdl_t_scale_idx(int val) {_t_scale_idx = val;}

  /// vector of start time points in each sample period
  const VECP taus() const {
    return _taus;
  }
  /// set vector of start times
  void set_taus(const VECP n_taus) {
    // first use _taus as intermediate storage for initialization of _ts
    v_copy_elements(n_taus, _taus); // this also checks the dimension of n_taus
    if (_t_scale_idx >= 0) {
      for (int kk = 1; kk < (int)_taus->dim; kk++)
        _taus[kk] = _taus[kk-1] + ((n_taus[kk] - n_taus[kk-1])
                                   / _mdl_us[kk][_t_scale_idx]);
    }
    // take over unscaled ts
    set_ts(_taus);
    // finally take over n_taus
    if (_t_scale_idx >= 0)
      v_copy_elements(n_taus, _taus);
  }	
  //@}

  /**
   * @name Read methods for model specific members (no If prefix).
   * Note that the prefix mdl_ is omitted in the description of the
   * detailed mathematical program.
   */
  //@{
  /// initial states
  const VECP mdl_x0() const {return _mdl_x0;}

  /// free initial states (default: 0)
  const IVECP mdl_x0_active() const {return _mdl_x0_active;}

  /// lower bounds for initial states
  const VECP mdl_x0_min() const {return _mdl_x0.min;}

  /// upper bounds for initial states
  const VECP mdl_x0_max() const {return _mdl_x0.max;}

  /// model inputs at initial time
  const VECP mdl_u0() const {return _mdl_u0;}

  /// lower bounds for model inputs at initial time
  const VECP mdl_u0_min() const {return _mdl_u0.min;}

  /// upper bounds for model inputs at initial time
  const VECP mdl_u0_max() const {return _mdl_u0.max;}

  /// model outputs at initial time
  const VECP mdl_y0() const {return _mdl_y0;}

  /// lower bounds for model outputs at initial time
  const VECP mdl_y0_min() const {return _mdl_y0.min;}

  /// upper bounds for model outputs at initial time
  const VECP mdl_y0_max() const {return _mdl_y0.max;}

  /// weight for linear objective term at initial time (default: 0)
  const VECP mdl_y0_weight1() const {return _mdl_y0.weight1;}

  /// weight for quadratic objective term at initial time (default: 0)
  const VECP mdl_y0_weight2() const {return _mdl_y0.weight2;}

  /// soft lower bounds for model outputs at initial time
  const VECP mdl_y0_soft_min() const {return _mdl_y0_soft.min;}

  /// soft upper bounds for model outputs at initial time
  const VECP mdl_y0_soft_max() const {return _mdl_y0_soft.max;}

  /// soft weight for linear objective term at initial time (default: 0)
  const VECP mdl_y0_soft_weight1() const {return _mdl_y0_soft.weight1;}

  /// soft weight for quadratic objective term at initial time (default: 0)
  const VECP mdl_y0_soft_weight2() const {return _mdl_y0_soft.weight2;}

  /// interpolation order (0 (constant) or 1 (linear), default: 1)
  const IVECP mdl_u_order() const {return _mdl_u_order;}

  /// indicate optimized inputs
  const IVECP mdl_u_active() const {return _mdl_u.active;}

  /// indicate integers
  const IVECP mdl_u_integer() const {return _mdl_u.integer;}

  /// numbers of fixed control inputs at begin of time horizon (default: 1)
  const IVECP mdl_u0_nfixed() const {return _mdl_u0_nfixed;}

  /// decimation for optimized model inputs (default: 1)
  const IVECP mdl_u_decimation() const {return _mdl_u_decimation;}

  /// first and last value are equal (default: false)
  const IVECP mdl_u_periodic() const {return _mdl_u_periodic;}

  /// lower bounds for optimized model inputs
  const VECP mdl_u_min() const {return _mdl_u.min;}

  /// upper bounds for optimized model inputs
  const VECP mdl_u_max() const {return _mdl_u.max;}

  /// reference values for optimized model inputs (default: 0)
  const VECP mdl_u_ref() const {return _mdl_u.ref;}

  /// weight for linear objective term (default: 0)
  const VECP mdl_u_weight1() const {return _mdl_u.weight1;}

  /// weight for quadratic objective term (default: 0)
  const VECP mdl_u_weight2() const {return _mdl_u.weight2;}

  /// lower bounds for rates of change of optimized model inputs
  const VECP mdl_der_u_min() const {return _mdl_der_u.min;}

  /// upper bounds for rates of change of optimized model inputs
  const VECP mdl_der_u_max() const {return _mdl_der_u.max;}

  /// reference values for rates of change of optimized inputs (default: 0)
  const VECP mdl_der_u_ref() const {return _mdl_der_u.ref;}

  /// weight for linear objective term (default: 0)
  const VECP mdl_der_u_weight1() const {return _mdl_der_u.weight1;}

  /// weight for quadratic objective term (default: 0)
  const VECP mdl_der_u_weight2() const {return _mdl_der_u.weight2;}

  /// soft lower bounds for rates of change of optimized model inputs
  const VECP mdl_der_u_soft_min() const {return _mdl_der_u_soft.min;}

  /// soft upper bounds for rates of change of optimized model inputs
  const VECP mdl_der_u_soft_max() const {return _mdl_der_u_soft.max;}

  /// weight for linear objective term (default: 0)
  const VECP mdl_der_u_soft_weight1() const {return _mdl_der_u_soft.weight1;}

  /// weight for quadratic objective term (default: 0)
  const VECP mdl_der_u_soft_weight2() const {return _mdl_der_u_soft.weight2;}

  /// indicate integer valued states
  const IVECP mdl_x_integer() const {return _mdl_x.integer;}

  /// first and last value are equal (default: false)
  const IVECP mdl_x_periodic() const {return _mdl_x_periodic;}

  /// lower bounds on states at stage boundaries (default: -Inf)
  const VECP mdl_x_min() const {return _mdl_x.min;}

  /// upper bounds on states at stage boundaries (default: Inf)
  const VECP mdl_x_max() const {return _mdl_x.max;}

  /// interpolation order (0 (constant) or 1 (linear), default: 1)
  const IVECP mdl_y_order() const {return _mdl_y_order;}

  /// bias for outputs
  const VECP mdl_y_bias() const {return _mdl_y_bias;}

  /// max Lagrange multiplier of all constraints for the respective output
  const VECP mdl_y_lambda() const {return _mdl_y_lambda;}

  /// lower bounds for model outputs
  const VECP mdl_y_min() const {return _mdl_y.min;}

  /// upper bounds for model outputs
  const VECP mdl_y_max() const {return _mdl_y.max;}

  /// reference values for model outputs (default: 0)
  const VECP mdl_y_ref() const {return _mdl_y.ref;}

  /// weight for linear objective term (default: 0)
  const VECP mdl_y_weight1() const {return _mdl_y.weight1;}

  /// weight for quadratic objective term (default: 0)
  const VECP mdl_y_weight2() const {return _mdl_y.weight2;}

  /// soft lower bounds for model outputs
  const VECP mdl_y_soft_min() const {return _mdl_y_soft.min;}

  /// soft upper bounds for model outputs
  const VECP mdl_y_soft_max() const {return _mdl_y_soft.max;}

  /// weight for linear objective term (default: 0)
  const VECP mdl_y_soft_weight1() const {return _mdl_y_soft.weight1;}

  /// weight for quadratic objective term (default: 0)
  const VECP mdl_y_soft_weight2() const {return _mdl_y_soft.weight2;}

  /// model inputs at final time
  const VECP mdl_uf() const {return _mdl_uf;}

  /// lower bounds for model inputs at final time
  const VECP mdl_uf_min() const {return _mdl_uf.min;}

  /// upper bounds for model inputs at final time
  const VECP mdl_uf_max() const {return _mdl_uf.max;}

  /// model outputs at final time
  const VECP mdl_yf() const {return _mdl_yf;}

  /// lower bounds for model outputs at final time
  const VECP mdl_yf_min() const {return _mdl_yf.min;}

  /// upper bounds for model outputs at final time
  const VECP mdl_yf_max() const {return _mdl_yf.max;}

  /// weight for linear objective term (default: 0)
  const VECP mdl_yf_weight1() const {return _mdl_yf.weight1;}

  /// weight for quadratic objective term (default: 0)
  const VECP mdl_yf_weight2() const {return _mdl_yf.weight2;}

  /// soft lower bounds for model outputs at final time
  const VECP mdl_yf_soft_min() const {return _mdl_yf_soft.min;}

  /// soft upper bounds for model outputs at final time
  const VECP mdl_yf_soft_max() const {return _mdl_yf_soft.max;}

  /// weight for linear objective term (default: 0)
  const VECP mdl_yf_soft_weight1() const {return _mdl_yf_soft.weight1;}

  /// weight for quadratic objective term (default: 0)
  const VECP mdl_yf_soft_weight2() const {return _mdl_yf_soft.weight2;}

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
  /// set initial states and copy them to all mdl_xs
  void set_mdl_x0(const VECP v) {
    v_copy_elements(v, _mdl_x0);
    // copy initial states to all sample periods
    // to achieve multiple shooting behavior when no simulation is performed
    for (int k = 0; k <= _K; k++)
      for (int i = 0; i < _mdl_nx; i++)
	_mdl_xs[k][i] = v[i];
  }

  /// set free initial states
  void set_mdl_x0_active(const IVECP v) {iv_copy_elements(v, _mdl_x0_active);}

  /// set lower bounds for initial states
  void set_mdl_x0_min(const VECP v) {v_copy_elements(v, _mdl_x0.min);}

  /// set upper bounds for initial states
  void set_mdl_x0_max(const VECP v) {v_copy_elements(v, _mdl_x0.max);}

  /// set initial values for model inputs and copy them to all mdl_us
  void set_mdl_u0(const VECP v) {
    v_copy_elements(v, _mdl_u0);
    for (int k = 0; k <= _K; k++)
      for (int i = 0; i < _mdl_nu; i++)
	_mdl_us[k][i] = v[i];
  }

  /// set lower bounds for initial model inputs
  void set_mdl_u0_min(const VECP v) {v_copy_elements(v, _mdl_u0.min);}

  /// set upper bounds for initial model inputs
  void set_mdl_u0_max(const VECP v) {v_copy_elements(v, _mdl_u0.max);}

  /// set lower bounds for initial model outputs
  void set_mdl_y0_min(const VECP v) {v_copy_elements(v, _mdl_y0.min);}

  /// set upper bounds for initial model outputs
  void set_mdl_y0_max(const VECP v) {v_copy_elements(v, _mdl_y0.max);}

  /// set linear weight
  void set_mdl_y0_weight1(const VECP v) {v_copy_elements(v, _mdl_y0.weight1);}

  /// set quadratic weight
  void set_mdl_y0_weight2(const VECP v) {v_copy_elements(v, _mdl_y0.weight2);}

  /// set soft lower bounds for model outputs at initial time
  void set_mdl_y0_soft_min(const VECP v) {v_copy_elements(v, _mdl_y0_soft.min);}

  /// set soft upper bounds for model outputs at initial time
  void set_mdl_y0_soft_max(const VECP v) {v_copy_elements(v, _mdl_y0_soft.max);}

  /// set soft weight for linear objective term at initial time (default: 0)
  void set_mdl_y0_soft_weight1(const VECP v) {v_copy_elements(v, _mdl_y0_soft.weight1);}

  /// set soft weight for quadratic objective term at initial time (default: 0)
  void set_mdl_y0_soft_weight2(const VECP v) {v_copy_elements(v, _mdl_y0_soft.weight2);}

  /// set interpolation order
  void set_mdl_u_order(const IVECP v) {iv_copy_elements(v, _mdl_u_order);}

  /// set optimized inputs
  void set_mdl_u_active(const IVECP v) {iv_copy_elements(v, _mdl_u.active);}

  /// set integers
  void set_mdl_u_integer(const IVECP v) {iv_copy_elements(v, _mdl_u.integer);}

  /// set numbers of fixed control inputs
  void set_mdl_u0_nfixed(const IVECP v) {iv_copy_elements(v, _mdl_u0_nfixed);}

  /// set decimation for optimized model inputs
  void set_mdl_u_decimation(const IVECP v)
  {iv_copy_elements(v, _mdl_u_decimation);}

  /// set periodic congtrols
  void set_mdl_u_periodic(const IVECP v) {iv_copy_elements(v, _mdl_u_periodic);}

  /// set lower bounds for model inputs
  void set_mdl_u_min(const VECP v) {v_copy_elements(v, _mdl_u.min);}

  /// set upper bounds for model inputs
  void set_mdl_u_max(const VECP v) {v_copy_elements(v, _mdl_u.max);}

  /// set reference model inputs
  void set_mdl_u_ref(const VECP v) {v_copy_elements(v, _mdl_u.ref);}

  /// set linear weight
  void set_mdl_u_weight1(const VECP v) {v_copy_elements(v, _mdl_u.weight1);}

  /// set quadratic weight
  void set_mdl_u_weight2(const VECP v) {v_copy_elements(v, _mdl_u.weight2);}

  /// set lower bounds for rates of change of model inputs
  void set_mdl_der_u_min(const VECP v) {v_copy_elements(v, _mdl_der_u.min);}

  /// set upper bounds for rates of change of model inputs
  void set_mdl_der_u_max(const VECP v) {v_copy_elements(v, _mdl_der_u.max);}

  /// set reference rates of change of model inputs
  void set_mdl_der_u_ref(const VECP v) {v_copy_elements(v, _mdl_der_u.ref);}

  /// set linear weight
  void set_mdl_der_u_weight1(const VECP v)
  {v_copy_elements(v, _mdl_der_u.weight1);}

  /// set quadratic weight
  void set_mdl_der_u_weight2(const VECP v)
  {v_copy_elements(v, _mdl_der_u.weight2);}

  /// set soft lower bounds for rates of change of model inputs
  void set_mdl_der_u_soft_min(const VECP v)
  {v_copy_elements(v, _mdl_der_u_soft.min);}

  /// set soft upper bounds for rates of change of model inputs
  void set_mdl_der_u_soft_max(const VECP v)
  {v_copy_elements(v, _mdl_der_u_soft.max);}

  /// set linear weight
  void set_mdl_der_u_soft_weight1(const VECP v)
  {v_copy_elements(v, _mdl_der_u_soft.weight1);}

  /// set quadratic weight
  void set_mdl_der_u_soft_weight2(const VECP v)
  {v_copy_elements(v, _mdl_der_u_soft.weight2);}

  /// set integer valued states
  void set_mdl_x_integer(const IVECP v) {iv_copy_elements(v, _mdl_x.integer);}

  /// set periodic states
  void set_mdl_x_periodic(const IVECP v) {iv_copy_elements(v, _mdl_x_periodic);}

  /// set lower bounds on states
  void set_mdl_x_min(const VECP v)
  {v_copy_elements(v, _mdl_x.min);}

  /// set upper bounds on states
  void set_mdl_x_max(const VECP v)
  {v_copy_elements(v, _mdl_x.max);}

  /// set interpolation order
  void set_mdl_y_order(const IVECP v) {iv_copy_elements(v, _mdl_y_order);}

  /// set output bias
  void set_mdl_y_bias(const VECP v) {v_copy_elements(v, _mdl_y_bias);}

  /// set lower bounds for model outputs
  void set_mdl_y_min(const VECP v) {v_copy_elements(v, _mdl_y.min);}

  /// set upper bounds for model outputs
  void set_mdl_y_max(const VECP v) {v_copy_elements(v, _mdl_y.max);}

  /// set reference model outputs
  void set_mdl_y_ref(const VECP v) {v_copy_elements(v, _mdl_y.ref);}

  /// set linear weight
  void set_mdl_y_weight1(const VECP v) {v_copy_elements(v, _mdl_y.weight1);}

  /// set quadratic weight
  void set_mdl_y_weight2(const VECP v) {v_copy_elements(v, _mdl_y.weight2);}

  /// set soft lower bounds for model outputs
  void set_mdl_y_soft_min(const VECP v) {v_copy_elements(v, _mdl_y_soft.min);}

  /// set soft upper bounds for model outputs
  void set_mdl_y_soft_max(const VECP v) {v_copy_elements(v, _mdl_y_soft.max);}

  /// set linear weight
  void set_mdl_y_soft_weight1(const VECP v)
  {v_copy_elements(v, _mdl_y_soft.weight1);}

  /// set quadratic weight
  void set_mdl_y_soft_weight2(const VECP v)
  {v_copy_elements(v, _mdl_y_soft.weight2);}

  /// set lower bounds for final model inputs
  void set_mdl_uf_min(const VECP v) {v_copy_elements(v, _mdl_uf.min);}

  /// set upper bounds for final model inputs
  void set_mdl_uf_max(const VECP v) {v_copy_elements(v, _mdl_uf.max);}

  /// set lower bounds for final model outputs
  void set_mdl_yf_min(const VECP v) {v_copy_elements(v, _mdl_yf.min);}

  /// set upper bounds for final model outputs
  void set_mdl_yf_max(const VECP v) {v_copy_elements(v, _mdl_yf.max);}

  /// set linear weight
  void set_mdl_yf_weight1(const VECP v) {v_copy_elements(v, _mdl_yf.weight1);}

  /// set quadratic weight
  void set_mdl_yf_weight2(const VECP v) {v_copy_elements(v, _mdl_yf.weight2);}

  /// set soft lower bounds for model outputs at final time
  void set_mdl_yf_soft_min(const VECP v)
  {v_copy_elements(v, _mdl_yf_soft.min);}

  /// set soft upper bounds for model outputs at final time
  void set_mdl_yf_soft_max(const VECP v)
  {v_copy_elements(v, _mdl_yf_soft.max);}

  /// set linear weight
  void set_mdl_yf_soft_weight1(const VECP v)
  {v_copy_elements(v, _mdl_yf_soft.weight1);}

  /// set quadratic weight
  void set_mdl_yf_soft_weight2(const VECP v)
  {v_copy_elements(v, _mdl_yf_soft.weight2);}

  /// set model inputs
  void set_mdl_us(const MATP v) {
    m_copy_elements(v, _mdl_us);
    // also set additionally treated mdl_u0
    for (int i = 0; i < _mdl_nu; i++)
      _mdl_u0[i] = _mdl_us[0][i];
  }

  /// set model states
  void set_mdl_xs(const MATP v) {
    m_copy_elements(v, _mdl_xs);
    // also set additionally treated mdl_x0
    for (int i = 0; i < _mdl_nx; i++)
      _mdl_x0[i] = _mdl_xs[0][i];
  }

  /// model outputs are read only

  //@}
};  

#endif

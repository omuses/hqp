/**
 * @file Prg_SFunctionEst.h
 *    Estimation of initial states and parameters for a model given
 *    as MEX S-function.
 *
 * rf, 7/25/00
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

#ifndef Prg_SFunctionEst_H
#define Prg_SFunctionEst_H

#include "Prg_SFunction.h"
#include "Omu_Variables.h"

/**
   Estimation of parameters and initial states for a model given
   as MEX S-function. The estimation problem is formulated 
   at the discrete time points @f$t^{0,l}<t^{1,l}<\ldots<t^{K_l}@f$ with
   @f$l=1,\ldots,N_{ex}@f$ a number of experiments. The data of all
   experiments is concatenated into combined vectors of time points,
   model inputs, and reference outputs, i.e. @f$K_l+1=k_{0,l+1}@f$.
   In order to distinguish multiple experiments, it is assumed that
   @f$t^{K_l}>t^{0,l+1}@f$; for instance each experiment could start
   with @f$t^{0,l}=0@f$. Parameter values are constrained to be the
   same for all experiments, whereas initial states are individual for
   each experiment. Note that this is not a restriction as generally
   states can be initialized with parameters and as parameters can be
   modeled as constant states.

   In the following all vector operations are defined element wise.
   The treated estimation problem reads
   @f[
   \begin{array}{l}
    \displaystyle{}
    J\ =\ \sum_{l=1}^{N_{ex}} \sum_{k_l=k_{l,0}}^{K_l} \sum_{i=1}^{dim(I)} 
      \left\{
      \left[\frac{ys^{k_l}(I)-ys^{k_l}_{ref}}{y_{nominal}(I)}\right]^2
      \right\}_i
      \quad\to\quad \min, \quad I = \mbox{find}(y_{active}),
   \end{array}
   @f]
   subject to the model given with the S-function methods
   mdlDerivatives @f$f@f$ and mdlOutputs  @f$g@f$
   @f[
   \begin{array}{l}
    \displaystyle \dot{x}(t) = 
     \frac{f[p,\ x_{nominal}\,x(t),\ u(t)]}
          {x_{nominal}}, \\[3ex]
    \displaystyle y(t) = 
     \frac{g[p,\ x_{nominal}\,x(t),\ u(t)]}
          {y_{nominal}}, \\[3ex]
    ys^{k_l} = y_{nominal}y(t^{k_l}),
        \quad k_l=k_{l,0},\ldots,K_l, \quad l=1,\ldots,N_{ex}
   \end{array}
   @f]
   with piecewise linear interpolation of the inputs
   @f[
   \begin{array}{l}
      u(t) = \displaystyle \frac{t^{k_l+1}-t}{t^{k_l+1}-t^{k_l}}\ us^{k_l} 
             + \frac{t-t^{k_l}}{t^{k_l+1}-t^{k_l}}\ us^{k_l+1},
       \quad t\in[t^{k_l},t^{k_l+1}), \quad k_l=k_{l,0},\ldots,K_l-1,
       \quad l=1,\ldots,N_{ex},
   \end{array}
   @f]
   and subject to the constraints
   @f[
   \begin{array}{rcccll}
    \displaystyle \left\{ \frac{p_{min}}{p_{nominal}} \right. 
        &<& \displaystyle \frac{p}{p_{nominal}}
        &<& \left. \displaystyle \frac{p_{max}}{p_{nominal}} \right\}_i,
        \quad & i \in \mbox{find}(p_{active}), \\[3ex]
    \displaystyle \left\{ \frac{x^0_{min}}{x_{nominal}} \right. &<& x(t^{0,l}) 
        &<& \displaystyle \left. \frac{x^0_{max}}{x_{nominal}} \right\}_i, 
        \quad & i \in \mbox{find}(x^0_{active}),
        \quad l=1,\ldots,N_{ex},  \\[3ex]
    && \{ x(t^{0,l}) &=& \displaystyle \frac{x0s^l}{x_{nominal}}
       \}_i, \quad & i \notin \mbox{find}(x^0_{active}),
        \quad l=1,\ldots,N_{ex},  \\[3ex]
    \displaystyle \frac{der\_x^0_{min}}{x_{nominal}} 
        &<& \dot{x}(t^{0,l}) 
        &<& \displaystyle \frac{der\_x^0_{max}}{x_{nominal}},
	\quad & l=1,\ldots,N_{ex}.
   \end{array}
   @f]
   The problem is treated as multistage problem with @f$K_{N_{ex}}@f$ stages
   and one sample period per stage per default. 
   Discrete-time state variables with unknown initial value are introduced
   for the estimated parameters p, in order to preserve a sparse structure.
   Additional @f$K_{N_{ex}}@f$ junction conditions (equality constraints)
   are introduced for the estimated parameters and @f$K_{N_{ex}}-N_{ex}+1@f$
   junction conditions are introduced for the state variables x.
   Alternatively the problem can be treated without stages hiding model
   states from the optimizer.

   In addition to the estimated parameters and initial states, the 
   measurement matrix M is calculated in each optimization iteration.
   It is defined as
   @f[
   M = \left[\begin{array}{cc}
       M_p^0, & M_{x^0}^0 \\[2ex]
       M_p^1, & M_{x^0}^1 \\[2ex]
       \vdots, & \vdots \\[2ex]
       M_p^{K_{N_{ex}}}, & M_{x^0}^{K_{N_{ex}}}
       \end{array}\right]
   @f]
   with the sub-matrices
   @f[
   \begin{array}{r}
     \displaystyle M_p^{k,l} = \frac{d\,\displaystyle\frac{ys^{k,l}}
                                    {y_{nominal}}(I_y)}
                 {d\,\displaystyle\frac{p}{p_{nominal}}(I_p)}, \quad
     M_{x^0}^{k,l} = \frac{d\,\displaystyle\frac{ys^{k,l}}{y_{nominal}}(I_y)}
                   {d\,\displaystyle\frac{x0s^l}{x_{nominal}}(I_{x^0})}, \quad
     k_l=k_{l,0},\ldots,K_l, \quad l=1,\ldots,N_{ex}, \\[7ex]
     I_y = \mbox{find}(y_{active}),\ I_p = \mbox{find}(p_{active}),
     \ I_{x^0} = \mbox{find}(x^0_{active}).
   \end{array}
   @f]
   The measurement matrix can be used to obtain confidences for estimation
   results. Note that the sub-matrices for estimated initial states
   of all experiments are stored in compressed form in one column, even
   though individual initial states are estimated for each experiment.
 */
class Prg_SFunctionEst: public Prg_SFunction {

 protected:
  int		_mdl_np;	///< number of model parameters

  Omu_VariableVec _mdl_p;	///< model parameters (note: default min is 0)
  Omu_VariableVec _mdl_x0;	///< initial states

  IVECP		_mdl_p_active;	///< indicate estimated parameters
  IVECP		_mdl_x0_active; ///< indicate estimated states
  VECP		_mdl_der_x0_min;///< minimum for time derivative of x0
  VECP		_mdl_der_x0_max;///< maximum for time derivative of x0
  IVECP		_mdl_y_active; 	///< indicate measured outputs
  VECP 		_mdl_p_nominal;	///< nominal parameter values (for scaling)
  VECP 		_mdl_x_nominal; ///< nominal state values (for scaling)
  VECP		_mdl_y_nominal;	///< nominal output values (for scaling)

  int		_nx;	///< number of states for optimizer
  int		_np;	///< number of estimated parameters
  int		_nx0;	///< number of estimated initial states
  int		_ny;	///< number of reference outputs for estimation
  bool 		_multistage; 	///< treat as multistage problem

  int		_nex;	///< number of experiments used for estimation
  MATP 		_mdl_x0s;///< initial states for each experiment
  IVECP		_exs; 	///< index identifying experiment in estimation data

  // vectors for inputs, states, and outputs
  MATP	_mdl_us;	///< given model inputs
  MATP	_mdl_ys;	///< calculated model outputs
  MATP	_ys_ref;	///< reference values for active outputs

  // confidence things
  MATP 	_M;		///< measurement matrix M=dy/d(p,x0)
  MATP 	_dydx;		///< help variable dy/dx
  MATP 	_dxdpx0;	///< help variable dx/d(p,x0)
  MATP 	_dydpx0;	///< help variable dy/d(p,x0)
  MATP 	_dxfdx;		///< help variable dxf/dx
  MATP 	_dxfdpx0;	///< help variable dxf/d(p,x0)

  /**
   * @name Implementation of predefined methods.
   * @see Omu_Program
   */
  //@{
  void setup_stages(IVECP ks, VECP ts);

  void setup(int k,
	     Omu_VariableVec &x, Omu_VariableVec &u, Omu_VariableVec &c);

  void setup_struct(int k,
		    const Omu_VariableVec &x, const Omu_VariableVec &u,
		    Omu_DependentVec &xt, Omu_DependentVec &F,
		    Omu_DependentVec &f,
		    Omu_Dependent &f0, Omu_DependentVec &c);

  void init_simulation(int k,
		       Omu_VariableVec &x, Omu_VariableVec &u);

  void update(int kk,
	      const Omu_StateVec &x, const Omu_Vec &u,
	      const Omu_StateVec &xf,
	      Omu_DependentVec &f, Omu_Dependent &f0,
	      Omu_DependentVec &c);

  void consistic(int kk, double t,
		 const Omu_StateVec &x, const Omu_Vec &u,
		 Omu_DependentVec &xt);

  void continuous(int kk, double t,
		  const Omu_StateVec &x, const Omu_Vec &u,
		  const Omu_StateVec &xp, Omu_DependentVec &F);
  //@}

  /**
   * @name Helper methods for reading and writing S-function arguments.
   */
  //@{
  void read_mx_args(VECP p); ///< read _mx_args into vector p
  void write_mx_args(VECP p); ///< write vector p to _mx_args
  /// write active parameters that are packed in p to _mx_args
  void write_active_mx_args(VECP p);
  //@}

 public:

  Prg_SFunctionEst();		///< constructor
  ~Prg_SFunctionEst();		///< destructor

  char *name() {return "SFunctionEst";} ///< name SFunctionEst

  /**
   * @name Access methods for program specific members (If prefix: prg_)
   */
  //@{

  /// number of experiments used for estimation
  int nex() const {return _nex;}
  /// set number of experiments
  void set_nex(int val) {_nex = val;}

  /// indicate if problem is treated with one stage per time interval
  bool multistage() const {return _multistage;}
  /// set multistage flag
  void set_multistage(bool val) {_multistage = val;}

  /// reference values for active model outputs (size: KK+1 . ny)
  const MATP ys_ref() const {return _ys_ref;}
  /// set reference outputs
  void set_ys_ref(const MATP val) {m_copy_elements(val, _ys_ref);}

  /// measurement matrix M=dy/d(p,x0) (size: (KK+1)*ny . np+nx0)
  const MATP M() const {return _M;}

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

  /// nominal parameter values (for scaling)
  const VECP mdl_p_nominal() const {return _mdl_p_nominal;}

  /// model initial states
  /// (read only, note: write initial states of experiments via mdl_x0s)
  const VECP mdl_x0() const {return _mdl_x0;}

  /// lower bounds for model parameters (note: default is 0)
  const VECP mdl_x0_min() const {return _mdl_x0.min;}

  /// upper bounds for model parameters
  const VECP mdl_x0_max() const {return _mdl_x0.max;}

  /// indicate estimated initial states
  const IVECP mdl_x0_active() const {return _mdl_x0_active;}

  /// minimum for time derivative of x0  
  const VECP mdl_der_x0_min() const {return _mdl_der_x0_min;}

  /// maximum for time derivative of x0  
  const VECP mdl_der_x0_max() const {return _mdl_der_x0_max;}

  /// nominal state values (for scaling)
  const VECP mdl_x_nominal() const {return _mdl_x_nominal;}

  /// indicate measured outputs
  const IVECP mdl_y_active() const {return _mdl_y_active;}

  /// nominal output values (for scaling)
  const VECP mdl_y_nominal() const {return _mdl_y_nominal;}

  /// initial states for each experiment (size: nex . mdl_nx)
  const MATP mdl_x0s() const {return _mdl_x0s;}

  /// model inputs (size: KK+1 . mdl_nu)
  const MATP mdl_us() const {return _mdl_us;}

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

  /// set nominal parameters
  void set_mdl_p_nominal(const VECP v) {v_copy_elements(v, _mdl_p_nominal);}

  // model initial states are read only

  /// set lower bounds for initial states
  void set_mdl_x0_min(const VECP v) {v_copy_elements(v, _mdl_x0.min);}

  /// set upper bounds for initial states
  void set_mdl_x0_max(const VECP v) {v_copy_elements(v, _mdl_x0.max);}

  /// set estimated initial states
  void set_mdl_x0_active(const IVECP v) {iv_copy_elements(v, _mdl_x0_active);}

  /// set minimum for dx0dt
  void set_mdl_der_x0_min(const VECP v) {v_copy_elements(v, _mdl_der_x0_min);}

  /// set maximum for dx0dt
  void set_mdl_der_x0_max(const VECP v) {v_copy_elements(v, _mdl_der_x0_max);}

  /// set nominal states
  void set_mdl_x_nominal(const VECP v) {v_copy_elements(v, _mdl_x_nominal);}

  /// set measured outputs
  void set_mdl_y_active(const IVECP v) {iv_copy_elements(v, _mdl_y_active);}

  /// set nominal outputs
  void set_mdl_y_nominal(const VECP v) {v_copy_elements(v, _mdl_y_nominal);}

  /// set initial states
  void set_mdl_x0s(const MATP v) {m_copy_elements(v, _mdl_x0s);}

  /// set model inputs
  void set_mdl_us(const MATP v) {m_copy_elements(v, _mdl_us);}

  // model outputs are read only

  //@}
};  

#endif


/**
 * @file Prg_SFunctionOpt.h
 *    Optimal control problem for a model given as MEX S-function.
 *
 * rf, 7/25/00
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

#ifndef Prg_SFunctionOpt_H
#define Prg_SFunctionOpt_H

#include "Prg_SFunction.h"
#include "Omu_Variables.h"

/** Extend Omu_VariableVec with attributes for optimization criterion. */
class Omu_OptVarVec: public Omu_VariableVec {
public:
  VECP 	weight1;	///< weight for linear objective term (default: 0.0)
  VECP 	weight2;	///< weight for quadratic objective term (default: 0.0)
  VECP 	ref; 		///< reference value for quadratic term (default: 0.0)
  IVECP active; 	///< indicate used variables (default: 0 -- not used)

  Omu_OptVarVec(); 		///< allocate empty vectors
  virtual ~Omu_OptVarVec(); 	///< destroy vectors
  virtual void resize(int n); 	///< resize and initialize vectors
};

/**
   Optimal control problem for a model given as MEX S-function treated with
   multi-stage control vector parameterization.
   The optimization time horizon @f$[t_0,t_f]@f$ is split into @f$k=0,...,K@f$ 
   stages with time points @f$t_0=t^0<t^1<\ldots<t^K=t_f@f$. Each stage may 
   be further subdivided into @f$sps@f$ sample periods per stage. This leads 
   to the overall number of @f$KK=sps\,K@f$ sample periods with the sample
   time points @f$t^{kk}, kk=0,...,KK@f$ (sample time points within
   a stage are for instance useful to better treat path constraints).
   Sought control trajectories are described piecewise as constant or
   linear functions of control parameters in each stage.

   In the following all vector operations are defined element wise.
   The treated optimization problem reads
   @f[
   \begin{array}{l}
    \displaystyle{}
    J\ =\ \sum_{kk=0}^{KK} \sum_{i=1}^{n_u} \Delta t_{u,i}^{kk} \left\{
      u_{weight1}\,u(t^{kk})
      \ +\ u_{weight2}\left[u(t^{kk})-\frac{u_{ref}}{u_{nominal}}\right]^2
      \right\}_i
    \\[4ex] \displaystyle \qquad
    \ + \ \sum_{k=0}^{K-1} (t^{k+1}-t^{k}) \sum_{i=1}^{n_u} \left\{
      {der\_u}_{weight1}\,du^k
      \ +\ {der\_u}_{weight2}\left[du^k
                                -\frac{{der\_u}_{ref}}{u_{nominal}}\right]^2
      \right\}_i
    \\[4ex] \displaystyle \qquad
    \ + \ \sum_{kk=0}^{KK} \Delta t^{kk} \sum_{i=1}^{n_y} \left\{
      y_{weight1}\,y(t^{kk})
      \ +\ y_{weight2}\left[y(t^{kk})-\frac{y_{ref}}{y_{nominal}}\right]^2
    \right\}_i
    \\[4ex] \displaystyle \qquad
    \ + \ \sum_{kk=0}^{KK} \Delta t^{kk} \sum_{i=1}^{n_y} \left\{
         y_{soft\_weight1}\,s^{kk} + y_{soft\_weight2}\,s^{kk}s^{kk}
    \right\}_i
    \\[4ex] \displaystyle \qquad
    \ + \ \sum_{i=1}^{n_y} \left\{
      y_{0\_weight1}\,y(t_0)
      \ +\ y_{0\_weight2}\left[y(t_0)-\frac{y_{ref}}{y_{nominal}}\right]^2
    \right\}_i
    \\[4ex] \displaystyle \qquad
    \ + \ \sum_{i=1}^{n_y} \left\{
      y_{f\_weight1}\,y(t_f)
      \ +\ y_{f\_weight2}\left[y(t_f)-\frac{y_{ref}}{y_{nominal}}\right]^2
    \right\}_i
    \quad\to\quad \min
   \end{array}
   @f]
   with
   @f[
   \begin{array}{l}
    \displaystyle \Delta t_{u,i}^{kk} = \left\{\begin{array}{ll}
      t^{kk+1} - t^{kk}, & u_{order,i} = 0 \ \ \mbox{and}\ \ kk < KK, \\
      0, & u_{order,i} = 0 \ \ \mbox{and}\ \ kk = KK, \\
      \Delta t^{kk}, & u_{order,i} = 1,
    \end{array}\right. \\[5ex]
    \displaystyle \Delta t^{kk} = \frac{1}{2}\left\{\begin{array}{ll}
      t^{kk+1} - t^{kk}, & kk=0, \\
      t^{kk} - t^{kk-1}, & kk=KK, \\
      t^{kk+1} - t^{kk-1}, & \mbox{else},
    \end{array}\right.
   \end{array}
   @f]
   subject to the model given with the S-function methods
   mdlDerivatives @f$f@f$ and mdlOutputs @f$g,\ t\in[t_0,t_f]@f$
   @f[
   \begin{array}{l}
    \displaystyle \dot{x}(t) = 
     \frac{f[x_{nominal}\,x(t),\ u_{nominal}\,u(t)]}{x_{nominal}}, \\[3ex]
    \displaystyle y(t) = 
     \frac{g[x_{nominal}\,x(t),\ u_{nominal}\,u(t)]}{y_{nominal}}
      \ +\ \frac{y_{bias}}{y_{nominal}},
   \end{array}
   @f]
   with piecewise constant or linear approximation of @f$u(t)@f$ either
   using optimized control parameters @f$du^k@f$ or given inputs @f$us@f$ 
   @f[
   \begin{array}{ll}
    \left\{ u(t) = u(t^{k-1}) + (t^{k}-t^{k-1})du^{k} \right\}_i,
    & i \in \mbox{find}(u_{active}\ \mbox{and}\ u_{order}=0), \\[1ex]
    & t\in[t^{k},t^{k+1}),\ k=1,\ldots,K-1, \\[3ex]
    \left\{ \dot{u}(t) = du^{k} \right\}_i,
    & i \in \mbox{find}(u_{active}\ \mbox{and}\ u_{order}=1), \\[1ex]
    & t\in[t^{k},t^{k+1}),\ k=0,\ldots,K-1, \\[3ex]
    \left\{ u(t) = \displaystyle \frac{us^{kk}}{u_{nominal}} \right\}_i, &
    i \in \mbox{find}(\mbox{not}\ u_{active}\ \mbox{and}\ u_{order}=0),\\[3ex]
    \left\{ u(t) = \displaystyle
      \frac{t^{kk+1}-t}{t^{kk+1}-t^{kk}}\ \frac{us^{kk}}{u_{nominal}}
      + \frac{t-t^{kk}}{t^{kk+1}-t^{kk}} \ \frac{us^{kk+1}}{u_{nominal}}
    \right\}_i, &
    i \in \mbox{find}(\mbox{not}\ u_{active}\ \mbox{and}\ u_{order}=1), \\[3ex]
    & t\in[t^{kk},t^{kk+1}),\ kk=0,\ldots,KK-1,
   \end{array}
   @f]
   and subject to the constraints
   @f[
   \begin{array}{rcccll}
    \displaystyle && \{ x(t^0) &=& \displaystyle \frac{x^0}{x_{nominal}} \}_i, 
        \quad & i\notin\mbox{find}(x_{0\_active}), \\[3ex]
    \displaystyle \frac{y_{min}}{y_{nominal}} &<& y(t^{0})
        &<& \displaystyle \frac{y_{max}}{y_{nominal}}, \\[3ex]
    \displaystyle && u(t^0) &=& \displaystyle \frac{us^0}{u_{nominal}},
    \quad & nus_{fixed}>0, \\[3ex]
    \displaystyle \frac{u_{min}}{u_{nominal}} &<& u(t^{k})
        &<& \displaystyle \frac{u_{max}}{u_{nominal}}, \quad &
        k=0,\ldots,K \ \ \mbox{and}\ \ sps\,k \ge nus_{fixed}, \\[3ex]
    \displaystyle && du^k &=& \displaystyle du^k_{initial}, \quad &
        k = 0,\ldots,K-1 \ \ \mbox{and}\ \ sps\,k < nus_{fixed}-1, \\[3ex]
    \displaystyle \frac{{der\_u}_{min}}{u_{nominal}} &<& du^{k}
        &<& \displaystyle \frac{{der\_u}_{max}}{u_{nominal}}, &
        k=0,\ldots,K-1 \ \ \mbox{and}\ \ sps\,k \ge nus_{fixed}-1, \\[3ex]
    \displaystyle \frac{y_{min}}{y_{nominal}} &<& y(t^{kk})
        &<& \displaystyle \displaystyle \frac{y_{max}}{y_{nominal}}, \quad &
        kk=0,\ldots,KK, \\[3ex]
    \displaystyle \frac{y_{soft\_min}}{y_{nominal}} - s^{kk} &<& y(t^{kk})
        &<& \displaystyle \frac{y_{soft\_max}}{y_{nominal}} + s^{kk}, \\[3ex]
    \displaystyle && s^{kk} &>& 0, \quad & kk=0,\ldots,KK, \\[3ex]
    \displaystyle \frac{y_{f\_min}}{y_{nominal}} &<& y(t_f)
        &<& \displaystyle \displaystyle \frac{y_{f\_max}}{y_{nominal}}.
   \end{array}
   @f]
   The initial guess is taken from given initial states and model inputs
   @f[
   \begin{array}{rcll}
    x_{initial}(t^0) &=& \displaystyle \frac{x^0}{x_{nominal}}, \\[3ex]
    u_{initial}(t^0) &=& \displaystyle \frac{us^0}{u_{nominal}}, \\[3ex]
    du^k_{initial} &=& \displaystyle 
     \left\{\frac{us^{sps\,(k+1)} - us^{sps\,k}}
                 {(t^{sps\,(k+1)}-t^{sps\,k})\,u_{nominal}}
     \right\}_i, & i \in \mbox{find}(u_{active}), \\[1ex]
               &&& k=0,\ldots,K-1.
   \end{array}
   @f]

   The problem is treated as multistage problem with K stages per default. 
   Consequently additional K junction conditions (equality constraints)
   are introduced for the state variables x and the piecewise linear
   approximated control trajectories u. Alternatively the problem can
   be treated without stages applying pure control vector parameterization
   and hiding model states from the optimizer.

   Model inputs and outputs can be accessed through
   @f[
   \begin{array}{l}
    us^{kk} = u_{nominal}\,u(t^{kk}), \\[1ex]
    ys^{kk} = y_{nominal}\,y(t^{kk}), \quad kk=0,\ldots,KK.
   \end{array}
   @f]
 */
class Prg_SFunctionOpt: public Prg_SFunction {

 protected:
  IVECP 	_mdl_x0_active;	///< free initial states (default: 0)
  Omu_OptVarVec _mdl_y0; 	///< model outputs at initial time
  Omu_OptVarVec _mdl_u; 	///< model inputs
  Omu_OptVarVec _mdl_der_u; 	///< rates of change of inputs
  Omu_OptVarVec _mdl_y; 	///< model outputs
  Omu_OptVarVec _mdl_y_soft; 	///< attributes for relaxed output constraints
  Omu_OptVarVec _mdl_yf; 	///< model outputs at final time

  IVECP 	_mdl_u_order; 	///< interpolation order (default: 1 (linear))

  double 	_t_nominal; 	///< nominal time (used internally for scaling)
  VECP 		_mdl_u_nominal;	///< nominal inputs (for scaling)
  VECP 		_mdl_x_nominal;	///< nominal states (for scaling)
  VECP 		_mdl_y_nominal;	///< nominal outputs (for scaling)

  VECP 		_mdl_y_bias;	///< bias correction (offset) for outputs 

  int		_nx;	///< number of states for optimizer
  int		_nu;	///< number of optimized control inputs
  int		_nc;	///< number of constrained outputs
  int		_nc0;	///< number of constrained/used outputs at initial time
  int		_ncf;	///< number of constrained/used outputs at final time
  int		_ns;	///< number of slack variables for soft constraints
  int		_nsc;	///< number of soft constraints
  bool 		_multistage; 	///< treat as multistage problem

  /**
   * Number of sample periods per stage (default: 1).
   * The value can be increased to devide each control interval into
   * multiple sample periods, e.g. for evaluating constraints and
   * the objective within control intervals. Currently _sps>1 is only
   * supported for multistage problems.
   */
  int 		_sps;

  MATP	_mdl_us;	///< given model inputs (controls and disturbances)
  MATP	_mdl_ys;	///< calculated model outputs

  /**
   * Number of fixed control inputs at begin of time horizon (default: 0).
   * The initial value is fixed for _nus_fixed=1, the initial and the 
   * second value are fixed for _nus_fixed=2, and so on).
   */
  int 	_nus_fixed;

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
   * @name Overloaded gradient routines
   * These routines call the S-function method mdlJacobian if available;
   * otherwise they direct calls to the according methods by Omu_Program.
   */
  //@{

  /**
   * Overloaded update routine for obtaining gradients.
   */
  void update_grds(int kk, 
		   const Omu_StateVec &x, const Omu_Vec &u,
		   const Omu_StateVec &xf,
		   Omu_DependentVec &f, Omu_Dependent &f0,
		   Omu_DependentVec  &c);

  /**
   * Overloaded continuous routine for obtaining gradients.
   */
  void continuous_grds(int kk, double t,
		       const Omu_StateVec &x, const Omu_Vec &u,
		       const Omu_StateVec &xp, Omu_DependentVec &F);
  //@}

 public:

  Prg_SFunctionOpt();		///< constructor
  ~Prg_SFunctionOpt();		///< destructor

  char *name() {return "SFunctionOpt";} ///< name SFunctionOpt

  /**
   * @name Access methods for program specific members (If prefix: prg_)
   */
  //@{
  /// Numbers of fixed control inputs at begin of time horizon (default: 0)
  int nus_fixed() const {return _nus_fixed;}
  /// set numbers of fixed constrol inputs
  void set_nus_fixed(int val) {_nus_fixed = val;}

  /// indicate if problem is treated with one stage per time interval
  bool multistage() const {return _multistage;}
  /// set multistage flag
  void set_multistage(bool val) {_multistage = val;}

  //@}

  /**
   * @name Read methods for model specific members (no If prefix).
   * Note that the prefix mdl_ is omitted in the detailed mathematical
   * problem description.
   */
  //@{
  ///< free initial states (default: 0)
  const IVECP mdl_x0_active() const {return _mdl_x0_active;}

  /// lower bounds for model outputs at initial time
  const VECP mdl_y0_min() const {return _mdl_y0.min;}

  /// upper bounds for model outputs at initial time
  const VECP mdl_y0_max() const {return _mdl_y0.max;}

  /// weight for linear objective term at initial time (default: 0)
  const VECP mdl_y0_weight1() const {return _mdl_y0.weight1;}

  /// weight for quadratic objective term at initial time (default: 0)
  const VECP mdl_y0_weight2() const {return _mdl_y0.weight2;}

  ///< interpolation order (0 (constant) or 1 (linear), default: 1)
  const IVECP mdl_u_order() const {return _mdl_u_order;}

  /// indicate optimized inputs
  const IVECP mdl_u_active() const {return _mdl_u.active;}

  /// nominal input values (for scaling)
  const VECP mdl_u_nominal() const {return _mdl_u_nominal;}

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

  /// bias for outputs
  const VECP mdl_y_bias() const {return _mdl_y_bias;}

  /// nominal output values (for scaling)
  const VECP mdl_y_nominal() const {return _mdl_y_nominal;}

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

  /// lower bounds for model outputs at final time
  const VECP mdl_yf_min() const {return _mdl_yf.min;}

  /// upper bounds for model outputs at final time
  const VECP mdl_yf_max() const {return _mdl_yf.max;}

  /// weight for linear objective term (default: 0)
  const VECP mdl_yf_weight1() const {return _mdl_yf.weight1;}

  /// weight for quadratic objective term (default: 0)
  const VECP mdl_yf_weight2() const {return _mdl_yf.weight2;}

  /// nominal state values (for scaling)
  const VECP mdl_x_nominal() const {return _mdl_x_nominal;}

  /// model inputs (size: KK+1 . mdl_nu)
  const MATP mdl_us() const {return _mdl_us;}

  /// model outputs (read only, size: KK+1 . mdl_ny)
  const MATP mdl_ys() const {return _mdl_ys;}

  //@}
  /**
   * @name Write methods for model specific members (no If prefix).
   */
  //@{
  /// set free inputs
  void set_mdl_x0_active(const IVECP v) {iv_copy_elements(v, _mdl_x0_active);}

  /// set lower bounds for final model outputs
  void set_mdl_y0_min(const VECP v) {v_copy_elements(v, _mdl_y0.min);}

  /// set upper bounds for final model outputs
  void set_mdl_y0_max(const VECP v) {v_copy_elements(v, _mdl_y0.max);}

  /// set linear weight
  void set_mdl_y0_weight1(const VECP v) {v_copy_elements(v, _mdl_y0.weight1);}

  /// set quadratic weight
  void set_mdl_y0_weight2(const VECP v) {v_copy_elements(v, _mdl_y0.weight2);}

  ///< set interpolation order
  void set_mdl_u_order(const IVECP v) {iv_copy_elements(v, _mdl_u_order);}

  /// set optimized inputs
  void set_mdl_u_active(const IVECP v) {iv_copy_elements(v, _mdl_u.active);}

  /// set nominal inputs
  void set_mdl_u_nominal(const VECP v) {v_copy_elements(v, _mdl_u_nominal);}

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

  /// set output bias
  void set_mdl_y_bias(const VECP v) {v_copy_elements(v, _mdl_y_bias);}

  /// set nominal outputs
  void set_mdl_y_nominal(const VECP v) {v_copy_elements(v, _mdl_y_nominal);}

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

  /// set lower bounds for final model outputs
  void set_mdl_yf_min(const VECP v) {v_copy_elements(v, _mdl_yf.min);}

  /// set upper bounds for final model outputs
  void set_mdl_yf_max(const VECP v) {v_copy_elements(v, _mdl_yf.max);}

  /// set linear weight
  void set_mdl_yf_weight1(const VECP v) {v_copy_elements(v, _mdl_yf.weight1);}

  /// set quadratic weight
  void set_mdl_yf_weight2(const VECP v) {v_copy_elements(v, _mdl_yf.weight2);}

  /// set nominal states
  void set_mdl_x_nominal(const VECP v) {v_copy_elements(v, _mdl_x_nominal);}

  /// set model inputs
  void set_mdl_us(const MATP v) {m_copy_elements(v, _mdl_us);}

  /// model outputs are read only

  //@}
};  

#endif


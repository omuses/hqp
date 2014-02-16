/**
 * @file Prg_SFunctionOpt.h
 *    Optimal control problem for a model given as MEX S-function.
 *
 * rf, 7/25/00
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
   The model time may be parameterized by defining a piecewise constant
   model input as scaling factor.
   Sought control trajectories are described piecewise as constant or
   linear functions of control parameters in each stage.

   In the following all vector operations are defined element wise.
   The treated optimization problem reads
   @f[
   \begin{array}{l}
    \displaystyle 
    J\ =\ \sum_{i=1}^{n_y} \left\{
      y_{0\_weight1}\left[y(t_0)-\frac{y_{ref}}{y_{nominal}}\right]
      \ +\ y_{0\_weight2}\left[y(t_0)-\frac{y_{ref}}{y_{nominal}}\right]^2
    \right\}_i
    \\[4ex] \displaystyle \qquad 
    \ + \ \sum_{kk=0}^{KK} \sum_{i=1}^{n_u} \Delta t_{u,i}^{kk} \left\{
      u_{weight1}\left[u(t^{kk})-\frac{u_{ref}}{u_{nominal}}\right]
      \ +\ u_{weight2}\left[u(t^{kk})-\frac{u_{ref}}{u_{nominal}}\right]^2
      \right\}_i
    \\[4ex] \displaystyle \qquad
    \ + \ \sum_{k=0}^{K-1} (t^{k+1}-t^{k}) \sum_{i=1}^{n_u} \left\{
      {der\_u}_{weight1}\left[du^k
                                -\frac{{der\_u}_{ref}}{u_{nominal}}\right]
      \ +\ {der\_u}_{weight2}\left[du^k
                                -\frac{{der\_u}_{ref}}{u_{nominal}}\right]^2
      \right\}_i
    \\[4ex] \displaystyle \qquad
    \ + \ \sum_{kk=0}^{KK} \Delta t_{y,i}^{kk} \sum_{i=1}^{n_y} \left\{
      y_{weight1}\left[y(t^{kk})-\frac{y_{ref}}{y_{nominal}}\right]
      \ +\ y_{weight2}\left[y(t^{kk})-\frac{y_{ref}}{y_{nominal}}\right]^2
    \right\}_i
    \\[4ex] \displaystyle \qquad
    \ + \ \sum_{i=1}^{n_y} \left\{
      y_{f\_weight1}\left[y(t_f)-\frac{y_{ref}}{y_{nominal}}\right]
      \ +\ y_{f\_weight2}\left[y(t_f)-\frac{y_{ref}}{y_{nominal}}\right]^2
    \right\}_i
    \\[4ex] \displaystyle \qquad
    \ + \ \sum_{kk=0}^{KK} \Delta t_{y,i}^{kk} \sum_{i=1}^{n_y} \left\{
      y_{soft\_weight1}\,s^{kk} + y_{soft\_weight2}\,s^{kk}s^{kk}
    \right\}_i
    \\[4ex] \displaystyle \qquad
    \ + \ \sum_{i=1}^{n_y} \left\{
      y_{f\_soft\_weight1}\,s_f + y_{f\_soft\_weight2}\,s_fs_f
    \right\}_i
    \\[4ex] \displaystyle \qquad
    \ \to\quad \min
   \end{array}
   @f]
   with
   @f[
   \begin{array}{l}
    \displaystyle \Delta t_{u,i}^{kk} = 
    \left\{\begin{array}{ll}
      \Delta t_0^{kk}, & u_{order,i} = 0, \\[1ex]
      \Delta t^{kk}, & u_{order,i} = 1,
    \end{array}\right. \\[5ex]
    \displaystyle \Delta t_{y,i}^{kk} = 
    \left\{\begin{array}{ll}
      \Delta t_0^{kk}, & y_{order,i} = 0, \\[1ex]
      \Delta t^{kk}, & y_{order,i} = 1,
    \end{array}\right. \\[5ex]
    \displaystyle \Delta t_0^{kk} = 
    \left\{\begin{array}{ll}
      t_{scale}^{kk+1}(t^{kk+1} - t^{kk}), & kk < KK, \\[1ex]
      0, & kk = KK,
    \end{array}\right. \\[5ex]
    \displaystyle \Delta t^{kk} = \frac{1}{2}
    \left\{\begin{array}{ll}
      t_{scale}^{kk+1}(t^{kk+1} - t^{kk}), & kk=0, \\[1ex]
      t_{scale}^{kk}(t^{kk} - t^{kk-1}), & kk=KK, \\[1ex]
      t_{scale}^{kk+1}(t^{kk+1} - t^{kk}) + t_{scale}^{kk}(t^{kk} - t^{kk-1}), & \mbox{else},
    \end{array}\right. \\[5ex]
    \displaystyle t_{scale} = 
    \left\{\begin{array}{ll}
      u_{t\_scale\_idx}, & t\_scale\_idx \ge 0, \\
      1, & \mbox{else},
    \end{array}\right.
   \end{array}
   @f]
   subject to the model given with the S-function methods
   mdlDerivatives @f$f@f$, mdlOutputs @f$g@f$ and mdlUpdate
   @f$h,\ t\in[t_0,t_f]@f$, as well as with parameterized time @f$\tau@f$
   @f[
   \begin{array}{l}
    \displaystyle \quad\ \frac{d\tau(t)}{dt} = t_{scale}(t), \\[3ex]
    \displaystyle x(\tau(t)) =
    \left[\begin{array}{ll}
      x_d^{kk}, & t\in[t^{kk},t^{kk+1}),\ kk=0,\ldots,KK \\ x_c(\tau(t))
    \end{array}\right], \\[3ex]
    \displaystyle x_d^{kk} =
     \frac{h[x_{nominal}\,x(\tau(t^-)),\ u_{nominal}\,u(\tau(t^-))]}
          {x_{nominal_d}},\quad t=t^{kk},\ kk=1,\ldots,KK,\\[3ex]
    \displaystyle \frac{dx_c(\tau(t))}{d\tau(t)} =
     \frac{f[x_{nominal}\,x(\tau(t)),\ u_{nominal}\,u(\tau(t))]}
          {x_{nominal_c}}, \\[3ex]
    \displaystyle y(\tau(t)) = 
     \frac{g[x_{nominal}\,x(\tau(t)),\ u_{nominal}\,u(\tau(t))]}{y_{nominal}}
      \ +\ \frac{y_{bias}}{y_{nominal}},
   \end{array}
   @f]
   with piecewise constant or linear approximation of @f$u(t)@f$ either
   using optimized control parameters @f$du^k@f$ or given inputs @f$us@f$ 
   @f[
   \begin{array}{ll}
    \left\{ u(t) = u(t^{k-1}) + \Delta u(t^{k-1}) \right\}_i,
    & i \in \mbox{find}(u_{active}\ \mbox{and}\ u_{order}=0), \\[1ex]
    & t\in[t^{k},t^{k+1}),\ k=1,\ldots,K-1, \\[1ex]
    \ \ \Delta u(t^{k-1}) =
      \left\{\begin{array}{ll}
        (t^{k}-t^{k-1})du^{k-1}, & i = t\_scale\_idx \\[1ex]
        (t^{k}-t^{k-1})du^{k-1}t_{scale}^k, & \mbox{else}
      \end{array}\right.
    & k=1,\ldots,K-1, \\[3ex]
    \left\{ \frac{du(\tau(t))}{d\tau(t)} = du^{k} \right\}_i,
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
   and subject to the constraints at initial and final time
   @f[
   \begin{array}{rcccll}
    \displaystyle && \{ x(t^0) &=& \displaystyle \frac{x^0}{x_{nominal}}\}_i,
        \quad & i\notin\mbox{find}(x_{0\_active}), \\[3ex]
    \displaystyle \left\{ \frac{x^0_{min}}{x_{nominal}} \right.
        &\le& x(t^{0}) &\le& 
        \displaystyle \left. \frac{x^0_{max}}{x_{nominal}} \right\}_i,
	\quad & i\in\mbox{find}(x_{0\_active}), \\[3ex]
    \displaystyle \left\{ \frac{der\_x^0_{min}}{x_{nominal}} \right.
        &\le& \dot{x}(t^{0}) &\le& 
        \displaystyle \left. \frac{der\_x^0_{max}}{x_{nominal}} \right\}_i,
	\quad & i\in\mbox{find}(x_{0\_active}), \\[3ex]
    \displaystyle \left\{ \frac{u^0_{min}}{u_{nominal}} \right.
        &\le& u(t^{0}) &\le& 
        \displaystyle \left. \frac{u^0_{max}}{u_{nominal}} \right\}_i,
	\quad & i\in\mbox{find}(u_{active}\ \mbox{and}\ u_{0,nfixed}=0), \\[3ex]
    \displaystyle && \{ u(t^0) &=& \displaystyle \frac{us^0}{u_{nominal}} \}_i,
        \quad & i \in \mbox{find}(u_{0,nfixed}>1-u_{order}), \\[3ex]
    \displaystyle\left\{\frac{us^0+der\_u_{min}}{u_{nominal}}\right. &\le& u(t^0)
        &\le& \displaystyle\left.\frac{us^0+der\_u_{max}}{u_{nominal}}\right\}_i,
        \ & i \in \mbox{find}(u_{0,nfixed}=1\ \mbox{and}\ u_{order}=0), \\[3ex]
    \displaystyle \frac{y_{0,min}}{y_{nominal}} &\le& y(t^{0})
        &\le& \displaystyle \frac{y_{0,max}}{y_{nominal}}, \\[3ex]
    \displaystyle \frac{u_{f\_min}}{u_{nominal}} &\le& u(t_f)
        &\le& \displaystyle \frac{u_{f\_max}}{u_{nominal}}, \\[3ex]
    \displaystyle \frac{y_{f\_min}}{y_{nominal}} &\le& y(t_f)
        &\le& \displaystyle \frac{y_{f\_max}}{y_{nominal}}, \\[3ex]
    \displaystyle \frac{y_{f\_soft\_min}}{y_{nominal}} - s_f &\le& y(t_f)
        &\le& \displaystyle \frac{y_{f\_soft\_max}}{y_{nominal}} + s_f, \\[3ex]
    \displaystyle && s_f &\ge& 0,
   \end{array}
   @f]
   as well as at all time points
   @f[
   \begin{array}{rcccll}
    \displaystyle\left\{ \frac{u_{min}}{u_{nominal}} \right. &\le& u(t^{k})
        &\le& \displaystyle\left. \frac{u_{max}}{u_{nominal}} \right\}_i, 
        \quad & i \in \mbox{find}(k \ge u_{0,nfixed} + u_{order} - 1), \\[1ex]
	&& && & k=0,\ldots,K, \\[2ex]
    \displaystyle \frac{{der\_u}_{min}}{u_{nominal}} &\le& du^{k}
        &\le& \displaystyle \frac{{der\_u}_{max}}{u_{nominal}}, \quad &
	k=0,\ldots,K-1, \\[3ex]
    \displaystyle \frac{x_{min}}{x_{nominal}} &\le& x(t^{k})
        &\le& \displaystyle \frac{x_{max}}{x_{nominal}}, \quad &
        k=0,\ldots,K, \\[3ex]
    \displaystyle \frac{y_{min}}{y_{nominal}} &\le& y(t^{kk})
        &\le& \displaystyle \displaystyle \frac{y_{max}}{y_{nominal}}, \quad &
        kk=0,\ldots,KK, \\[3ex]
    \displaystyle \frac{y_{soft\_min}}{y_{nominal}} - s^{kk} &\le& y(t^{kk})
        &\le& \displaystyle \frac{y_{soft\_max}}{y_{nominal}} + s^{kk}, \\[3ex]
    \displaystyle && s^{kk} &\ge& 0, \quad & kk=0,\ldots,KK.
   \end{array}
   @f]
   Note that path constraints over continuous time intervals can be 
   approximated by specifying @f$sps>1@f$, leading to output constraints 
   at interior sample time points.

   The actually optimized rates of changes @f$du^k, k=0,\ldots,K-1@f$ for
   active inputs may be set to fixed values by specifying @f$u_{0,nfixed}@f$
   and @f$u_{decimation}@f$ (defaults: 1), fixing initial values and holding
   optimized inputs constant over multiple stages, respectively
   @f[
    du^k_i = \left\{\begin{array}{ll}
      \displaystyle 0 , \quad &
        i \in \mbox{find}(\mbox{mod}(k+1, u_{decimation}) \ne 0), \\[2ex]
      \displaystyle du^k_{initial,i}, \quad &
        i \in \mbox{find}(\mbox{mod}(k+1, u_{decimation}) = 0\ \mbox{and}\ 
	k < u_{0,nfixed}+u_{order}-2), \\[2ex]
      \displaystyle \mbox{free} , \quad & \mbox{else}.
      \end{array}\right.
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
   are introduced for the state variables x and the control trajectories u. 
   Alternatively the problem can be treated without stages applying pure 
   control vector parameterization and hiding model states from the optimizer.

   When treated as multistage problem, the additional optimization
   variables introduced for states are normally initialized with the results
   of an initial-value simulation for the given initial states and inputs 
   (multistage=1). Alternatively, with multistage=2, the multistage problem
   may be treated with multiple shooting, i.e. all states are initialized
   with the initial states @f$x^0@f$
   @f[
   \begin{array}{rcll}
    x_{initial}(t^k) &=& \displaystyle \frac{x^0}{x_{nominal}},
      & k=1,\ldots,K
   \end{array}
   @f]
   or with an explicitly given initial guess for states @f$xs@f$ in all stages
   @f[
   \begin{array}{rcll}
    x_{initial}(t^k) &=& \displaystyle \frac{xs^{sps\,k}}{x_{nominal}},
      & k=1,\ldots,K.
   \end{array}
   @f]
   The multiple shooting method is advantageous if the expected state
   trajectories are known (e.g. constant at a set point or following a
   ramp), but if an initial guess for the control trajectories causing
   those state trajectories is unknown.

   Model inputs, states and outputs can be accessed through
   @f[
   \begin{array}{l}
    us^{kk} = u_{nominal}\,u(t^{kk}), \\[1ex]
    xs^{kk} = x_{nominal}\,x(t^{kk}), \\[1ex]
    ys^{kk} = y_{nominal}\,y(t^{kk}), \quad kk=0,\ldots,KK.
   \end{array}
   @f]
 */
class Prg_SFunctionOpt: public Prg_SFunction {

 protected:
  Omu_VariableVec _mdl_x0;	///< initial states for optimization
  IVECP 	_mdl_x0_active;	///< free initial states (default: 0)
  VECP		_mdl_der_x0_min;///< minimum for time derivative of x0
  VECP		_mdl_der_x0_max;///< maximum for time derivative of x0
  Omu_VariableVec _mdl_u0;	///< initial inputs
  Omu_OptVarVec _mdl_y0; 	///< model outputs at initial time
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

  VECP		_taus;		///< scaled time communicated to outside as prg_ts
  int 		_t_scale_idx; 	///< index into mdl_u vector for variable used for scaling of time
  int 		_t_active; 	///< time is being scaled
  int 		_t_scale_i; 	///< index of internal optimization variable for scaling of time
  double 	_t_scale_nominal;///< nominal value for time scaling (default: 1)
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
  int		_nsu;	///< number of slack variables for soft cns on inputs
  int		_nsuc;	///< number of soft constraints on inputs
  int		_nsf;	///< number of slacks for soft constraints at final time
  int		_nscf;	///< number of soft constraints at final time
  int 		_multistage; 	///< treat as multistage problem

  /**
   * Number of sample periods per stage (default: 1).
   * The value can be increased to devide each control interval into
   * multiple sample periods, e.g. for evaluating constraints and
   * the objective within control intervals. Currently _sps>1 is only
   * supported for multistage problems.
   */
  int 		_sps;

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

  /**
   * @name Implementation of predefined methods.
   * @see Omu_Program
   */
  //@{
  void setup_model();

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
		  const Omu_StateVec &dx, Omu_DependentVec &F);
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
		       const Omu_StateVec &dx, Omu_DependentVec &F);
  //@}

 public:

  Prg_SFunctionOpt();		///< constructor
  ~Prg_SFunctionOpt();		///< destructor

  const char *name() {return "SFunctionOpt";} ///< name SFunctionOpt

  /**
   * @name Access methods for program specific members (If prefix: prg_)
   */
  //@{
  /// indicate if problem is treated with one stage per time interval
  int multistage() const {return _multistage;}
  /// set multistage flag
  void set_multistage(int val) {_multistage = val;}

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
    Omu_Program::set_ts(_taus);
    // finally take over n_taus
    if (_t_scale_idx >= 0)
      v_copy_elements(n_taus, _taus);
  }	
  //@}

  /**
   * @name Read methods for model specific members (no If prefix).
   * Note that the prefix mdl_ is omitted in the detailed mathematical
   * problem description.
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

  /// minimum for time derivative of x0  
  const VECP mdl_der_x0_min() const {return _mdl_der_x0_min;}

  /// maximum for time derivative of x0  
  const VECP mdl_der_x0_max() const {return _mdl_der_x0_max;}

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

  /// nominal input values (for scaling)
  const VECP mdl_u_nominal() const {return _mdl_u_nominal;}

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

  /// nominal state values (for scaling)
  const VECP mdl_x_nominal() const {return _mdl_x_nominal;}

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
    for (int kk = 0; kk <= _KK; kk++)
      for (int i = 0; i < _mdl_nx; i++)
	_mdl_xs[kk][i] = v[i];
  }

  /// set free initial states
  void set_mdl_x0_active(const IVECP v) {iv_copy_elements(v, _mdl_x0_active);}

  /// set lower bounds for initial states
  void set_mdl_x0_min(const VECP v) {v_copy_elements(v, _mdl_x0.min);}

  /// set upper bounds for initial states
  void set_mdl_x0_max(const VECP v) {v_copy_elements(v, _mdl_x0.max);}

  /// set minimum for dx0dt
  void set_mdl_der_x0_min(const VECP v) {v_copy_elements(v, _mdl_der_x0_min);}

  /// set maximum for dx0dt
  void set_mdl_der_x0_max(const VECP v) {v_copy_elements(v, _mdl_der_x0_max);}

  /// set initial values for model inputs and copy them to all mdl_us
  void set_mdl_u0(const VECP v) {
    v_copy_elements(v, _mdl_u0);
    for (int kk = 0; kk <= _KK; kk++)
      for (int i = 0; i < _mdl_nu; i++)
	_mdl_us[kk][i] = v[i];
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

  /// set nominal inputs
  void set_mdl_u_nominal(const VECP v) {v_copy_elements(v, _mdl_u_nominal);}

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

  /// set nominal states
  void set_mdl_x_nominal(const VECP v) {v_copy_elements(v, _mdl_x_nominal);}

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


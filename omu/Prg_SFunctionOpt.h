/**
 * @file Prg_SFunctionOpt.h
 *    Optimal control problem for a model given as MEX S-function.
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

#ifndef Prg_SFunctionOpt_H
#define Prg_SFunctionOpt_H

#include "Prg_SFunction.h"

/**
   Optimal control problem for a model given as MEX S-function.
   The optimization time horizon @f$[t_0,t_f]@f$ is split into @f$k=0,...,K@f$ 
   stages with time points @f$t_0=t^0<t^1<\ldots<t^K=t_f@f$. Each stage may 
   be further subdivided. This leads to @f$KK=sps\ K@f$ sample periods, 
   where sps is the number of samples periods per stage, and with the sample
   time points @f$t^{kk}, kk=0,...,KK@f$. Additional sample time points within
   a stage are for instance useful for better treating path constraints.

   In the following all vector operations are defined element wise.
   The treated optimization problem reads
   @f[
   \begin{array}{l}
    \displaystyle{}
    J\ =\ \sum_{kk=0}^{KK} \Delta t^{kk} \sum_{i=1}^{n_u} \left\{
      u_{weight2}\left[u(t^{kk})-\frac{u_{ref}}{u_{nominal}}\right]^2
      \right\}_i
    \\[4ex] \displaystyle \qquad
    \ + \ \sum_{kk=0}^{KK-1} \Delta t^{kk} \sum_{i=1}^{n_u} \left\{
      {der\_u}_{weight2}\left[\dot{u}(t^{kk})\right]^2
      \right\}_i
    \\[4ex] \displaystyle \qquad
    \ + \ \sum_{kk=0}^{KK} \Delta t^{kk} \sum_{i=1}^{n_y} \left\{
      y_{weight2}\left[y(t^{kk})-\frac{y_{ref}}{y_{nominal}}\right]^2
    \right\}_i
    \\[4ex] \displaystyle \qquad
    \ + \ \sum_{kk=0}^{KK} \Delta t^{kk} \sum_{i=1}^{n_y} \left\{
         y_{soft\_weight1}\,s^{kk} + y_{soft\_weight2}\,s^{kk}s^{kk}
        \right\}_i
      \quad\to\quad \min
   \end{array}
   @f]
   with
   @f[
    \Delta t^{kk} = \frac{1}{2}\left\{\begin{array}{ll}
      t^{kk+1} - t^{kk}, & kk=0, \\
      t^{kk} - t^{kk-1}, & kk=KK, \\
      t^{kk+1} - t^{kk-1}, & \mbox{else},
    \end{array}\right.
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
   with piecewise linear approximation of @f$u(t)@f$ 
   either using optimized control parameters @f$u^k@f$ or given inputs 
   @f$us@f$ 
   @f[
   \begin{array}{ll}
    \left\{ \dot{u}(t) = u^{k}
    \right\}_i, & i \in \mbox{find}(u_{active}), \\[1ex]
    & t\in[t^{k},t^{k+1}),\ k=0,\ldots,K-1, \\[3ex]
    \left\{ u(t) = \displaystyle 
        \frac{t^{kk+1}-t}{t^{kk+1}-t^{kk}}\ \frac{us^{kk}}{u_{nominal}}
      + \frac{t-t^{kk}}{t^{kk+1}-t^{kk}} \ \frac{us^{kk+1}}{u_{nominal}}
    \right\}_i, & i \notin \mbox{find}(u_{active}), \\[1ex]
    & t\in[t^{kk},t^{kk+1}),\ kk=0,\ldots,KK-1,
   \end{array}
   @f]
   and subject to the constraints
   @f[
   \begin{array}{rcccll}
    \displaystyle &&x(t^0) &=& \displaystyle \frac{x^0}{x_{nominal}}, \\[3ex]
    \displaystyle &&u(t^0) &=& \displaystyle \frac{us^0}{u_{nominal}}, \\[3ex]
    \displaystyle \frac{u_{min}}{u_{nominal}} &<& u(t^{k})
        &<& \displaystyle \frac{u_{max}}{u_{nominal}}, \quad &
        k=1,\ldots,K, \\[3ex]
    \displaystyle \frac{{der\_u}_{min}}{u_{nominal}} &<& u^{k}
        &<& \displaystyle \frac{{der\_u}_{max}}{u_{nominal}}, \quad &
        k=0,\ldots,K-1, \\[3ex]
    \displaystyle \frac{y_{min}}{y_{nominal}} &<& y(t^{kk})
        &<& \displaystyle \displaystyle \frac{y_{max}}{y_{nominal}}, \quad &
        kk=0,\ldots,KK, \\[3ex]
    \displaystyle \frac{y_{min\_soft}}{y_{nominal}} - s^{kk} &<& y(t^{kk})
        &<& \displaystyle \frac{y_{max\_soft}}{y_{nominal}} + s^{kk}, 
        \ \ s^{kk} > 0, \quad & kk=0,\ldots,KK.
   \end{array}
   @f]
   The problem is treated as multi-stage problem with K stages. 
   Consequently additional K junction conditions (equality constraints)
   are introduced for the state variables x and the piecewise linear
   approximated control trajectories u.

   Model inputs and outputs can be accessed through
   @f[
   \begin{array}{l}
    us^{kk} = u_{nominal}\,u(t^{kk}), \quad kk=1,\ldots,KK, \\[1ex]
    ys^{kk} = y_{nominal}\,y(t^{kk}), \quad kk=0,\ldots,KK.
   \end{array}
   @f]
 */
class Prg_SFunctionOpt: public Prg_SFunction {

 protected:
  VECP 		_mdl_x_nominal;	///< nominal initial states (for scaling)

  IVECP		_mdl_u_active; 	///< indicate optimized model inputs
  VECP 		_mdl_u_nominal;	///< nominal inputs (for scaling)
  VECP 		_mdl_u_min;	///< lower bounds for optimized inputs
  VECP 		_mdl_u_max;	///< upper bounds for optimized inputs
  VECP 		_mdl_u_ref;	///< reference to be reached for control inputs
  VECP 		_mdl_u_weight2;	///< weight for quadratic objective
  VECP 		_mdl_der_u_min;	///< lower bounds for rate of change
  VECP 		_mdl_der_u_max;	///< upper bounds for rate of change
  VECP 		_mdl_der_u_weight2;///< weight for rate of change of controls

  /// indicate outputs used in criterion or hard constraints
  IVECP		_mdl_y_active;
  VECP 		_mdl_y_nominal;	///< nominal outputs (for scaling)
  VECP 		_mdl_y_bias;	///< bias correction (offset) for outputs 
  VECP 		_mdl_y_min;	///< lower bounds for outputs
  VECP 		_mdl_y_max;	///< upper bounds for outputs
  VECP 		_mdl_y_ref;	///< reference to be reached for active outputs
  VECP 		_mdl_y_weight2;	///< weight for quadratic objective terms
  VECP 		_mdl_y_min_soft;      ///< lower soft bound
  VECP 		_mdl_y_max_soft;      ///< upper soft bound
  VECP 		_mdl_y_soft_weight1;  ///< linear weight for bound violation
  VECP 		_mdl_y_soft_weight2;  ///< quadratic weight for bound violation

  int		_nx;	///< number of states for optimizer
  int		_nu;	///< number of optimized control inputs
  int		_nc;	///< number of constrained outputs
  int		_ns;	///< number of slack variables for soft constraints
  int		_nsc;	///< number of soft constraints

  int 		_sps;	///< number of sample periods per stage

  // vectors for inputs, states, and outputs
  MATP	_mdl_us;	///< given model inputs
  MATP	_mdl_ys;	///< calculated model outputs

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

 public:

  Prg_SFunctionOpt();		///< constructor
  ~Prg_SFunctionOpt();		///< destructor

  char *name() {return "SFunctionOpt";} ///< name SFunctionOpt
};  

#endif


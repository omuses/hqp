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

/**
   Estimation of initial states and parameters for a model given
   as MEX S-function. The estimation problem is formulated 
   at the discrete time points @f$t^0<t^1<\ldots<t^K@f$.
   In the following all vector operations are defined element wise.
   The treated estimation problem reads
   @f[
   \begin{array}{l}
    \displaystyle{}
    J\ =\ \sum_{k=0}^{K} \sum_{i=1}^{dim(I)} \left\{
      \left[\frac{ys^k(I)-ys^k_{ref}}{y_{nominal}(I)}\right]^2
      \right\}_i
      \quad\to\quad \min, \quad I = \mbox{find}(y_{active}),
   \end{array}
   @f]
   subject to the model given with the S-function methods
   mdlDerivatives @f$f@f$ and mdlOutputs  @f$g,\ t\in[t_0,t_f]@f$
   @f[
   \begin{array}{l}
    \displaystyle \dot{x}(t) = 
     \frac{f[p,\ x_{nominal}\,x(t),\ u(t)]}
          {x_{nominal}}, \\[3ex]
    \displaystyle y(t) = 
     \frac{g[p,\ x_{nominal}\,x(t),\ u(t)]}
          {y_{nominal}}, \\[3ex]
    ys^{k} = y_{nominal}y(t^{k}), \quad k=0,\ldots,K,
   \end{array}
   @f]
   with piecewise linear interpolation of the inputs
   @f[
   \begin{array}{l}
      u(t) = \displaystyle \frac{t^{k+1}-t}{t^{k+1}-t^k}\ us^{k} 
             + \frac{t-t^k}{t^{k+1}-t^k}\ us^{k+1},
       \quad t\in[t^{k},t^{k+1}), \quad k=0,\ldots,K-1,
   \end{array}
   @f]
   and subject to the constraints
   @f[
   \begin{array}{rcccll}
    \displaystyle \left\{ \frac{p_{min}}{p_{nominal}} \right. 
        &<& \displaystyle \frac{p}{p_{nominal}}
        &<& \left. \displaystyle \frac{p_{max}}{p_{nominal}} \right\}_i,
        \quad & i \in \mbox{find}(p_{active}), \\[3ex]
    \displaystyle \left\{ \frac{x^0_{min}}{x_{nominal}} \right. &<& x(t^0) 
        &<& \displaystyle \left. \frac{x^0_{max}}{x_{nominal}} \right\}_i, 
        \quad & i \in \mbox{find}(x^0_{active}), \\[3ex]
    && \{ x(t^0) &=& \displaystyle \frac{x^0}{x_{nominal}}
       \}_i, \quad & i \notin \mbox{find}(x^0_{active}), \\[3ex]
    && \left\{ \dot{x}(t^0) \right. &=& \left. 0 \right\}_i, 
        \quad & i \in \mbox{find}(x^0_{steady}).
   \end{array}
   @f]
   The problem is either treated with one stage and K sample periods or 
   with K stages one sample period per stage. In the latter case, i.e.
   treatment as multistage problem, additional K junction conditions
   (equality constraints) are introduced for the state variables x and
   the estimated parameters p. Discrete-time state variables with unknown
   initial value are introduced for the estimated parameters, in order to
   preserve a sparse structure.
   The treatment as multistage problem often leads to a faster convergence
   of the optimization solver.

   In addition to the estimated parameters and initial states, the 
   measurement matrix is calculated, which is defined as
   @f[
   M = \left[\begin{array}{cc}
       M_p^0, & M_{x^0}^0 \\[2ex]
       M_p^1, & M_{x^0}^1 \\[2ex]
       \vdots, & \vdots \\[2ex]
       M_p^K, & M_{x^0}^K
       \end{array}\right]
   @f]
   with
   @f[
   \begin{array}{r}
     \displaystyle M_p^k = \frac{d\,\displaystyle\frac{ys^k}{y_{nominal}}(I_y)}
                 {d\,\displaystyle\frac{p}{p_{nominal}}(I_p)}, \quad
     M_{x^0}^k = \frac{d\,\displaystyle\frac{ys^k}{y_{nominal}}(I_y)}
                 {d\,\displaystyle\frac{x^0}{x_{nominal}}(I_{x^0})}, \quad
     k=0,\ldots,K, \\[7ex]
     I_y = \mbox{find}(y_{active}),\ I_p = \mbox{find}(p_{active}),
     \ I_{x^0} = \mbox{find}(x^0_{active}).
   \end{array}
   @f]
   The measurement matrix can be used to obtain confidences for estimation
   results.
 */
class Prg_SFunctionEst: public Prg_SFunction {

 protected:
  int		_mdl_args_p_idx; ///< index of p in S-function parameters
  mxArray 	*_mx_p;  ///< model parameters for S-function

  int		_mdl_np;///< number of model parameters

  IVECP		_mdl_p_active;	///< indicate estimated parameters
  IVECP		_mdl_x0_active; ///< indicate estimated states
  IVECP		_mdl_x0_steady; ///< initial state should be in steady-state
  IVECP		_mdl_y_active; 	///< indicate measured outputs
  VECP		_mdl_y_nominal;	///< nominal output values (for scaling)

  VECP 		_mdl_p;		///< initial parameter values
  VECP 		_mdl_p_nominal;	///< nominal parameter values (for scaling)
  VECP 		_mdl_p_min;	///< lower bounds for estimated parameters
  VECP 		_mdl_p_max;	///< upper bounds for estimated parameters
  VECP 		_mdl_x_nominal; ///< nominal state values (for scaling)
  VECP 		_mdl_x0_min;	///< lower bounds for estimated initial states
  VECP 		_mdl_x0_max;	///< upper bounds for estimated initial states

  int		_nx;	///< number of states for optimizer
  int		_np;	///< number of estimated parameters
  int		_nx0;	///< number of estimated initial states
  int		_ny;	///< number of reference outputs for estimation
  bool 		_multistage; 	///< treat as multistage problem

  // vectors for inputs, states, and outputs
  MATP	_mdl_us;	///< given model inputs
  MATP	_mdl_ys;	///< calculated model outputs
  MATP	_prg_ys_ref;	///< reference values for active outputs

  // confidence things
  MATP 	_M;		///< measurement matrix M=dy/d(p,x0)
  MATP 	_dxdpx0;	///< help variable dx/d(p,x0)
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

 public:

  Prg_SFunctionEst();		///< constructor
  ~Prg_SFunctionEst();		///< destructor

  char *name() {return "SFunctionEst";} ///< name SFunctionEst
};  

#endif


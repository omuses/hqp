/*
 * Omu_IntODE.h --
 *   -- base class for integrating an ODE over a stage
 *   -- derived classes add a specific integration rule
 *
 * rf, 11/13/96
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

#ifndef Omu_IntODE_H
#define Omu_IntODE_H

#include "Omu_Integrator.h"

/**
 * Interface for standard solvers for Ordinary Differential Equations (ODE)
 * hiding details of sensitivity analysis. Sensitivity equations are
 * appended to the model equations.
 * This base class bypasses the high-level Omu_Integrator::solve method and
 * directly calls sys->continuous (optionally exploiting ADOL-C).
 * Nevertheless discrete-time states and control parameters are treated
 * in one vector of model parameters for derived classes.
 */
class Omu_IntODE: public Omu_Integrator {

 public:

  Omu_IntODE(); 	///< constructor
  ~Omu_IntODE(); 	///< destructor

  /**
   * @name Implementation of predefined methods.
   * @see Omu_Integrator
   */

  //@{

  void init(int k,
	    const Omu_StateVec &xc, const Omu_Vec &q,
	    const Omu_DependentVec &Fc, bool sa);

  void solve(int kk, double tstart, double tend,
	     const Omu_VariableVec &x, const Omu_VariableVec &u,
	     Omu_Program *sys, Omu_DependentVec &Ft, Omu_StateVec &xt);

  //@}

  /**
   * Interface to be implemented by derived ODE solvers.
   * y contains the initial solution and the result for values and
   * sensitivities (dim(y)=_n=_nxt-_nd without sensitivities,
   * dim(y)=_n*(1+_nxt-_nv+_nu) with sensitivities). u contains
   * discrete-time states and control parameters (dim(u)=_m=_nd+_nu).
   * @see Omu_Integrator for more details about the dimensions.
   */
  virtual void ode_solve(double tstart, VECP y, const VECP u, double tend) = 0;

  /** Callback routine for evaluating the model by derived ODE solvers. */
  void syseq(double t, const VECP y, const VECP u, VECP f);

 protected:

  /**
   * @name Depreciated methods (provided for old versions of derived classes).
   */
  //@{
  void init_stage(int k,
		  const Omu_States &x, const Omu_Vector &u,
		  bool sa = false) {}
  //@}

 private:

  /**
   * alternative syseq implementation calling high-level
   * continuous and using forward for derivatives
   */
  void syseq_forward(double t, const VECP y, const VECP u, VECP f);

  // backing store sys and current stage
  Omu_Program	*_sys;
  Omu_StateVec 	*_xt_ptr;
  Omu_DependentVec *_Ft_ptr;

  VECP		_y;
  VECP		_u;

  void		resize();

  // vectors and matrices for low level _sys->continuous callback
  Omu_Vec	_ut;
  Omu_SVec	_dxt;
  MATP		_Yx;
  MATP		_Yu;

  // variables for ADOL-C
  VECP		_x;
  VECP		_v;
  MATP		_X2;
  MATP		_Y2;
};  

#endif

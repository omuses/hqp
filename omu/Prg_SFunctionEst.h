/*
 * Prg_SFunctionEst.h -- 
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

//--------------------------------------------------------------------------
class Prg_SFunctionEst: public Prg_SFunction {

 protected:
  int		_mdl_args_p_idx; // index of p in S-function parameters
  mxArray 	*_mx_p;  // model parameters for S-function

  int		_mdl_np;// number of model parameters

  IVECP		_mdl_p_active;	// indicate estimated parameters
  IVECP		_mdl_x0_active; // indicate estimated states
  IVECP		_mdl_x0_steady; // initial state should be in steady-state
  IVECP		_mdl_y_active; 	// indicate measured outputs

  VECP 		_mdl_p;		// initial parameter values
  VECP 		_mdl_p_nominal;	// nominal parameter values (for scaling)
  VECP 		_mdl_p_min;	// lower bounds for estimated parameters
  VECP 		_mdl_p_max;	// upper bounds for estimated parameters
  VECP 		_mdl_x0;	// initial states
  VECP 		_mdl_x0_nominal;// nominal initial states (for scaling)
  VECP 		_mdl_x0_min;	// lower bounds for estimated initial states
  VECP 		_mdl_x0_max;	// upper bounds for estimated initial states

  int		_nx;	// number of states for optimizer
  int		_np;	// number of estimated parameters
  int		_nx0;	// number of estimated initial states
  int		_ny;	// number of reference outputs for estimation
  bool 		_multistage; 	// treat as multistage problem

  // vectors for inputs, states, and outputs
  MATP	_mdl_us;	// given model inputs
  MATP	_mdl_ys;	// calculated model outputs
  MATP	_prg_ys_ref;	// reference values for active outputs
  VECP	_prg_y_nominal;	// nominal values for active outputs (for scaling)

  // confidence things
  MATP 	_M;		// measurement matrix M=dy/d(p,x0)
  MATP 	_dxdpx0;
  MATP 	_dxfdpx0;

  // methods
  void setup_stages(IVECP ks, VECP ts);

  void setup(int k,
	     Omu_Vector &x, Omu_Vector &u, Omu_Vector &c);

  void setup_struct(int k,
		    const Omu_Vector &x, const Omu_Vector &u,
		    Omu_DependentVec &xt, Omu_DependentVec &F,
		    Omu_DependentVec &f,
		    Omu_Dependent &f0, Omu_DependentVec &c);

  void init_simulation(int k,
		       Omu_Vector &x, Omu_Vector &u);

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

 public:

  Prg_SFunctionEst();
  ~Prg_SFunctionEst();

  char *name() {return "SFunctionEst";}
};  

#endif


/*
 * Prg_SFunctionOpt.h -- 
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

//--------------------------------------------------------------------------
class Prg_SFunctionOpt: public Prg_SFunction {

 protected:
  VECP 		_mdl_x_nominal;	// nominal initial states (for scaling)
  VECP 		_mdl_x_min;	// lower bounds for estimated initial states
  VECP 		_mdl_x_max;	// upper bounds for estimated initial states

  IVECP		_mdl_u_active; 	// indicate optimized model inputs
  VECP 		_mdl_u_nominal;	// nominal inputs (for scaling)
  VECP 		_mdl_u_min;	// lower bounds for optimized inputs
  VECP 		_mdl_u_max;	// upper bounds for optimized inputs
  VECP 		_mdl_u_ref;	// reference to be reached for control inputs
  VECP 		_mdl_u_weight;	// weight for quadratic objective
  VECP 		_mdl_der_u_min;	// lower bounds for rate of change
  VECP 		_mdl_der_u_max;	// upper bounds for rate of change
  VECP 		_mdl_der_u_weight;// weight for rate of change of controls

  IVECP		_mdl_y_active; 	// indicate used outputs
  VECP 		_mdl_y_nominal;	// nominal outputs (for scaling)
  VECP 		_mdl_y_bias;	// bias correction (offset) for outputs 
  VECP 		_mdl_y_min;	// lower bounds for outputs
  VECP 		_mdl_y_max;	// upper bounds for outputs
  VECP 		_mdl_y_ref;	// reference to be reached for active outputs
  VECP 		_mdl_y_weight;	// weight for quadratic objective
  VECP 		_mdl_y_min_soft;	// lower soft bound
  VECP 		_mdl_y_min_weight1;	// weight for linear bound violation
  VECP 		_mdl_y_min_weight2;	// weight for quadratic bound violation
  VECP 		_mdl_y_max_soft;	// upper soft bound
  VECP 		_mdl_y_max_weight1;	// weight for linear bound violation
  VECP 		_mdl_y_max_weight2;	// weight for quadratic bound violation

  int		_nx;	// number of states for optimizer
  int		_nu;	// number of optimized control inputs
  int		_nc;	// number of constrained outputs
  int		_ns;	// number of slack variables for soft constraints

  int 		_sps;	// number of sample periods per stage

  // vectors for inputs, states, and outputs
  MATP	_mdl_us;	// given model inputs
  MATP	_mdl_ys;	// calculated model outputs

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

  Prg_SFunctionOpt();
  ~Prg_SFunctionOpt();

  char *name() {return "SFunctionOpt";}
};  

#endif


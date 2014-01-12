/*
 * Prg_CUTE_ST.h -- 
 *   - convert a CUTE problem into a Hqp_SqpProgram
 *
 * rf, 2/10/97
 */

/*
    Copyright (C) 1994--2014  Ruediger Franke

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

#ifndef Prg_CUTE_ST_H
#define Prg_CUTE_ST_H

#include "Hqp_SqpProgram.h"

/*
 * define types for external functions, 
 * that have to be provided for a problem
 */

#ifndef fint
#define fint int
#endif

#ifndef freal
#define freal double
#endif

#ifndef fbool
#define fbool int
#endif

/*
 * class declaration
 */
class Prg_CUTE_ST: public Hqp_SqpProgram {

 protected:
  // variables for the CUTE interface
  fint _N;		// number of variables
  fint _M;		// number of constraints
  fint _NNZJ;		// number of Jacobian non-zeros
  fint _NNZJMAX;	// dimension of _NNZJ
  fint *_INDVAR;	// variable index for constraints
  fint *_INDFUN;	// constraint index
  fint _NHI;		// number of Hessian entries
  fint _NEL;		// number of Hessian entries
  fint _NHIMAX;		// dimension of _HI
  fint _NELMAX;		// dimension of _IPRHI and _IPRNHI
  fint *_IPRHI;		// first Hessian entry of each finite-element in _H
  fint *_IRNHI;		// variable indices of Hessian entries
  fint *_IPRNHI;	// first variable index of each element in _IRNHI
  IVECP _iprhi;		// first Hessian entry of each finite-element in _H
  IVECP _irnhi;		// variable indices of Hessian entries
  IVECP _iprnhi;	// first variable index of each element in _irnhi

  double _xscale;// scaling of the variables
  double _fscale;// scaling of the objective function
  VECP _x0;	// initial variable values
  VECP _xs;	// variable vector for SIF
  VECP _var_lb;	// lower bounds for variables
  VECP _var_ub;	// upper bounds for variables
  VECP _cns_lb;	// lower bounds for constraints
  VECP _cns_ub;	// upper bounds for constraints
  double _Inf;	// number used for Infinity

  VECP _cns;	// values of constraints
  VECP _J;	// elements of sparse Jacobian
  VECP _H;	// elements of sparse Hessian
  VECP _v;	// Lagrangian multipliers

  // store for each variable and constraint the corresponding row in A or C
  IVECP _var_ridx;
  IVECP _cns_ridx;
  void parse_bounds(const VECP lb, const VECP ub, int n,
		    IVECP ridx, int &me, int &m);
  void update_bounds(const VECP val,
		     const VECP lb, const VECP ub, int n,
		     const IVECP ridx);

  int _nstages;	// number of stages for HQP
  int _n;	// number of variables for HQP
  int _me_bnd;	// number of variable bounds
  int _me_lin;	// number of linear equality constr. that come first in A
  int _m_lin;	// number of linear inequality constr. that come first in C

  IVECP _h2s;	// the SIF index for each HQP variable
  IVECP _s2h0;	// the HQP index for objective variables
  IVECP _s2h;	// the HQP index for each SIF variable per stage
  IVECP _s2h_stage; // the according stage
  IVECP _s2h_begin; // begin of a new variable in _s2h and _s2h_stage
  IVECP _equ_idx0;	// original variable index for blow up
  IVECP _equ_idx1;	// introduced variable index for blow up
  IVECP _grp_stage; 	// the stage that belongs to each nonlinear group
  IVECP _cns_stage; 	// the stage that belongs to each constraint

  bool _hela;
  bool _hela_init;
  bool _stretch;
  int _fbd_evals;

 public:

  Prg_CUTE_ST();
  ~Prg_CUTE_ST();

  void	setup();
  void	init_x();

  void	update_fbd();
  void	update(const VECP y, const VECP z);

  int	write_soln(IF_DEF_ARGS);	// write the solution

  const char *name() {return "CUTE_ST";}
};  


#endif



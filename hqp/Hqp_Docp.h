/*
 * Hqp_Docp.h --
 *   -- superclass: Hqp_SqpProgram
 *   -- superclass for formulating Discrete time Optimal Control Programs
 * rf, 11/12/94
 *
 * rf, 9/15/96
 *   -- introduce update_stage() for vals, grds, and hela
 *   -- separate update_bounds()
 *
 * rf, 1/31/97
 *   -- variable number of states/controls per stage
 *
 * rf, 11/01/00
 *   -- introduce new virtual method setup_horizon
 *   -- delete method init_vars
 *      (introduce additional arguments to setup_vars instead)
 *   -- provide default implementation for method init_simulation
 */

/*
    Copyright (C) 1994--2000  Ruediger Franke

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

#ifndef Hqp_Docp_H
#define Hqp_Docp_H

#include <math.h>

#include "Hqp_SqpProgram.h"

class Hqp_DocpAssoc;

//--------------------------------------
class Hqp_Docp: public Hqp_SqpProgram {

 public:

  Hqp_Docp();
  ~Hqp_Docp();

  void	horizon(int k0, int kf);

  // types of states for automatic setup of constraints
  static const int Periodical;	// equal initial and final states

  void	setup();		// init task for given _k0 and _kf
  void	init_x();		// set start variables

  int simulate(IF_DEF_ARGS);	// perform a simulation for current _x

  // provide abstract interface for an optimal control task
  // (default implementation of update_grds() uses finite differences)
  //------------------------------------------------------------------
  virtual void setup_horizon(int &k0, int &kf);

  virtual void setup_vars(int k,
			  VECP x, VECP xmin, VECP xmax,
			  VECP u, VECP umin, VECP umax,
			  VECP c, VECP cmin, VECP cmax) = 0;

  virtual void setup_struct(int k,
			    VECP f0x, VECP f0u, int &f0_lin,
			    MATP fx, MATP fu, IVECP f_lin,
			    MATP cx, MATP cu, IVECP c_lin,
			    MATP Lxx, MATP Luu, MATP Lxu);

  virtual void init_simulation(int k,
			       VECP x, VECP u);

  virtual void update_vals(int k, const VECP x, const VECP u,
			   VECP f, Real &f0, VECP c) = 0;

  virtual void update_grds(int k, const VECP x, const VECP u,
			   MATP fx, MATP fu, VECP f0x, VECP f0u,
			   MATP cx, MATP cu);

  virtual void update_hela(int k, const VECP x, const VECP u,
			   const VECP vx, const VECP vc,
			   MATP Lxx, MATP Luu, MATP Lxu);

  virtual void update_stage(int k, const VECP x, const VECP u,
			    VECP f, Real &f0, VECP c,
			    MATP fx, MATP fu, VECP f0x, VECP f0u,
			    MATP cx, MATP cu,
			    const VECP rf, const VECP rc,
			    MATP Lxx, MATP Luu, MATP Lxu);

  // setup_vars() calls alloc_vars() to allocate variables and bounds
  //----------------------------------------------------------------------
  void	alloc_vars(VECP v, VECP vmin, VECP vmax, int n);

 private:

  int	_k0;
  int	_kf;
  VECP	_xk;			// subvector into _x
  VECP	_uk;			// subvector into _x
  VECP	_fk;			// subvector into _x
  VECP	_s1;			// subvector 1
  VECP	_s2;			// subvector 2
  VEC	_xk_head;
  VEC	_uk_head;
  VEC	_fk_head;
  VEC	_s1_head;
  VEC	_s2_head;
  IVECP _x_type;

  int 	_fbd_evals;		// number of function evaluation

  // implement methods to define a Hqp_SqpProgram
  //---------------------------------------------
  void	setup_x();		// called by setup()
  void	setup_qp();		// called by setup()
  void	update_fbd();		// update values
  void	update(const VECP y, const VECP z); // update all
  void	reinit_bd();		// re-initialize bounds

  void	update_bounds();	// update state/control bounds and
				// initial/final state constraints

  void	parse_constr(const VECP cmin, const VECP cmax,
		     int idx, Hqp_DocpAssoc *eq,
		     Hqp_DocpAssoc *lb, Hqp_DocpAssoc *ub);

  // variables
  //----------
  IVECP _nxs;			// number of states per k
  IVECP _nus;			// number of controls per k
  VECP	_xus_init;		// initial values for x and u per k
  int	_granul;		// amount for growing _x_start during setup
  
  // objective
  //----------
  IVECP	_f0_lin;		// mark linear objective terms

  // state/control bounds
  //---------------------
  Hqp_DocpAssoc *_xu_eq;	// equality constraints
  Hqp_DocpAssoc *_xu_lb;	// lower bounds
  Hqp_DocpAssoc *_xu_ub;	// upper bounds
  IVECP	_f_lin;			// mark linear state equations
  IVECP	_f_start;		// start indizes into _f_lin for each k

  // additional constraints
  //-----------------------
  Hqp_DocpAssoc *_cns_eq;	// equality constraints
  Hqp_DocpAssoc *_cns_lb;	// lower bounds
  Hqp_DocpAssoc *_cns_ub;	// upper bounds
  IVECP	_cns_lin;		// mark linear constraints
  IVECP	_cns_start;		// start indizes into _cns_lin for each k
};  


#endif

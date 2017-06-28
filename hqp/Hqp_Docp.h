/**
 * @file Hqp_Docp.h 
 *   base class for formulating Discrete-time Optimal Control Programs
 *
 * rf, 11/12/94
 *
 */

/*
    Copyright (C) 1994--2017  Ruediger Franke

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

/** Base class for Discrete-time Optimal Control Programs */
class HQP_API Hqp_Docp: public Hqp_SqpProgram {

 public:

  Hqp_Docp(int ncpu = 1); 	///< constructor
  ~Hqp_Docp(); 			///< destructor

  void	horizon(int k0, int kf);///< setup

  // types of states for automatic setup of constraints
  static const int Periodical;	///< equal initial and final states

  void	setup();		///< init task for given _k0 and _kf
  void	init_x();		///< set start variables

  void simulate();		///< perform a simulation for current _x

  /**
   * @name Interface for an optimal control problem.
   */
  //@{
  virtual void setup_horizon(int &k0, int &kf);

  virtual void setup_vars(int k,
			  VECP x, VECP x_min, VECP x_max, IVECP x_int,
			  VECP u, VECP u_min, VECP u_max, IVECP u_int,
			  VECP c, VECP c_min, VECP c_max) = 0;

  virtual void setup_struct(int k, const VECP x, const VECP u,
			    MATP fx, MATP fu, IVECP f_lin,
			    VECP f0x, VECP f0u, int &f0_lin,
			    MATP cx, MATP cu, IVECP c_lin,
			    MATP Lxx, MATP Luu, MATP Lxu);

  virtual void init_simulation(int k,
			       VECP x, VECP u);

  virtual void update_vals(int k, const VECP x, const VECP u,
			   VECP f, Real &f0, VECP c) = 0;

  /** Update Jacobians.
   *  The default implementation obtains finite differences.
   */
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

  //@}

  /// setup_vars() should call alloc_vars() to allocate variables and bounds
  void	alloc_vars(VECP v, VECP v_min, VECP v_max, IVECP v_int, int n);

  /**
   * @name Member access methods.
   */
  //@{

  /** number of states per stage k */
  const IVECP nxs() const {return _nxs;}

  /** number of control parameters per stage k */
  const IVECP nus() const {return _nus;}

  /** number of function evaluations */
  int fbd_evals() const {return _fbd_evals;}

  /** number of cpus */
  int ncpu() const {return _ncpu;}
  void set_ncpu(int val) {
    if (val < 1 || val > (int)_f0->dim)
      m_error(E_SIZES, "ncpu must not exceed constructor value");
    _ncpu = val;
  }

  //@}

 private:

  int	_k0;
  int	_kf;
  int	_ncpu;
  VECP	*_xk;			// subvector into _x
  VECP	*_uk;			// subvector into _x
  VECP	*_fk;			// subvector into _x
  VECP	*_s1;			// subvector 1
  VECP	*_s2;			// subvector 2
  VECP	*_ck;
  VECP	_f0;
  MATP	*_fkx;
  MATP	*_fku;
  MATP	*_ckx;
  MATP	*_cku;
  MATP	*_Lkxx;
  MATP	*_Lkuu;
  MATP	*_Lkxu;
  VECP	*_vfk;
  VECP	*_vck;
  VEC	*_xk_head;
  VEC	*_uk_head;
  VEC	*_fk_head;
  VEC	*_s1_head;
  VEC	*_s2_head;
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
  IVECP _xus_integer; 		// integer settings for x and u per k
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
  IVECP	_xu_start;		// start indizes into _x for each k

  // additional constraints
  //-----------------------
  Hqp_DocpAssoc *_cns_eq;	// equality constraints
  Hqp_DocpAssoc *_cns_lb;	// lower bounds
  Hqp_DocpAssoc *_cns_ub;	// upper bounds
  IVECP	_cns_lin;		// mark linear constraints
  IVECP	_cns_start;		// start indizes into _cns_lin for each k
};  


#endif

/*
 * Hqp_SqpSolver.h -- 
 *   - base class for solving nonlinear programs with 
 *     sequential quadratic programming
 *
 * rf, 6/6/94
 *
 * rf, 12/23/99
 *   - extended member access methods
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

#ifndef Hqp_SqpSolver_H
#define Hqp_SqpSolver_H

#include <If_List.h>
#include <If_Class.h>

#include "Hqp_impl.h"

class Hqp_Program;
class Hqp_Solver;
class Hqp_SqpProgram;
class Hqp_HL;

IF_BASE_DECLARE(Hqp_SqpSolver);

class Hqp_SqpSolver {

 protected:
  If_List	_ifList;

  Hqp_SqpProgram *_prg;		// program to work with

  Hqp_Solver	*_solver;
  Hqp_HL	*_hela;
  VEC		*_x0;		// for line search
  VEC		*_xk;		// for line search
  VEC		*_d;		// x^(k+1) - x^k (including _alpha)
  VEC		*_y;		// multipiers for equality constraints
  VEC		*_z;		// multipiers for inequality constraints
  VEC		*_dL;		// grd_L^(k+1) - grd_L^k
  VEC		*_grd_L;
  Real		_phi;		// current value of penalty function
  Real		_dphi;		// (descending) direction of penalty function
  Real		_sQs;
  Real		_xQx;

  int		_iter;
  int		_max_iters;
  int		_inf_iters;
  int		_max_inf_iters;
  Real		_eps;
  Real		_norm_dx;
  Real		_norm_x;
  Real		_norm_inf;	// infeasibility norm
  Real		_norm_grd_L;
  Real		_norm_df;	// change of objective
  Real		_f_bak;		// objective of last iterate respectively
  Hqp_Result	_status;
  bool		_logging;	// print status messages

  Real		_alpha;		// step length
  Real		_min_alpha;	// step length, that identifies a stall

  bool		_hot_started;	// indicate a hot start
  SPMAT		*_qp_Q_hot;	// backing store Hessian of last cold start

  virtual void	feasible_vals();// update_xyz for suboptimal
  virtual void	update_vals()=0;// update _alpha,_d,_y,_z, and prg->x,f,b,d
  virtual VEC	*grd_L(const Hqp_Program *, VEC *out);
  Real		norm_inf(const Hqp_Program *);

 public:
  Hqp_SqpSolver();
  virtual ~Hqp_SqpSolver();

  Hqp_SqpProgram *prg() {return _prg;}
  void	set_prg(Hqp_SqpProgram *);

  virtual void	init();
  virtual void	qp_update();
  virtual void	qp_solve();
  virtual void	step();
  virtual void	hela_restart();
  virtual void	qp_reinit_bd();

  virtual void	solve();

  virtual const char *name() = 0;

  // member access methods

  int iter() {return _iter;}
  int max_iters() {return _max_iters;}
  void set_max_iters(int val) {_max_iters = val;}

  Real norm_dx() {return _norm_dx;}
  Real norm_inf() {return _norm_inf;}
  Real norm_grd_L() {return _norm_grd_L;}
};  


#endif

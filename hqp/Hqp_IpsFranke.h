/*
 * Hqp_IpsFranke.h -- 
 *   - interior point algorithm for solving quadratic programs
 *   - uses a Hqp_IpMatrix for performing iterates
 * rf, 5/28/94
 *
 * rf, 9/14/96
 *   - extended interface to Hqp_IpMatrix
 *
 * rf, 8/13/98
 *   - rename Hqp_IpSolver to Hqp_IpsFranke
 *   - make Hqp_IpsFranke an exchangeable interface class
 */

/*
    Copyright (C) 1994--1998  Ruediger Franke

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

#ifndef Hqp_IpsFranke_H
#define Hqp_IpsFranke_H

#include "Hqp_Solver.h"

class Hqp_Program;
class Hqp_IpMatrix;

class Hqp_IpsFranke: public Hqp_Solver {

 protected:
  int		_n, _me, _m;
  VEC		*_w;		// vector of slacks
  VEC		*_r1;		// right hand sides
  VEC		*_r2;
  VEC		*_r3;
  VEC		*_r4;
  VEC		*_a1;		// correcture vectors
  VEC		*_a2;
  VEC		*_a3;
  VEC		*_dx;		// modifications for one step
  VEC		*_dy;
  VEC		*_dz;
  VEC		*_dw;
  Real		_mu0;		// used for cold start
  Real		_Ltilde;
  Real		_zeta;
  Real		_alpha;		// step width
  Real		_alphabar;
  Real		_beta;
  Real		_rhomin;
  Real		_gap;		// duality gap

  Hqp_IpMatrix	*_matrix;

  int		_hot_started;	// flag
  int		_fail_iters;	// lost iters after a failed warm start
  int		_max_warm_iters;

 public:

  Hqp_IpsFranke();
  ~Hqp_IpsFranke();

  // initializing a program
  int		init(IF_DEF_ARGS);
  int		update(IF_DEF_ARGS);

  // solving a program
  int		cold_start(IF_DEF_ARGS);
  int		hot_start(IF_DEF_ARGS);
  int		step(IF_DEF_ARGS);
  int	 	solve(IF_DEF_ARGS);

  // member access
  Real		gap() {return _gap;}
  Real		zeta() {return _zeta;}

  char	*name() {return "Franke";}
};  

#endif



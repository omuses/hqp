/*
 * Hqp_IpsMehrotra.h -- 
 *   - interior point algorithm for solving quadratic programs
 *   - uses a Hqp_IpMatrix for performing iterates
 * rf, 5/28/94
 *
 * rf, 9/14/96
 *   - extended interface to Hqp_IpMatrix
 * ea, 5/9/97
 *   - Mehrotra's primal-dual predictor-corrector method for QP problems.
 *     S. J. Wright: Primal-dual interior-point methods.
 *                   SIAM, Philadelphia, 1997.
 */

/*
    Copyright (C) 1994--2014  Eckhard Arnold and Ruediger Franke

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

#ifndef Hqp_IpsMehrotra_H
#define Hqp_IpsMehrotra_H

#include "Hqp_Solver.h"

class Hqp_Program;
class Hqp_IpMatrix;

class Hqp_IpsMehrotra: public Hqp_Solver {

 protected:
  int		_n, _me, _m;
  VEC		*_w;		// vector of slacks
  VEC		*_r1;		// right hand sides
  VEC		*_r2;
  VEC		*_r3;
  VEC		*_r4;
  VEC		*_dx;		// modifications for one step
  VEC		*_dy;
  VEC		*_dz;
  VEC		*_dw;
  VEC		*_dxa;		// affine step
  VEC		*_dya;
  VEC		*_dza;
  VEC		*_dwa;
  VEC           *_z_hot;        // hot start
  VEC           *_w_hot;        
  VEC           *_phimin;       // convergence check
  VEC		*_d1;
  VEC		*_d2;
  Real		_mu0;		// used for cold start
  Real		_alpha;		// step length
  Real          _gammaf;
  Real		_gap;		// duality gap
  Real          _test;          // KKT residual
  Real          _norm_r0;
  Real          _norm_data;

  Hqp_IpMatrix	*_matrix;

  int		_hot_started;	// flag
  int		_fail_iters;	// lost iters after a failed warm start
  int		_max_warm_iters;
  int           _logging;       // output
  int           _init_method;   // initialization method

 public:

  Hqp_IpsMehrotra();
  ~Hqp_IpsMehrotra();

  // initializing a program
  void		init();
  void		update();

  // solving a program
  void		cold_start();
  void		hot_start();
  void		step();
  void	 	solve();

  // member access
  Real		gap() {return _gap;}
  //  Real		zeta() {return _zeta;}

  const char *name() {return "Mehrotra";}
};  

#endif

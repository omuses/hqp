/*
 * Hqp_SqpPowell.h -- 
 *   - Powells SQP algorithm
 *
 * rf, 6/8/94
 */

/*
    Copyright (C) 1994--2002  Ruediger Franke

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

#ifndef Hqp_SqpPowell_H
#define Hqp_SqpPowell_H

#include "Hqp_SqpSolver.h"

class Hqp_SqpProgram;


class Hqp_SqpPowell: public Hqp_SqpSolver {

 protected:
  bool	_relaxed;
  int	_watchdog_start;  // first allowed watchdog iteration
  int	_watchdog_credit; // allowed number of bad iterations
  int	_watchdog_iter;   // iteration for backtracking
  bool  _watchdog_logging;
  bool	_damped_multipliers;
  Real	_phil;
  Real	_phil_test;
  VEC	*_re;	// penalty coeffizients for equality constraints
  VEC	*_r;	// penalty coeffizients for inequality constraints
  VEC	*_sy_y;
  VEC	*_sz_z;
  VEC	*_x0;
  VEC	*_y0;
  VEC	*_z0;
  VEC	*_xk;
  VEC	*_xl;	// backing store for watchdog
  VEC	*_qp_xl;
  VEC	*_yl;
  VEC	*_zl;

  Real	phi();
  Real	phi1();
  VEC	*update_r(const VEC *z, VEC *r);

  void	update_vals();
  Real	val_L();

 public:
  Hqp_SqpPowell();
  ~Hqp_SqpPowell();

  int	init(IF_DEF_ARGS);

  char	*name() {return "Powell";}
};  


#endif

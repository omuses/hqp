/*
 * Hqp_Solver.h -- 
 *   - abstract base class for a quadratic solver
 *
 * rf, 7/19/94
 *
 * rf, 8/13/98
 *   - make Hqp_Solver a visible interface class
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

#ifndef Hqp_Solver_H
#define Hqp_Solver_H

#include <If_List.h>
#include <If_Class.h>

#include "Hqp_impl.h"

class Hqp_Program;

class Hqp_Solver {

 protected:

  If_List	_ifList;

  Hqp_Program	*_qp;		// program to work with
  Hqp_Result	_result;	// current processing state

  int		_iter;
  int		_max_iters;
  Real		_eps;

  VEC		*_y;		// multiplicators for equalities
  VEC		*_z;		// multiplicators for inequalities

 public:

  Hqp_Solver();
  virtual ~Hqp_Solver();

  void		qp(Hqp_Program *);

  // initializing a program (updating for SQP integration)
  virtual void	init() = 0;
  virtual void	update() = 0;

  // solving a program (hot start for SQP integration)
  virtual void	cold_start() = 0;
  virtual void	hot_start() = 0;
  virtual void	step() = 0;
  virtual void	solve() = 0;

  // member access
  const char 	*qp_result() const;
  Hqp_Result	result() {return _result;}
  int		iter() {return _iter;}
  int		max_iters() {return _max_iters;}
  void		max_iters(int n_max_iters) {_max_iters = n_max_iters;}
  Real		eps() {return _eps;}
  void		eps(Real n_eps) {_eps = n_eps;}
  const VEC	*y() {return _y;}
  const VEC	*z() {return _z;}

  // interface name
  virtual const char *name()=0;
};  

HQP_API IF_BASE_DECLARE(Hqp_Solver);

#endif

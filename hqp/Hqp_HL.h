/*
 * Hqp_HL.h -- 
 *   - base class for approximating sparse Hessian of Lagrangian
 *
 * rf, 7/19/94
 *
 * rf, 9/15/96
 *   - extended init() to additionally backup f, b, d
 *
 * rf, 2/13/97
 *   - new method posdef() that performs Gerschgorin modification
 *
 * rf, 2/18/98
 *   - update_Q() gets sqp_alpha as additional argument
 *
 * E. Arnold, 2001-08-28
 *   - different variants of initial scaling (_scale, _init_multipliers)
 *
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

#ifndef Hqp_HL_H
#define Hqp_HL_H

#include <If_List.h>
#include <If_Class.h>

#include "Hqp_impl.h"

IF_BASE_DECLARE(Hqp_HL);

class Hqp_Program;
class Hqp_SqpProgram;

class Hqp_HL {

 protected:
  If_List	_ifList;
  int		_scale;         // initial scaling (0,1,2,3)
  bool 		_init_multipliers;       
  Real		_eps;		// ensure positive definiteness
  bool          _logging;	// print status messages
  VEC		*_rowsum;				     

  VEC 	*grd_L(const VEC *y, const VEC *z, const Hqp_Program *qp,
	       VEC *out);
  
 public:
  Hqp_HL();
  virtual ~Hqp_HL();

  virtual void	setup(Hqp_SqpProgram *) = 0;
  virtual void	init(const VEC *y, const VEC *z, Hqp_SqpProgram *);
  virtual void	update(const VEC *s, const VEC *u, Real alpha,
		       Hqp_SqpProgram *) = 0;

  void est_y(Hqp_SqpProgram *, VEC *);
  virtual void	posdef(Hqp_SqpProgram *);

  virtual char *name() = 0;
};  


#endif

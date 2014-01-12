/*
 * Hqp_IpMatrix.h -- abstract base class
 *   - manage the Jacobian matrix of Interior Point algorithms
 * rf, 6/1/94
 *
 * rf, 9/14/96
 *   - extended solve-interface
 *
 * rf, 3/3/97
 *   - introduce iterative refinement of solutions
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

#ifndef Hqp_IpMatrix_H
#define Hqp_IpMatrix_H

#include <If_List.h>
#include <If_Class.h>

#include "Hqp_impl.h"

class Hqp_Program;

IF_BASE_DECLARE(Hqp_IpMatrix);

class Hqp_IpMatrix {

 protected:
  If_List _ifList;

  VEC	*_r1;
  VEC	*_r2;
  VEC	*_r3;
  VEC	*_r4;

  VEC	*_dx;
  VEC	*_dy;
  VEC	*_dz;
  VEC	*_dw;

  Real 	_eps;

 public:
  Hqp_IpMatrix();
  virtual ~Hqp_IpMatrix();
  
  virtual void	init(const Hqp_Program *) = 0;
  virtual void	update(const Hqp_Program *) = 0;

  virtual void	factor(const Hqp_Program *,
		       const VEC *z, const VEC *w) = 0;

  // default solve() implements iterative refinement
  virtual Real	solve(const Hqp_Program *,
		      const VEC *z, const VEC *w,
		      const VEC *r1, const VEC *r2, const VEC *r3,
		      const VEC *r4,
		      VEC *dx, VEC *dy, VEC *dz, VEC *dw);

  virtual Real	residuum(const Hqp_Program *,
			 const VEC *z, const VEC *w,
			 const VEC *r1, const VEC *r2, const VEC *r3,
			 const VEC *r4,
			 VEC *dx, VEC *dy, VEC *dz, VEC *dw);

  virtual void	step(const Hqp_Program *,
		     const VEC *z, const VEC *w,
		     const VEC *r1, const VEC *r2, const VEC *r3,
		     const VEC *r4,
		     VEC *dx, VEC *dy, VEC *dz, VEC *dw);

  virtual const char *name()=0;
};

#endif

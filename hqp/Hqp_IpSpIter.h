/*
 * Hqp_IpSpIter.h --
 *   - manage the Jacobian matrix of Interior Point algorithms
 *   - eliminate inequality constraints
 *   - use sparse matrix and iterative solution
 *
 * rf, 9/13/94
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

#ifndef Hqp_IpSpIter_H
#define Hqp_IpSpIter_H

extern "C" {
#include <iter.h>
}

#include "Hqp_IpMatrix.h"

class Hqp_Structure;

class Hqp_IpSpIter: public Hqp_IpMatrix {

 protected:
  int		_n, _me, _m;	// dimensions
  SPMAT		*_CT;
  SPMAT		*_J;
  SPMAT		*_J_raw;
  SPMAT		*_J_fct;
  PERM		*_QP2J;
  PERM		*_J2QP;
  PERM		*_pivot;
  PERM		*_blocks;
  VEC		*_zw;
  VEC		*_scale;
  VEC		*_r12;
  VEC		*_xy;
  VEC		*_test;

  Hqp_Structure	*_CTC_structure;
  SPMAT		*sub_CTC(const PERM *, SPMAT *);

  ITER		*_iter;
  static VEC	*iter_solve(void *cld, VEC *x, VEC *out);

 public:
  		Hqp_IpSpIter();
  virtual 	~Hqp_IpSpIter();
  
  virtual void	init(const Hqp_Program *);
  virtual void	update(const Hqp_Program *);

  virtual void	factor(const VEC *z, const VEC *w);
  virtual Real	solve(const VEC *r1, const VEC *r2, const VEC *r3,
		      VEC *dx, VEC *dy, VEC *dz);

  const char *name() {return "SpIter";}
};

#endif

/*
 * Hqp_IpFullBdLU.h --
 *   - manage the Jacobian matrix of Interior Point algorithms
 *   - use band matrix and LU factorization
 *   - only one "solve" is supported per "factorize" (because of bdLUsolve)
 *
 * rf, 6/1/94
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

#ifndef Hqp_IpFullBdLU_H
#define Hqp_IpFullBdLU_H

#include "Hqp_IpMatrix.h"


class Hqp_IpFullBdLU: public Hqp_IpMatrix {

 protected:
  int		_n, _me, _m;	// dimensions
  BAND		*_J;
  BAND		*_J_raw;
  BAND		*_J_fct;
  PERM		*_QP2J;
  PERM		*_J2QP;
  PERM		*_pivot;
  VEC		*_scale;
  VEC		*_r123;
  VEC		*_xyz;
  VEC		*_test;

 public:
  		Hqp_IpFullBdLU();
  virtual 	~Hqp_IpFullBdLU();
  
  virtual void	init(const Hqp_Program *);
  virtual void	update(const Hqp_Program *);

  virtual void	factor(const VEC *z, const VEC *w);
  virtual Real	solve(const VEC *r1, const VEC *r2, const VEC *r3,
		      VEC *dx, VEC *dy, VEC *dz);

  const char *name() {return "FullBdLU";}
};

#endif
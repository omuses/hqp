/*
 * Hqp_IpRedSpBKP.h --
 *   - manage the Jacobian matrix of Interior Point algorithms
 *   - eliminate inequality constraints
 *   - use sparse matrix and BKP factorization
 *
 * rf, 9/13/94
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

#ifndef Hqp_IpRedSpBKP_H
#define Hqp_IpRedSpBKP_H

#include "Hqp_IpMatrix.h"


class Hqp_IpRedSpBKP: public Hqp_IpMatrix {

 protected:
  int		_n, _me, _m;	// dimensions
  int		_sbw;
  Real		_tol;		// tolerance for fill-in vs. stability
  SPMAT		*_CT;
  SPMAT		*_J;
  SPMAT		*_J_raw;
  PERM		*_QP2J;
  PERM		*_J2QP;
  PERM		*_pivot;
  PERM		*_blocks;
  VEC		*_zw;
  VEC		*_scale;
  VEC		*_r12;
  VEC		*_xy;

  IVEC		*_CTC_degree;
  IVEC		*_CTC_neigh_start;
  IVEC		*_CTC_neighs;

  SPMAT		*sub_CTC(const PERM *, SPMAT *);

 public:
  	Hqp_IpRedSpBKP();
   	~Hqp_IpRedSpBKP();
  
  void	init(const Hqp_Program *);
  void	update(const Hqp_Program *);

  void	factor(const Hqp_Program *, const VEC *z, const VEC *w);
  void	step(const Hqp_Program *, const VEC *z, const VEC *w,
	     const VEC *r1, const VEC *r2, const VEC *r3, const VEC *r4,
	     VEC *dx, VEC *dy, VEC *dz, VEC *dw);

  char	*name() {return "RedSpBKP";}
};

#endif

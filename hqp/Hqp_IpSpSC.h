/*
 * Hqp_IpSC.h --
 *   - manage the Jacobian matrix of Interior Point algorithms
 *   - eliminate inequality constraints
 *   - use sparse matrix and solution of schur complement system
 *
 * hl, 96/10/08
 */

/*
    Copyright (C) 1996--2014  Hartmut Linke

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

#ifndef Hqp_IpSpSC_H
#define Hqp_IpSpSC_H

#include "Hqp_IpMatrix.h"


class Hqp_IpSpSC: public Hqp_IpMatrix {

 protected:
  int		_n, _me, _m;	    // dimensions
  int		_sbw;
  int		_piv_QCVC_flag;
  int		_piv_AQCVCA_flag;
  
  double        _macheps;

  SPMAT		*_Q;
  SPMAT		*_CT;
  SPMAT		*_AT;
  SPMAT		*_QCVC;
  SPMAT		*_QCVCA;
  SPMAT		*_AQCVC;
  SPMAT		*_AQCVCA;
  SPMAT		*_AQCVCA_fac;

  PERM	        *_piv_QCVC;
  PERM	        *_piv_AQCVCA;

  VEC           *_v1;
  VEC           *_v2;
  VEC		*_zw;
  VEC		*_scale;
  VEC		*_r12;
  VEC		*_xy;

  IVEC		*_CTC_degree;
  IVEC		*_CTC_neigh_start;
  IVEC		*_CTC_neighs;

  SPMAT		*sub_CTC(const PERM *, SPMAT *);

 public:

   Hqp_IpSpSC();
   ~Hqp_IpSpSC();
  
  void	init(const Hqp_Program *);
  void	update(const Hqp_Program *);

  void	factor(const Hqp_Program *, const VEC *z, const VEC *w);
  Real	solve(const Hqp_Program *, const VEC *z, const VEC *w,
	      const VEC *r1, const VEC *r2, const VEC *r3, const VEC *r4,
	      VEC *dx, VEC *dy, VEC *dz, VEC *dw);

  const char *name() {return "SpSC";}
};

#endif

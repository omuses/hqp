/*
 * Hqp_HL_BFGS.h -- 
 *   - BFGS Hessian approximation with Powell's modification
 *   - computes a separate update for every diagonal block
 *
 * rf, 7/19/94
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

#ifndef Hqp_HL_BFGS_H
#define Hqp_HL_BFGS_H

#include "Hqp_HL.h"


class Hqp_HL_BFGS: public Hqp_HL {

 protected:
  Real _gamma;
  bool _eigen_control;
  int  _bsize;
  int  _max_bsize;

  VEC	*_v;
  VEC	*_sQ;
  VEC	*_Qs;

  MAT	*_b_Q;
  int 	_b_begin;

  void	update_b_Q(const VEC *s, const VEC *u, Real alpha, MAT *Q);
  int	next_block(const SPMAT *, int *offs, int *size);
  
 public:
  Hqp_HL_BFGS();
  ~Hqp_HL_BFGS();

  void setup(Hqp_SqpProgram *);
  void update(const VEC *s, const VEC *u, Real alpha, Hqp_SqpProgram *);

  char *name() {return "BFGS";}
};  


#endif

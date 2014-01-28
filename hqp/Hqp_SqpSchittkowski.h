/*
 * Hqp_SqpSchittkowski.h -- 
 *   - Schittkowski's SQP algorithm
 *
 * rf, 6/8/94
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

#ifndef Hqp_SqpSchittkowski_H
#define Hqp_SqpSchittkowski_H

#include "Hqp_SqpSolver.h"

class Hqp_SqpProgram;


class Hqp_SqpSchittkowski: public Hqp_SqpSolver {

 protected:
  VEC		*_x0;	// start vector for line search
  VEC		*_xk;	// currently treated vector
  VEC		*_re;	// penalty coeffizients for equality constraints
  VEC		*_r;	// penalty coeffizients for inequality constraints
  VEC		*_ve;
  VEC		*_v;
  VEC		*_v0;
  VEC		*_ve0;
  VEC		*_ue_ve;
  VEC		*_u_v;
  VEC		*_sgme;
  VEC		*_sgm;

  Real		_mu;
  Real		_eps;
  Real		_beta;
  bool		_damped_multipliers;

  Real		phi();
  Real		dphi();
  VEC		*update_sgm(const VEC *r, VEC *sgm);
  VEC		*update_r(const VEC *u, const VEC *v, const VEC *sgm,
			  Real dQd, VEC *r);

  void		update_vals();

 public:
  Hqp_SqpSchittkowski();
  ~Hqp_SqpSchittkowski();

  void	init();
  
  const char *name() {return "Schittkowski";}
};  


#endif

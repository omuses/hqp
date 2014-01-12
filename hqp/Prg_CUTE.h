/*
 * Prg_CUTE.h -- 
 *   - convert a CUTE problem into a Hqp_SqpProgram
 *
 * rf, 10/27/95
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

#ifndef Prg_CUTE_H
#define Prg_CUTE_H

#include "Hqp_SqpProgram.h"

/*
 * define types for external functions, 
 * that have to be provided for a problem
 */

#ifndef fint
#define fint int
#endif

#ifndef freal
#define freal double
#endif

#ifndef fbool
#define fbool int
#endif


/*
 * class declaration
 */
class Prg_CUTE: public Hqp_SqpProgram {

 protected:
  // variables for the CUTE interface
  fint _N;
  fint _M;
  fint _NNZJ;
  fint _NNZH;
  fint _NNZJMAX;
  fint _NNZHMAX;
  fint *_INDVAR;
  fint *_INDFUN;
  fint *_IRNSH;
  fint *_ICNSH;

  double _fscale;// scaling of the objective function
  VEC *_x0;	// initial variable values
  VEC *_var_lb;	// lower bounds for variables
  VEC *_var_ub;	// upper bounds for variables
  VEC *_cns_lb;	// lower bounds for constraints
  VEC *_cns_ub;	// upper bounds for constraints
  Real _Inf;	// number used for Infinity

  VEC *_cns;	// values of constraints
  VEC *_J;	// elements of sparse Jacobian
  VEC *_H;	// elements of sparse Hessian
  VEC *_v;	// lagrangian multipliers

  // store for each variable and constraint corresponding row in A or C
  IVEC *_var_ridx;
  IVEC *_cns_ridx;
  void parse_bounds(const Real *l_ve, const Real *u_ve, int n,
		    int *ridx_ive, int *me, int *m);
  void update_bounds(const Real *val_ve,
		     const Real *l_ve, const Real *u_ve, int n,
		     const int *ridx_ive);

  // store number of linear constraints, that come first in A and C
  int _me_lin;
  int _m_lin;

  bool _hela;		// update Lagrangian Hessian
  bool _hela_init;	// initialize Lagrangian Hessian
  int _fbd_evals;

 public:

  Prg_CUTE();
  ~Prg_CUTE();

  void	setup();
  void	init_x();

  void	update_fbd();
  void	update(const VECP y, const VECP z);

  int	write_soln(IF_DEF_ARGS);	// write the solution

  const char *name() {return "CUTE";}
};  


#endif



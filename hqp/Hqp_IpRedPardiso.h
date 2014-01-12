/*
 * Hqp_IpRedPardiso.h --
 *   - manage the Jacobian matrix of Interior Point algorithms
 *   - eliminate inequality constraints
 *   - use Pardiso sparse matrix solver
 *
 * hl,  2006/11/22
 *
 * derived from Hqp_IpRedSpBKP.h (developed by Ruediger Franke)
 *
 */

/*
    Copyright (C) 2006    Hartmut Linke

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

#ifndef Hqp_IpRedPardiso_H
#define Hqp_IpRedPardiso_H

#include "Hqp_IpMatrix.h"


class Hqp_IpRedPardiso: public Hqp_IpMatrix {

 protected:
  int		_n, _me, _m, _nnz;	// dimensions
  int		_sbw;
  Real		_tol;		// tolerance for fill-in vs. stability
  SPMAT		*_CT;
  SPMAT		*_QCTC;
  PERM		*_pivot;
  PERM		*_blocks;
  VEC		*_zw;
  VEC		*_scale;
  VEC		*_r12;
  VEC		*_xy;

  IVEC		*_CTC_degree;
  IVEC		*_CTC_neigh_start;
  IVEC		*_CTC_neighs;

  // parameters of PARDISO solver

  void	        *_pardiso_pt[64];
  int           _pardiso_parm[64];

  int           _reinit;
  int           _pivot_strategy;

  IVEC          *_iv;
  IVEC          *_jv;
  IVEC          *_ipivot;
  VEC           *_v;
  VEC           *_v_raw;

  int           _maxfct;        // Maximum number of numerical factorizations
  int           _mnum;          // Which factorization to use
  int           _msglvl;        // Print statistical information in file
  int           _error;         // Initialize error flag
  int           _mtype;         // Real symmetric matrix
  int           _nrhs;          // Number of right hand sides
  int           _phase;         // solution phase of the PARDISO solver
  int           _dim;           // dimension of  equation system 

  void		    sub_CTC();
  void          free_pardiso();
  void          reinit_pardiso();

 public:
  	Hqp_IpRedPardiso();
   	~Hqp_IpRedPardiso();
  
  void	init(const Hqp_Program *);
  void	update(const Hqp_Program *);

  void	factor(const Hqp_Program *, const VEC *z, const VEC *w);
  void	step(const Hqp_Program *, const VEC *z, const VEC *w,
	     const VEC *r1, const VEC *r2, const VEC *r3, const VEC *r4,
	     VEC *dx, VEC *dy, VEC *dz, VEC *dw);

  const char *name() {return "RedPardiso";}
};

#endif

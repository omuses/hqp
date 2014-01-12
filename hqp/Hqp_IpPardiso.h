/**
 * @file Hqp_IpPardiso.h 
 *   solve the Jacobian matrix of Interior Point algorithms using
 *   the Parallel direct solver Pardiso
 *
 * hl, 2006/11/22
 *
 * derived from Hqp_IpSpBKP.h (developed by Ruediger Franke)
 *
 */

/*
    Copyright (C) 2006   Hartmut Linke

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

#ifndef Hqp_IpPardiso_H
#define Hqp_IpPardiso_H

#include "Hqp_IpMatrix.h"
#include "Hqp_DynLoad.h"

/** Signature of used PARDISO function (version 3) */
typedef int (pardiso_ft)
	(void *, int *, int *, int *, int *, int *,
	double *, int *, int *, int *, int *, int *,
	int *, double *, double *, int *);

/**
   Solve the Jacobian matrix of Interior Point algorithms 
   using the Pardiso solver. 

   The solver is loaded dynamically at runtime from the 
   Intel Math Kernel Library per default (tested with MKL 10.2.5.035).
   The library containing the solver and the function name can be
   configured using mat_pardiso_libname and mat_pardiso_funcname, respectively.
 */
class Hqp_IpPardiso: public Hqp_IpMatrix {

 protected:
  int		_n, _me, _m, _nnz; ///< dimensions
  int		_sbw;           ///< semi band width of _J
  Real		_tol;           ///< tolerance for fill-in vs. stability
  PERM		*_pivot;
  VEC		*_scale;
  VEC		*_r123;
  VEC		*_xyz;

  /**
   * @name Parameters of PARDISO solver
   */
  //@{

  Hqp_DynLoad 	_dl;            ///< dynamic load of solver library
  char          *_pardiso_libname; ///< name of library containing solver
  char          *_pardiso_funcname;///< name of Pardiso function
  pardiso_ft    *_pardiso_fp;   ///< pointer to Pardiso function
  int 		_nparallel;     ///< number of processor cores to use

  void	        *_pardiso_pt[64];
  int           _pardiso_parm[64];

  int           _reinit;

  IVEC          *_iv;
  IVEC          *_jv;
  VEC           *_v;
  VEC           *_v_raw;

  int           _maxfct;        ///< Maximum number of numerical factorizations
  int           _mnum;          ///< Which factorization to use
  int           _msglvl;        ///< Print statistical information in file
  int           _error;         ///< Initialize error flag
  int           _mtype;         ///< Real symmetric matrix
  int           _nrhs;          ///< Number of right hand sides
  int           _phase;         ///< solution phase of the PARDISO solver
  int           _dim;           ///< dimension of  equation system 

  //@}

  void          free_pardiso();
  void          reinit_pardiso();

 public:
        Hqp_IpPardiso();
   	~Hqp_IpPardiso();
  
  void	init(const Hqp_Program *);
  void	update(const Hqp_Program *);

  void	factor(const Hqp_Program *, const VEC *z, const VEC *w);
  void	step(const Hqp_Program *, const VEC *z, const VEC *w,
	     const VEC *r1, const VEC *r2, const VEC *r3, const VEC *r4,
	     VEC *dx, VEC *dy, VEC *dz, VEC *dw);

  /**
   * @name Member access methods
   */
  //@{

  /// library containing solver
  const char *pardiso_libname() const {return _pardiso_libname;}
  /// set library containing solver
  void set_pardiso_libname(const char *value);

  /// name of pardiso function
  const char *pardiso_funcname() const {return _pardiso_funcname;}
  /// set name of pardiso function
  void set_pardiso_funcname(const char *value);

  /// number of processor cores to use
  int nparallel() const {return _nparallel;}
  /// set number of processor cores to use
  void set_nparallel(int value) {_nparallel = value;}

  //@}

  const char *name() {return "Pardiso";}
};

#endif

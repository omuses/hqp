/*
 * Omu_IntDASPK.h --
 *   -- integrate DAE over a stage using DASPK3.0
 *      (derived from earlier Omu_IntDASPKSO implementation)
 *
 * rf, 10/2/96
 */

/*
    Copyright (C) 1996--2002  Ruediger Franke and Hartmut Linke

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

#ifndef Omu_IntDASPK_H
#define Omu_IntDASPK_H

#include "Omu_Integrator.h"

/*
 * FORTRAN data types
 */

#ifdef fint
#undef fint
#endif
#define fint int

#ifdef freal
#undef freal
#endif
#define freal double

/**
 * Solve differential-algebraic equation system using DASPK.
 * Currently DASPK version 3.0 is supported (distribution file
 * http://www.engineering.ucsb.edu/~cse/Software/ddaspk30.tar.gz).
 */
class Omu_IntDASPK: public Omu_Integrator {

 public:

  Omu_IntDASPK(); 	///< constructor
  ~Omu_IntDASPK(); 	///< destructor

  /**
   * @name Implementation of predefined methods.
   * @see Omu_Integrator
   */

  //@{

  char *name() {return "DASPK";}

  void init_stage(int k,
		  const Omu_States &x, const Omu_Vector &u,
		  bool sa);

  void solve(int kk, double tstart, double tend,
	     const Omu_States &x, const Omu_Vector &u,
	     Omu_Program *sys, Omu_DepVec &Fc, Omu_SVec &xc);

  //@}

  /**
   * @name Callback routines for DASPK.
   */
  //@{

  void res(freal *t, freal *x, freal *xprime,
	   freal *delta, fint *ires, freal *rpar, fint *ipar,
	   freal *senpar);

  void jac(freal *t, freal *y, freal *yprime,
	   freal *pd, freal *cj, freal *rpar, fint *ipar,
	   freal *senpar, fint *ijac);

  //@}

 protected:

  /**
   * User given semi-bandwidth of Jacobian (default: -1).
   * If a value >= 0 is specified, then it is used instead
   * of the automatic detection.
   */
  int		_jac_sbw;

  /**
   * User specification to allow banded solver (default: true).
   * Banded solvers are used if _banded is true, _jac_sbw < 0,
   * and if the automatic detection indicates that the problem
   * can be solved more efficiently in this way.
   */
  bool		_banded;
  bool		_banded_solver;	///< internal flag for using banded solver

  /**
   * User specification if Krylov iterative solver should be used
   * instead of direct solver (default: false).
   */
  bool		_krylov;

  /**
   * User specification if a preconditioner should be used
   * by the Krylov iterative solver (default: true).
   * A banded preconditioner is used if the problem is treated banded.
   * Otherwise an incomplete LU factorization is used.
   */
  bool		_krylov_prec;

  /**
   * User specification if externally provided Jacobian should be used
   * (default: true). Otherwise DASPK approximates it internally,
   * which is generally more efficient for not analytically given Jacobians.
   */
  bool 		_with_jac;

  /**
   * User defined fixed number of steps.
   * (default: 0, i.e. variable step size based on _rtol and _atol)
   */
  int		_nsteps;

 private:

  void		realloc();
  void		init_options(const Omu_States &x);

  // backing store sys and vector of dependent variables for callbacks
  Omu_Program	*_sys;
  Omu_SVec	*_xc_ptr;
  Omu_DepVec 	*_Fc_ptr;

  // variables for DASPK
  int		_mu;	// upper semi-bandwidth
  int		_ml;	// lower semi-bandwidth
  int    	_lwp;
  int     	_liwp;
  int    	_lwp_basic;

  VECP		_y;
  VECP		_yprime;
  IVECP		_info;
  VECP		_rwork;
  IVECP		_iwork;
  VECP		_rpar;
  IVECP		_ipar;
  VECP		_senpar;

  // arguments for low level _sys->continuous callback
  Omu_Vec	_uc;
  Omu_SVec	_xcp;
  Omu_SVec	_xc_jac;
  Omu_SVec	_xcp_jac;
  MATP		_Yx;
  MATP		_Yu;
};  

#endif

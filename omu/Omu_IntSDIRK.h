/*
 * Omu_IntSDIRK.h --
 *   -- class for integrating an Ode over a stage
 *   -- using singly diagonally implicit Runge Kutta method
 *
 * hl, 03/08/00
 */

/*
    Copyright (C) 2000  Hartmut Linke

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

#ifndef Omu_IntSDIRK_H
#define Omu_IntSDIRK_H

#include "Omu_Integrator.h"

/**
 * Singly diagonally implizit Runge Kutta method for (stiff)
 * ordinary differential equations (ODEs). The implementation
 * is based on Hairer/Wanner: "Solving ordinary differential
 * equations II".
 */

//--------------------------------------------------------------------------
class Omu_IntSDIRK: public Omu_Integrator {

 public:

  Omu_IntSDIRK(); 	///< constructor
  ~Omu_IntSDIRK(); 	///< destructor

  /**
   * @name Implementation of predefined methods.
   * @see Omu_Integrator
   */

  //@{ 

  char *name() {return "SDIRK";}

  // interface routine from Omu_Integrator
  virtual void init(int k,
		    const Omu_StateVec &xc, const Omu_Vec &q,
		    const Omu_DependentVec &Fc, bool sa);

  virtual void solve(int kk, double tstart, double tend,
		     Omu_StateVec &xc, Omu_StateVec &dxc, Omu_Vec &q,
		     Omu_DependentVec &Fc);

  //@}

 protected:


  /**
   * User given semi-bandwidth of Jacobian (default: -1).
   * If a value >= 0 is specified, then it is used instead
   * of the automatic detection.
   */
  int		_jac_sbw;

  /**
   * User specification to allow banded solver (default: false).
   * Banded solvers are used if _banded is true, _jac_sbw < 0,
   * and if the automatic detection indicates that the problem
   * can be solved more efficiently in this way.
   */
  bool		_banded;
  bool		_banded_solver;	  ///< internal flag for using banded solver

  /**
   * User specification to allow sparse solver (default: false).
   */
  bool		_sparse_solver;

  private:

  void init_yprime(int , double ,
		   const Omu_StateVec &,const Omu_Vec &, VECP );
  void jac(int ,double ,const VECP ,const VECP , const VECP , VECP );
  void res(double ,const VECP ,const VECP , const VECP , VECP);
  void resize();
  void init_method2();
  void init_method4();
  void init_method5();
  void init_yprime_pred(MATP );
  void solve_stage(int , VECP , VECP );
  void solve_final(VECP , VECP );
  void sensitivity_stage(int , const VECP , const VECP , const VECP );
  void sensitivity_update( );
  double  error_check();

  // backing store sys and vector of dependent variables for callbacks
  Omu_StateVec	    *_xc_ptr;
  Omu_DependentVec  *_Fc_ptr;

  bool          _ode;
  bool          _recalc_jac;
  bool          _stiffly_accurate;
  bool          _sens_adolc;

  int           _output;

  int		_mu;	     // upper semi-bandwidth
  int		_ml;	     // lower semi-bandwidth

  double        _t;

  int           _irk_stages;
  int           _maxiters;
  int           _modnewtonsteps;
  int           _newtonsteps;
  int           _nsteps;

  MATP          _irk_jac;
  BAND          *_irk_jac_bd;
  SPMATP        _irk_jac_sp;

  MATP          _irk_jacf;

  PERM          *_ppivot;

  IVECP         _x_algebraic;

  VECP          _irk_delta;
  VECP          _irk_res;
  VECP          _irk_y;
  VECP          _irk_yprime;

  VECP          _senpar;
  VECP          _y;
  VECP          _y0;
  VECP          _yprime0;
  VECP          _yprime1;
  VECP          _yprime2;

  MATP          _yprime;

  // parameters of irk method

  bool          _init_method;
  int           _method;

  double        _gamma0;
  VECP          _b;
  VECP          _b_err;
  VECP          _c;
  MATP          _a;
  MATP          _a_pred;

  // stepsize selection
  double        _hinit;
  double        _h_new;
  double        _h;
  VECP          _err;

  // newton parameters
  double        _eta;
  double        _kappa;

  // vectors and matrices for low level _sys->continuous callback
  Omu_Vec      	_q;
  Omu_SVec	_xcp;

  VECP		_Fch;
  MATP		_Fcxh;

  // sensitivity equations
  MATP          _Sh;
  MATP          _Sh1;
  MATP          _Sh2;
  MATP          _Sh3;
  PERM          *_Spivot;

  VECP          _vh;
  SPMATP        _smh;
  SPMATP        _smh1;

};

#endif








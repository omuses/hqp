/*
 * Omu_Integrator.h --
 *   -- abstract base class
 *   -- integrate system equations over a discrete stage
 *
 * rf, 10/2/96
 */

/*
    Copyright (C) 1997--2002  Ruediger Franke

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

#ifndef Omu_Integrator_H
#define Omu_Integrator_H

#include <If_List.h>
#include <If_Class.h>

#include "Omu_Program.h"
#include "Omu_Vars.h"
#include "Omu_Deps.h"

IF_BASE_DECLARE(Omu_Integrator);

/**
 * Interface for a solver for differential (algebraic) equations. The solver
 * is employed internally to numerically integrate continuous-time model
 * equations over sample periods and discrete-time stages.
 */
class Omu_Integrator {

 public:

  Omu_Integrator();		///< constructor
  virtual ~Omu_Integrator();	///< destructor

  /** Setup stages for integrator. This method is called once
      after Omu_Program's setup_stages has been performed. */
  virtual void setup_stages(const Omu_Program *sys) {}

  /** Initialize solution of differential equations for stage k.
      The numbers of variables and equation structure (given in F)
      shall be the same for all sample periods of a stage. */
  virtual void init_stage(int k,
			  const Omu_States &x, const Omu_Vector &u,
			  const Omu_DepVec &F, bool sa = false);

  /** Initialize solution for one sample period kk. */
  virtual void init_sample(int kk, double tstart, double tend);

  /** Solve differential equations over sample period kk from tstart to tend.
      The differential-algebraic equations can be evaluated by calling 
      Omu_Program::continuous for sys. cF contains pre-initialized data
      structures for resisuals and Jacobians. The initial solution and the
      result are stored in cx, including values and sensitivities
      (dim(cx)=_nxt, dim(cx.Sx)=_nxt._nx, dim(cx.Su)=_nxt._nu).
      Note that the first _nd elements of cx are discrete-time states
      that need to be treated like control parameters u by an integrator.
      Furthermore, the last _nv elements of cx are expansion variables that
      are not seen by the optimizer and accordingly no sensitivities are
      required with respect to them.
      The default implementation of solve calls the depreciated version. */
  virtual void solve(int kk, double tstart, double tend,
		     const Omu_States &x, const Omu_Vector &u,
		     Omu_Program *sys, Omu_DepVec &cF, Omu_SVec &cx);

  /** Name of a specific integrator. */
  virtual char *name() = 0;

 protected:

  /**
   * @name Depreciated methods.
   * These methods are called by the default implementations of hte current
   * versions.
   */
  //@{

  /** Initialize solution of differential equations for stage k. */
  virtual void init_stage(int k,
			  const Omu_States &x, const Omu_Vector &u,
			  bool sa = false);

  /** Solve differential equations over sample period. */
  virtual void solve(int kk, double tstart, double tend,
		     const Omu_States &x, const Omu_Vector &u,
		     Omu_Program *sys, VECP xt,
		     MATP Sx = MNULL, MATP Su = MNULL) {}

  //@}

  int		_kk; 	///< index of current sample period
  int		_nxt;	///< number of continuous-time variables
  int		_nd;	///< number of discrete-time states
  int		_nv;	///< number of expansion variables
  int		_nx;	///< number of states treated by optimizer (nx=nxt-nv)
  int		_nu;	///< number of control parameters

  /**
   * Indicate if sensitivity analysis is required together with
   * solution of differential equations.
   */
  bool		_sa;

  /**
   * Boolean to indicate if sensitivities should be considered
   * in error test (default: false).
   */
  bool		_serr;

  /**
   * Number of continuous-time states during integration (nxt-nd: 
   * total number minus discrete-time states).
   */
  int		_n;

  /**
   * Number of parameters during integration (nd+nu: 
   * discrete-time states and control parameters).
   */
  int		_m;

  If_List	_ifList; ///< container for interface elements

  /**
   * Constant stepsize to be taken by integator.
   * For stepsize>0, an integrator takes ceil((tend-tstart)/stepsize)
   * equidistant steps in the sampling interval [tstart,tend].
   * Default: 0.0, i.e. integrator chooses own stepsize.
   */
  double	_stepsize;

  /**
   * Relative integration error (only for variable stepsize).
   */
  double        _rtol;

  /**
   * Absolute integration error (only for variable stepsize).
   */
  double        _atol;

  int		_res_evals;	///< number of residual evaluations
  int		_sen_evals;	///< number of sensitivity evaluations
  int		_jac_evals;	///< number of Jacobian evaluations
};  

#endif


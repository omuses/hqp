/**
 * @file Omu_Integrator.h 
 *    interface for integrating continuous-time system equations
 *
 * rf, 10/2/96
 */

/*
    Copyright (C) 1997--2014  Ruediger Franke

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

/**
 * Interface for a solver for differential (algebraic) equations. The solver
 * is employed internally to numerically integrate continuous-time model
 * equations over sample periods and discrete-time stages.
 */
class Omu_Integrator {

 public:

  Omu_Integrator();		///< constructor
  virtual ~Omu_Integrator();	///< destructor

  /**
   * @name Low-level interface for Omuses employing an integrator.
   * The default implementations call the high-level interface.
   */
  //@{
  /** Setup stages for integrator. This method is called once
      after Omu_Program's setup_stages has been performed.
      @param sys pointer to problem definition
  */
  virtual void setup_stages(const Omu_Program *sys);

  /** Setup struct for one stage k. This method is called once during
      problem setup.
      Note: The numbers of variables and equation structure (given in Ft)
      shall be the same for all sample periods of a stage.
      @param k stage
      @param x state variables for optimization, including bounds etc.
      @param u control parameters for optimization
      @param Ft pre-initialized data structures for resisuals and Jacobians
  */
  virtual void setup_struct(int k,
			    const Omu_VariableVec &x, const Omu_VariableVec &u,
			    const Omu_DependentVec &Ft);

  /** Initialize solution of differential equations for stage k.
      Note: The numbers of variables and equation structure (given in F)
      shall be the same for all sample periods of a stage.
      @param k stage
      @param x state variables for optimization, including bounds etc.
      @param u control parameters for optimization
      @param Ft pre-initialized data structures for resisuals and Jacobians
      @param sa true if sensitivity analysis shall be performed
  */
  virtual void init_stage(int k,
			  const Omu_VariableVec &x, const Omu_VariableVec &u,
			  const Omu_DependentVec &Ft, bool sa = false);

  /** Solve differential equations over one sample period.
      The differential-algebraic equations can be evaluated by calling 
      Omu_Program::continuous for sys.
      Note that the first _nd elements of xt are discrete-time states
      that need to be treated like control parameters u by an integrator.
      Furthermore, the last _nv elements of xt are expansion variables that
      are not seen by the optimizer and accordingly no sensitivities are
      required with respect to them.
      The default implementation of solve calls the high-level version.
      @param kk sample period
      @param tstart start time for integration
      @param tend end time for integration
      @param x state variables for optimization, including bounds etc.
      @param u control parameters for optimization
      @param sys problem definition for calling Omu_Program::continuous
      @param Ft pre-initialized data structures for resisuals and Jacobians
      @param xt start values for integration at tstart, dim(xt)=_nxt;
        initial sensitivity matrices xt.Sx, dim(xt.Sx)=_nxt._nx, and 
        xt.Su, dim(xt.Su)=_nxt._nu, if _sa is true
      @retval xt result at tend; xt.Sx and xt.Su at tend if _sa is true
      @throw E_CONV indicate that the integrator failed to converge; 
        may be caught during optimization step length test
  */
  virtual void solve(int kk, double tstart, double tend,
		     const Omu_VariableVec &x, const Omu_VariableVec &u,
		     Omu_Program *sys, Omu_DependentVec &Ft, Omu_StateVec &xt);
  //@}

  /**
   * @name High-level interface for derived classes implementing an integrator.
   *  This interface hides details of Omuses, like the classification of
   *  variables as discrete-time states, expansion states, and
   *  control parameters. Instead it treats continuous-time DAE states
   *  and sensitivity parameters.
   */
  //@{

  /** Initialize integrator, like adaptation of array dimensions and
      test of pre-conditions, like explicit ODE. An explicit ODE is present
      if Fc.Jdx.is_scalar() and Fc.Jdx[0][0] != 0.0. An algebraic state
      j, 0<=j<_n can be identified with Fc.Jdx.is_zero_column[j].
      Note that the init method should not find consistent initial conditions
      of a DAE as parameters may jump between sample periods kk afterwards.
      @param k stage
      @param xc continuous-time states, dim(xc)=_n
      @param q sensitivity parameters, dim(q)=_nq
      @param Fc pre-initialized data structures for evaluating residuals,
        dim(Fc)=_n, including structurally analyzed Jacobians,
        dim(Fc.Jx)=_n._n, dim(Fc.Jxp)=_n._n, dim(Fc.Jq)=_n._nq
      @param sa indicates if sensitivity analysis is required 
  */
  virtual void init(int k,
		    const Omu_StateVec &xc, const Omu_Vec &q,
		    const Omu_DependentVec &Fc, bool sa) {}

  /** Solve differential equations over one sample period. Additionally
      solve sensitivity equations if _sa is true.
      The default implementation calls the depreciated solve method.
      @param kk sample period
      @param tstart start time for integration
      @param tend end time for integration
      @param xc continuous-time states at tstart, dim(xc)=_n; sensitivities
        xc.Sq at tstart, dim(xc.Sq)=_n._nq, if _sa is true
      @param dxc pre-allocated time derivatives of continuous-time states,
        dim(dxc)=_n
      @param q sensitivity parameters, dim(q)=_nq
      @param Fc pre-allocated argument for residual calls, dim(Fc)=_n,
        including structurally analyzed Jacobians Fc.Jx, Fc.Jdx, Fc.Jq
      @retval xc continuous-time states at tend;
        sensitivities xc.Sq at tend if _sa is true
      @throw E_CONV indicate that the integrator failed to converge; 
        may be caught during optimization step length test
  */
  virtual void solve(int kk, double tstart, double tend,
		     Omu_StateVec &xc, Omu_StateVec &dxc, Omu_Vec &q,
		     Omu_DependentVec &Fc);

  /** Callback to evaluate residuals of continuous-time equations.
      @param kk sample period
      @param t time for evaluation
      @param xc continuous-time states, dim(xc)=_n
      @param dxc derivatives of continuous-time states wrt time, dim(dxc)=_n
      @param q sensitivity parameters, dim(q)=_nq
      @retval Fc residuals for continuous-time equations, dim(Fc)=_n;
        Jacobians Fc.Jx, Fc.Jdx, Fc.Jq if Fc.is_required_J() is true
   */
  virtual void residual(int kk, double t,
			const Omu_StateVec &xc, const Omu_StateVec &dxc,
			const Omu_Vec &q, Omu_DependentVec &Fc);
  //@}

  /**
   * @name Member access methods
   */
  //@{

  /// get number of stages for which integrator has been set up
  int		K() const {return _K;}

  //@}

  /** Name of a specific integrator. */
  virtual const char *name() = 0;

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
			  bool sa = false) {}

  /** Solve differential equations over sample period. */
  virtual void solve(int kk, double tstart, double tend,
		     const Omu_States &x, const Omu_Vector &u,
		     Omu_Program *sys, VECP xt,
		     MATP Sx = MNULL, MATP Su = MNULL) {}

  //@}

  /**
   * Pointer to problem definition that is valid during solve
   */
  Omu_Program 	*_sys;

  int		_K; 	///< number of stages for which integrator was set up
  int		_k; 	///< index of current stage
  int		_kk; 	///< index of current sample period
  int		_nxt;	///< total number of DAE states
  int		_nd;	///< number of discrete-time states
  int		_nv;	///< number of expansion variables
  int		_na;	///< number of algebraic states (no time derivative)
  int		_nx;	///< number of states treated by optimizer (nx=nxt-nv)
  int		_nu;	///< number of control parameters
  int		_nq;	///< number of sensitivity parameters (nq=nx+nu)

  /**
   * Indicate if sensitivity analysis is required together with
   * solution of differential equations.
   */
  bool		_sa;

  /**
   * Boolean to indicate if sensitivities should be considered
   * for error tolerance (default: false).
   */
  bool		_serr;

  /**
   * Number of continuous-time states during integration (nxt-nd: 
   * total number minus discrete-time states).
   */
  int		_n;

  /**
   * Number of control parameters during integration (nd+nu: 
   * discrete-time states and actual control parameters).
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
   * Minimal stepsize to be taken by integator.
   * Default: 0.0, i.e. no limitation.
   */
  double	_min_stepsize;

  /**
   * Relative integration error tolerance (only for variable stepsize).
   */
  double        _rtol;

  /**
   * Absolute integration error tolerance (only for variable stepsize).
   */
  double        _atol;

  int		_res_evals;	///< number of residual evaluations
  int		_sen_evals;	///< number of sensitivity evaluations
  int		_jac_evals;	///< number of Jacobian evaluations

 private:

  /** initialize _nx, _nu, etc. for stage k */
  void init_dims(int k,
		 const Omu_VariableVec &x, const Omu_VariableVec &u,
		 const Omu_DependentVec &Ft);

  /**
   * Pointers to data available in solve() and required by residual()
   */
  //@{
  const Omu_VariableVec *_x_ptr;
  const Omu_VariableVec *_u_ptr;
  Omu_DependentVec 	*_Ft_ptr;
  Omu_StateVec 		*_xt_ptr;
  //@}

  Omu_Vec 	_ut; 	///< control parameterstime for sys->continuous
  Omu_SVec 	_dxt; 	///< time derivatives of states for sys->continuous

  Omu_SVec 	_xc; 	///< continuous-time states for integration
  Omu_SVec 	_dxc; 	///< time derivatives of continuous-time states
  Omu_Vec	_q; 	///< sensitivity parameters for integration
  Omu_DepVec	*_Fcs; 	///< residual arguments for all stages
};  

OMU_API IF_BASE_DECLARE(Omu_Integrator);

#endif

/*
 * Omu_Program.h --
 *   -- declare an interface for a multistage optimization problem
 *      described by differential and algebraic equations
 *
 * rf, 16/1/97
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

#ifndef Omu_Program_H
#define Omu_Program_H

#include <adouble.h>
#include <Meschach.h>
#include "Omu_Variables.h"
#include "Omu_Dependents.h"

#include <If_Class.h>

IF_BASE_DECLARE(Omu_Program);

/**
 * Interface for a multi-stage optimization problem. The problem may be
 * defined in a derived class by overloading the appropriate methods.
 * The problem may contain discrete-time and continuous-time equations.
 * It may be formulated using ADOL-C for automatic differentiation
 * (defining the high-level methods) or by providing structural information
 * and derivatives manually (defining the low-level methods).
 */
class Omu_Program {

 public:

  Omu_Program();		///< constructor
  virtual ~Omu_Program();	///< destructor

  /*
   * High-level methods to be primarily implemented by a derived class.
   */

  /**
   * Setup stages and time sampling. The default implementation sets up
   * an optimization problem without stages.
   */
  virtual void setup_stages(IVECP ks, VECP ts);

  /**
   * Allocate states x, control parameters u, and constraints c for stage k.
   * Additionally setup bounds and initial values.
   */
  virtual void setup(int k,
		     Omu_VariableVec &x, Omu_VariableVec &u,
		     Omu_VariableVec &c) = 0;

  /**
   * Initialize problem variables using simulation. 
   * A problem specification may set state variables x and/or control
   * parameters u prior to the evaluation of stage k.
   * The default version causes an initial value simulation based
   * on initial values given in setup.
   */
  virtual void init_simulation(int k,
			       Omu_VariableVec &x, Omu_VariableVec &u);

  /**
   * High-level update for the optimization criterion, constraints, and
   * discrete state equations.
   */
  virtual void update(int kk, 
		      const adoublev &x, const adoublev &u,
		      adoublev &f, adouble &f0, adoublev &c);

  /**
   * High-level consistic (consistent initial conditions) routine
   * for the initialization of continuous-time states from
   * optimization variables x and u.
   * The routine is called via the default low-level version
   * before the first call to continuous in each sample period.
   * The default implementation copies x to xt.
   */
  virtual void consistic(int kk, double t,
			 const adoublev &x, const adoublev &u,
			 adoublev &xt);

  /**
   * High-level continuous routine for specifying differential equations
   * of the form F(x,u,\dot{x})=0.
   * The default implementation is empty to indicate that no differential
   * equations are used.
   */
  virtual void continuous(int kk, double t,
			  const adoublev &x, const adoublev &u,
			  const adoublev &xp, adoublev &F);

  /**
   * @name Low-level methods
   * Low-level methods are primarily called from within Omuses to evaluate
   * the problem specification. Default implementations are provided
   * that call the corresponding high-level methods and use ADOL-C
   * to obtain required structural information and derivatives.
   * A derived class should only overload these low-level methods
   * in special cases, e.g. to achieve a speedup.
   */
  //@{

  /**
   * Setup stages and time sampling. The default implementation calls
   * setup_stages(IVECP ks, VECP ts).
   */
  virtual void setup_stages();

  /**
   * Setup sparsity patterns of Jacobians. 
   * Note: The sparsity pattern must be the same for all sample periods of
   * a stage, i.e. this method has to setup a superset of all sparsity
   * patterns of individual sample periods.
   */
  virtual void setup_struct(int k,
			    const Omu_VariableVec &x, const Omu_VariableVec &u,
			    Omu_DependentVec &xt, Omu_DependentVec &F,
			    Omu_DependentVec &f,
			    Omu_Dependent &f0, Omu_DependentVec &c);

  /**
   * Low-level consistic (consistent initial conditions) routine
   * for the user specification of derivative information.
   */
  virtual void consistic(int kk, double t,
			 const Omu_StateVec &x, const Omu_Vec &u,
			 Omu_DependentVec &xt);

  /**
   * Low-level continuous routine.
   */
  virtual void continuous(int kk, double t,
			  const Omu_StateVec &x, const Omu_Vec &u,
			  const Omu_StateVec &xp, Omu_DependentVec &F);

  /**
   * Low-level update. 
   */
  virtual void update(int kk, 
		      const Omu_StateVec &x, const Omu_Vec &u,
		      const Omu_StateVec &xf,
		      Omu_DependentVec &f, Omu_Dependent &f0,
		      Omu_DependentVec  &c);

  /**
   * This routine is intended for a calling module,
   * which might want to avoid calling the low-level continuous
   * if it is not overridden (e.g. as it has then a faster sensitivity
   * calculation).
   * @return true if low-level continuous is overloaded
   * (currently low-level continuous must have been called once before 
   *  false may be returned)
   */
  bool has_low_level_continuous();
  //@}

  /**
   * @name Member access routines
   */
  //@{

  int		K() {return _K;}	///< get number of stages
  void		set_K(int K) {_K = K;}	///< set number of stages
  int		KK() {return _KK;}	///< get number of sample periods
  void		set_KK(int KK) {_KK = KK;}///< set number of sample periods

  /// get vector of start indices for sample periods in each stage
  const IVECP 	ks() {return _ks;}

  /// get vector of start time points in each sample period
  const VECP 	ts() {return _ts;}

  /// get start index of sample periods in stage k
  int	 	ks(int k) {return _ks[k];}

  /// get start time point of sample period kk
  Real 		ts(int kk) {return _ts[kk];}

  //@}

  /**
   * Textual name of a problem definition.
   */
  virtual char *name() = 0;

  /**
   * @name Depreciated methods
   */
  //@{

  /**
   * Old low-level consistic (consistent initial conditions) routine
   * for the user specification of derivative information.
   * The default implementation calls the high-level version
   * and uses ADOL-C for the determination of xtx and xtu.
   * If a calling routine passes no xtx and xtu (or null pointers),
   * then no Jacobians are calculated at all.
   * If a calling routine passes xtx and xtu, then the matrices are
   * initialized with zeros.
   */
  virtual void consistic(int kk, double t,
			 const Omu_Vector &x, const Omu_Vector &u,
			 VECP xt, MATP xtx = MNULL, MATP xtu = MNULL);

  /**
   * Old low-level continuous routine,
   * which may be implemented in addition to high-level continuous
   * by a derived class to provide Jacobians for differential equations.
   * The default implementation calls the high-level continuous()
   * and employs ADOL-C for automatic Jacobian calculation.
   * If a calling routine passes no Fxp (or a null pointer),
   * then only Fx and Fu are calculated.
   * If a calling routine passes no Fx, Fu, and Fxp (or null pointers),
   * then no Jacobians are calculated at all.
   */
  virtual void continuous(int kk, double t,
			  const VECP x, const VECP u,
			  const VECP xp, VECP F,
			  MATP Fx = MNULL, MATP Fu = MNULL,
			  MATP Fxp = MNULL);

  //@}

 protected:

  If_List	_ifList; ///< container for interface elements

  //
  // predefined basic variables
  //

  int 	_K;	///< number of stages
  int 	_KK;	///< number of sample periods over all stages

  /**
   * Service routine that can be called from setup_stages for
   * initalizing ks and ts.
   */
  void 	stages_alloc(IVECP ks, VECP ts, int K, int sps,
		     double t0 = 0.0, double tf = 1.0);

  /**
   * @name Low-level gradient routines
   * These routines can be used to obtain Jacobians from low-level
   * methods employing numerical differentiation.
   * They can be called from the corresponding
   * low-level method after the values have been calculated.
   * The default versions implement finite forward differences.
   */
  //@{
   
  /**
   * Low-level consistic routine for obtaining gradients.
   */
  virtual void consistic_grds(int kk, double t,
			      const Omu_StateVec &x, const Omu_Vec &u,
			      Omu_DependentVec &xt);

  /**
   * Low-level continuous routine for obtaining gradients.
   */
  virtual void continuous_grds(int kk, double t,
			       const Omu_StateVec &x, const Omu_Vec &u,
			       const Omu_StateVec &xp, Omu_DependentVec &F);

  /**
   * Low-level update routine for obtaining gradients.
   */
  virtual void update_grds(int kk, 
			   const Omu_StateVec &x, const Omu_Vec &u,
			   const Omu_StateVec &xf,
			   Omu_DependentVec &f, Omu_Dependent &f0,
			   Omu_DependentVec  &c);
  //@}

private:
  IVECP _ks;	// starting indices of stages
  VECP 	_ts;	// time steps

  // allocation methods and memory for ADOL-C
  void		ad_alloc(int ndep, int nindep);
  void		ad_realloc(int ndep, int nindep);
  void		ad_free();
  int		_max_ndep;
  int		_max_nindep;
  MATP 		_U;
  MATP 		_Z;
  double	***_Z3;
  short	  	**_nz;

  bool 		_has_low_level_continuous;
};  

#endif

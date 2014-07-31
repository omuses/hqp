/**
 * @file Hqp_Omuses.h 
 *   multi stage optimal control problems described by DAE's
 *
 * rf, 7/27/96
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

#ifndef Hqp_Omuses_H
#define Hqp_Omuses_H

#include <Hqp_Docp.h>

#include <If_List.h>

class Omu_Program;
class Omu_Integrator;
class Omu_VarVec;
class Omu_SVarVec;
class Omu_SVec;
class Omu_DepVec;
class Omu_Dep;

/**
 * Extend DOCP interface with treatment of continuous-time differential
 * algebraic equations (DAE) and multiple sample periods per stage.
 */
class Hqp_Omuses: public Hqp_Docp {

 public:
  Hqp_Omuses();
  ~Hqp_Omuses();

  /**
   * @name Member access methods.
   */
  //@{
  /// optimization program
  Omu_Program *prg() const {return _prg;}
  /// Set optimization program.
  /// Remind to delete an existing prg before creating a new one!
  void set_prg(Omu_Program *prg) {_prg = prg;}

  /// solver for continuous-time differential equations
  Omu_Integrator *integrator() const {return _integrator;}
  /// Set solver for continuous-time differential equations.
  /// Remind to delete an existing integrator before creating a new one!
  void set_integrator(Omu_Integrator *integrator) {_integrator = integrator;}

  /// flag about use of automatic differentiation
  bool ad() const {return _ad;}
  /// Set flag about use of automatic differentiation.
  void set_ad(bool val) {_ad = val;}

  /// scaling of optimization criterion
  double fscale() const {return _fscale;}
  /// Set scaling of optimization criterion.
  void set_fscale(double val) {_fscale = val;}

  /// setup method for structure of Hessian of Lagrangian.
  /// 0: full Hessian, 1: setup based on structure information
  int hela_setup() const {return _hela_setup;}
  /// Set setup method for structure of Hessian of Lagrangian.
  void set_hela_setup(int val) {_hela_setup = val;}

  //@}

  const char *name() {return "Omuses";} ///< generic name Omuses

  /** setup stages of optimization program */
  void setup_stages();

 protected:
  If_List		_ifList; 	///< interface elements
  Omu_Program		*_prg; 		///< optimization program
  Omu_Integrator	*_integrator; 	///< integrator

  bool 	_stages_ok;	///< setup_stages() was called separately
  MATP 	_IS;		///< help matrix
  bool 	_ad;		///< flag about use of automatic differentiation
  double _fscale;	///< scaling of optimization criterion
  int 	_hela_setup;	///< setup structure of Hessian of Lagrangian

  /** Setup integrator for all stages.
      This is required if the integrator is exchanged after problem setup. */
  void resetup_integrator();

  /**
   * @name Implement interface of Hqp_Docp_stub
   */
  //@{
  void setup_horizon(int &k0, int &kf);

  void setup_vars(int k,
		  VECP x, VECP x_min, VECP x_max, IVECP x_int,
		  VECP u, VECP u_min, VECP u_max, IVECP u_int,
		  VECP c, VECP c_min, VECP c_max);

  void setup_struct(int k, const VECP x, const VECP u,
		    MATP fx, MATP fu, IVECP f_lin,
		    VECP f0x, VECP f0u, int &f0_lin,
		    MATP cx, MATP cu, IVECP c_lin,
		    MATP Lxx, MATP Luu, MATP Lxu);

  void init_simulation(int k,
		       VECP x, VECP u);

  void update_vals(int k, const VECP x, const VECP u,
		   VECP f, Real &f0, VECP c);

  void update_stage(int k, const VECP x, const VECP u,
		    VECP f, Real &f0, VECP c,
		    MATP fx, MATP fu, VECP f0x, VECP f0u,
		    MATP cx, MATP cu,
		    const VECP rf, const VECP rc,
		    MATP Lxx, MATP Luu, MATP Lxu);
  //@}

  Omu_SVarVec 	*_xs;	///< state information from problem setup
  Omu_VarVec	*_us;	///< control information from problem setup
  Omu_VarVec	*_css;	///< constraint information from problem setup

  Omu_SVec 	*_x0s;	///< initial states before integration
  Omu_SVec 	*_xfs;	///< final states after integration

  Omu_DepVec 	*_xts;	///< continuous time states from consistic
  Omu_DepVec 	*_Fs;	///< continuous time model equations from continuous
  Omu_DepVec 	*_fs;	///< discrete time states from update
  Omu_Dep 	*_f0s;  ///< objective from update
  Omu_DepVec 	*_cs; 	///< constraints from update
};  

#endif


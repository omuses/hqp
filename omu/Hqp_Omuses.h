/*
 * Hqp_Omuses.h --
 *   -- multi stage optimal control problems described by DAE's
 *
 * rf, 7/27/96
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

#ifndef Hqp_Omuses_H
#define Hqp_Omuses_H

#include "Hqp_DocpStub.h"

#include <adouble.h>
#include <If_List.h>
#include <If_Command.h>

class Omu_Program;
class Omu_Integrator;
class Omu_VarVec;
class Omu_DynVarVec;
class Omu_SVec;
class Omu_DepVec;
class Omu_Dep;

//--------------------------------------------------------------------------
class Hqp_Omuses: public Hqp_DocpStub {

 public:

  Hqp_Omuses();
  ~Hqp_Omuses();
  char *name() {return "Omuses";}

  int setup_stages(IF_DEF_ARGS);

 protected:

  If_List		_ifList;
  Omu_Program		*_prg;
  Omu_Integrator	*_integrator;
  void 			*_integrator_setup;

  bool _stages_ok;	// setup_stages() was called separately
  VECP _xt;		// initial states of a sample period
  MATP _xtxk;		// dxt/dxk
  MATP _xtuk;		// dxt/duk
  MATP _Sx;		// sensitivity for x after integration
  MATP _Sxk;		// sensitivity for x before integration
  MATP _IS;		// help matrix
  MATP _Su;		// sensitivity for u after integration
  MATP _Suk;		// sensitivity for u before integration

  VECP _fk;		// call argument for update

  MATP _fkxk;		// Jacobian delivered by update
  MATP _fkuk;		// Jacobian delivered by update
  MATP _fkfk;		// Jacobian delivered by update
  VECP _f0kxk;		// Jacobian delivered by update
  VECP _f0kuk;		// Jacobian delivered by update
  VECP _f0kfk;		// Jacobian delivered by update
  MATP _ckxk;		// Jacobian delivered by update
  MATP _ckuk;		// Jacobian delivered by update
  MATP _ckfk;		// Jacobian delivered by update

  bool _ad;		// flag about use of automatic differentiation
  double _fscale;	// scaling of the criterion

  // define interface of Hqp_Docp_stub
  void setup_horizon(int &k0, int &kf);

  void setup_vars(int k,
		  VECP x, VECP xmin, VECP xmax,
		  VECP u, VECP umin, VECP umax,
		  VECP c, VECP cmin, VECP cmax);

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

  // further methods
  void init_vars(int k,
		 VECP x, VECP u);

  void obtain_structure(int k,
			Omu_DynVarVec &xk, const Omu_VarVec &uk);

  Omu_DynVarVec	*_xs;	// state information from problem setup
  Omu_VarVec	*_us;	// control information from problem setup
  Omu_VarVec	*_css;	// constraint information from problem setup

  Omu_SVec 	*_x0s;	// initial states before integration
  Omu_SVec 	*_xfs;	// final states after integration

  Omu_DepVec 	*_xts;	// continuous time states from consistic
  Omu_DepVec 	*_Fs;	// continuous time model equations from continuous
  Omu_DepVec 	*_fs;	// discrete time states from update
  Omu_Dep 	*_f0s;  // objective from update
  Omu_DepVec 	*_cs; 	// constraints from update

  // variables for ADOL-C

  void		ad_alloc(int ndep, int nindep, int npar);
  void		ad_realloc(int ndep, int nindep, int npar);
  void		ad_free();
  int		_max_ndep;
  int		_max_nindep;
  int		_max_npar;

  VECP		_x;
  VECP		_y;
  MATP		_X;
  MATP		_Y;
  MATP		_U;
  double	***_Z3;
  short	  	**_nz;
};  

#endif


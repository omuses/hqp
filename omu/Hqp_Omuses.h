/*
 * Hqp_Omuses.h --
 *   -- multi stage optimal control problems described by DAE's
 *
 * rf, 7/27/96
 */

/*
    Copyright (C) 1997--2000  Ruediger Franke

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

#include "Hqp_Docp_stub.h"

#include <adouble.h>
#include <If_List.h>
#include <If_Command.h>

class Omu_Program;
class Omu_Integrator;
class Omu_Vars;
class Omu_Vector;
class Omu_States;

//--------------------------------------------------------------------------
class Hqp_Omuses: public Hqp_Docp_stub {

 public:

  Hqp_Omuses();
  ~Hqp_Omuses();
  char *name() {return "Omuses";}

  int setup_stages(IF_DEF_ARGS);

 protected:

  If_List		_ifList;
  Omu_Program		*_prg;
  Omu_Integrator	*_integrator;

  bool _stages_ok;	// setup_stages() was called separately
  VECP _xt;		// initial states of a sample period
  MATP _xtxk;		// dxt/dxk
  MATP _xtuk;		// dxt/duk
  MATP _Sx;		// sensitivity for x
  MATP _IS;		// help matrix
  MATP _Su;		// sensitivity for u
  bool _ad;		// flag about use of automatic differentiation
  double _fscale;	// scaling of the criterion

  // define interface of Hqp_Docp_stub
  void setup_horizon(int &k0, int &kf);

  void setup_vars(int k,
		  VECP x, VECP xmin, VECP xmax,
		  VECP u, VECP umin, VECP umax,
		  VECP c, VECP cmin, VECP cmax);

  void setup_struct(int k,
		    VECP f0x, VECP f0u, int &f0_lin,
		    MATP fx, MATP fu, IVECP f_lin,
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
			Omu_States &xk, const Omu_Vector &uk);

  Omu_States	*_xs;	// state information from problem setup
  Omu_Vars	*_us;	// control information from problem setup
  Omu_Vars	*_cs;	// constraint information from problem setup

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


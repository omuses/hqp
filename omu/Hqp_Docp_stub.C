/*
 * Hqp_Docp_stub.C -- class definition
 *
 * rf, 1/11/00
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

#include "Hqp_Docp_stub.h"

//-------------------------------------------------------------------------
static void setup_horizon(void *clientdata, int &k0, int &kf)
{
  ((Hqp_Docp_stub *)clientdata)->setup_horizon(k0, kf);
}

//-------------------------------------------------------------------------
static void setup_vars(void *clientdata, int k,
		       VECP x, VECP xmin, VECP xmax,
		       VECP u, VECP umin, VECP umax,
		       VECP c, VECP cmin, VECP cmax)
{
  ((Hqp_Docp_stub *)clientdata)->setup_vars(k, x, xmin, xmax,
					    u, umin, umax,
					    c, cmin, cmax);
}

//-------------------------------------------------------------------------
static void setup_struct(void *clientdata, int k,
			 VECP f0x, VECP f0u, int &f0_lin,
			 MATP fx, MATP fu, IVECP f_lin,
			 MATP cx, MATP cu, IVECP c_lin,
			 MATP Lxx, MATP Luu, MATP Lxu)
{
  ((Hqp_Docp_stub *)clientdata)->setup_struct(k,
					      f0x, f0u, f0_lin,
					      fx, fu, f_lin,
					      cx, cu, c_lin,
					      Lxx, Luu, Lxu);
}

//-------------------------------------------------------------------------
static void init_simulation(void *clientdata, int k,
			    VECP x, VECP u)
{
  ((Hqp_Docp_stub *)clientdata)->init_simulation(k, x, u);
}

//-------------------------------------------------------------------------
static void update_vals(void *clientdata, int k,
			const VECP x, const VECP u,
			VECP f, Real &f0, VECP c)
{
  ((Hqp_Docp_stub *)clientdata)->update_vals(k, x, u, f, f0, c);
}

//-------------------------------------------------------------------------
static void update_stage(void *clientdata, int k,
			 const VECP x, const VECP u,
			 VECP f, Real &f0, VECP c,
			 MATP fx, MATP fu,
			 VECP f0x, VECP f0u,
			 MATP cx, MATP cu,
			 const VECP rf, const VECP rc,
			 MATP Lxx, MATP Luu, MATP Lxu)
{
  ((Hqp_Docp_stub *)clientdata)->update_stage(k, x, u, f, f0, c,
					      fx, fu, f0x, f0u, cx, cu,
					      rf, rc, Lxx, Luu, Lxu);
}

//-------------------------------------------------------------------------
Hqp_Docp_stub::Hqp_Docp_stub()
{
  _k0 = 0;
  _kf = 0;

  Hqp_Docp_spec spec;
  spec.setup_horizon = ::setup_horizon;
  spec.setup_vars = ::setup_vars;
  spec.setup_struct = ::setup_struct;
  spec.init_simulation = ::init_simulation;
  spec.update_vals = ::update_vals;
  spec.update_stage = ::update_stage;

  _handle = Hqp_Docp_create(spec, this);
}

//-------------------------------------------------------------------------
Hqp_Docp_stub::~Hqp_Docp_stub()
{
  Hqp_Docp_destroy(_handle);
}

//-------------------------------------------------------------------------
void Hqp_Docp_stub::horizon(int k0, int kf)
{
  _k0 = k0;
  _kf = kf;
}

//-------------------------------------------------------------------------
void Hqp_Docp_stub::setup_horizon(int &k0, int &kf)
{
  k0 = _k0;
  kf = _kf;
}

//-------------------------------------------------------------------------
void Hqp_Docp_stub::alloc_vars(VECP v, VECP vmin, VECP vmax, int n)
{
  Hqp_Docp_alloc_vars(_handle, v, vmin, vmax, n);
}


//=========================================================================

/**
 * @file Hqp_DocpWrapper.h
 *   Implement an Hqp_Docp and connect it to Hqp_DocpSpec.
 *
 * rf, 10/31/00
 *
 */

/*
    Copyright (C) 1994--2009  Ruediger Franke

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

#ifndef Hqp_DocpWrapper_H
#define Hqp_DocpWrapper_H

#include <assert.h>

#include "Hqp_Docp.h"
#include "Hqp.h"

//---------------------------------------
class Hqp_DocpWrapper: public Hqp_Docp {
public:
  //----------------------------------------------------------
  Hqp_DocpWrapper(Hqp_DocpSpec &spec, void *clientdata)
  {
    // check for mandatory functions
    assert(spec.setup_horizon);
    assert(spec.setup_vars);
    assert(spec.update_vals);

    // initialize data
    _spec = spec;
    _clientdata = clientdata;

    // export reference to this object
    extern Hqp_SqpProgram *theSqpProgram;
    theSqpProgram = this;
  }

  //----------------------------------------------------------
  ~Hqp_DocpWrapper()
  {
    // clear reference to this object
    extern Hqp_SqpProgram *theSqpProgram;
    if (theSqpProgram == this)
      theSqpProgram = NULL;
  }

  // implement interface of Hqp_Docp by calling Hqp_DocpSpec
  //----------------------------------------------------------
  void setup_horizon(int &k0, int &kf)
  {
    (*_spec.setup_horizon)(_clientdata, k0, kf);
  }

  //----------------------------------------------------------
  void setup_vars(int k,
		  VECP x, VECP x_min, VECP x_max, IVECP x_int,
		  VECP u, VECP u_min, VECP u_max, IVECP u_int,
		  VECP c, VECP c_min, VECP c_max)
  {
    (*_spec.setup_vars)(_clientdata, k,
			x, x_min, x_max, x_int,
			u, u_min, u_max, u_int,
			c, c_min, c_max);
  }

  //----------------------------------------------------------
  void setup_struct(int k, const VECP x, const VECP u,
		    MATP fx, MATP fu, IVECP f_lin,
		    VECP f0x, VECP f0u, int &f0_lin,
		    MATP cx, MATP cu, IVECP c_lin,
		    MATP Lxx, MATP Luu, MATP Lxu)
  {
    if (_spec.setup_struct) {
      (*_spec.setup_struct)(_clientdata, k, x, u,
			    fx, fu, f_lin,
			    f0x, f0u, f0_lin,
			    cx, cu, c_lin,
			    Lxx, Luu, Lxu);
    }
    else {
      Hqp_Docp::setup_struct(k, x, u,
			     fx, fu, f_lin,
			     f0x, f0u, f0_lin,
			     cx, cu, c_lin,
			     Lxx, Luu, Lxu);
    }
  }

  //----------------------------------------------------------
  void init_simulation(int k, VECP x, VECP u)
  {
    if (_spec.init_simulation) {
      (*_spec.init_simulation)(_clientdata, k, x, u);
    }
    else {
      Hqp_Docp::init_simulation(k, x, u);
    }
  }

  //----------------------------------------------------------
  void update_vals(int k, const VECP x, const VECP u,
		   VECP f, double &f0, VECP c)
  {
    (*_spec.update_vals)(_clientdata, k, x, u, f, f0, c);
  }

  //----------------------------------------------------------
  void update_stage(int k, const VECP x, const VECP u,
		    VECP f, double &f0, VECP c,
		    MATP fx, MATP fu, VECP f0x, VECP f0u,
		    MATP cx, MATP cu,
		    const VECP rf, const VECP rc,
		    MATP Lxx, MATP Luu, MATP Lxu)
  {
    if (_spec.update_stage) {
      (*_spec.update_stage)(_clientdata, k, x, u, f, f0, c,
			    fx, fu, f0x, f0u, cx, cu,
			    rf, rc, Lxx, Luu, Lxu);
    }
    else {
      Hqp_Docp::update_stage(k, x, u, f, f0, c,
			     fx, fu, f0x, f0u, cx, cu,
			     rf, rc, Lxx, Luu, Lxu);
    }
  }

  //----------------------------------------------------------
  char *name() {return "DocpWrapper";}

protected:
  Hqp_DocpSpec		_spec;
  void			*_clientdata;
};  


#endif

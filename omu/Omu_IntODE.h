/*
 * Omu_IntODE.h --
 *   -- base class for integrating an ODE over a stage
 *   -- derived classes add a specific integration rule
 *
 * rf, 11/13/96
 */

/*
    Copyright (C) 1997--2001  Ruediger Franke

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

#ifndef Omu_IntODE_H
#define Omu_IntODE_H

#include "Omu_Integrator.h"

//--------------------------------------------------------------------------
class Omu_IntODE: public Omu_Integrator {

 public:

  Omu_IntODE();
  ~Omu_IntODE();

  // interface routine from Omu_Integrator
  void init_stage(int k,
		  const Omu_States &x, const Omu_Vector &u,
		  bool sa = false);
  void solve(int kk, Real tstart, Real tend,
	     const Omu_States &x, const Omu_Vector &u,
	     Omu_Program *sys, Omu_DepVec &Fc, Omu_SVec &xc);

  // routines provided by derived classes
  virtual void ode_solve(Real tstart, VECP y, VECP u, Real tend) = 0;

  // routines to be called by derived classes
  void syseq(Real t, const VECP y, const VECP u, VECP f);

 private:

  /**
   * alternative syseq implementation calling high-level
   * continuous and using forward for derivatives
   */
  void syseq_forward(Real t, const VECP y, const VECP u, VECP f);

  // backing store sys and current stage
  Omu_Program	*_sys;
  Omu_SVec 	*_xc_ptr;
  Omu_DepVec 	*_Fc_ptr;

  VECP		_y;
  VECP		_u;

  void		realloc();

  // vectors and matrices for low level _sys->continuous callback
  Omu_Vec	_uc;
  Omu_SVec	_xcp;
  MATP		_Yx;
  MATP		_Yu;

  // variables for ADOL-C
  VECP		_x;
  VECP		_v;
  MATP		_X;
  MATP		_Y;
};  

#endif

/*
 * Omu_IntRKsuite.h --
 *   -- integrate ODE over a stage using RKsuite
 *
 * rf, 11/13/96
 */

/*
    Copyright (C) 1997--1998  Ruediger Franke

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

#ifndef Omu_IntRKsuite_H
#define Omu_IntRKsuite_H

#include "Omu_IntODE.h"

//--------------------------------------------------------------------------
class Omu_IntRKsuite: public Omu_IntODE {

 public:

  Omu_IntRKsuite();
  ~Omu_IntRKsuite();

  char *name() {return "RKsuite";}

  // interface routine
  void ode_solve(Real tstart, VECP y, const VECP u, Real tend);

  // callback for calculation of system equations
  void F(Real *T, Real *Y, Real *YP);

 private:

  int  _neq;
  int  _npar;
/*    Real _rtol; */
/*    Real _atol; */
  Real _hnext;
  Real _tlast;
  int  _method;
  VECP _thres;
  VECP _work;
  VECP _u;
  VECP _yp;
  VEC _y_head;
  VEC _yp_head;
};  

#endif


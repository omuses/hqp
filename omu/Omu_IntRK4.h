/*
 * Omu_IntRK4.h --
 *   -- integrate ODE over a stage using a simple RK4 rule
 *
 * rf, 10/2/96
 *
 * rf, 11/13/96
 *   -- use base class Omu_IntODE
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

#ifndef Omu_IntRK4_H
#define Omu_IntRK4_H

#include "Omu_IntODE.h"

/**
 * Explicit, fixed step size integrator implementing a Runge Kutta
 * algorithm of 4th order.
 * If _stepsize is not specified, then the integrator takes one step
 * per sample period.
 */
class Omu_IntRK4: public Omu_IntODE {

 public:

  Omu_IntRK4(); 	///< constructor
  ~Omu_IntRK4(); 	///< destructor

  /**
   * @name Implementation of predefined methods.
   * @see Omu_IntODE
   */

  //@{

  const char *name() {return "RK4";}

  void ode_solve(double tstart, VECP y, const VECP u, double tend);

  //@}

 private:

  VECP		_y1;
  VECP		_yh;
  VECP		_fh;
};  

#endif


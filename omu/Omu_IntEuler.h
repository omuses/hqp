/*
 * Omu_IntEuler.h --
 *   -- integrate ODE over a stage using fixed step size Euler rule
 *
 * rf, 02/07/01
 *
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

#ifndef Omu_IntEuler_H
#define Omu_IntEuler_H

#include "Omu_IntODE.h"

/**
 * Explicit, fixed step size integrator implementing the Euler algorithm.
 * If _stepsize is not specified, then the integrator takes one step
 * per sample period.
 */
class Omu_IntEuler: public Omu_IntODE {

 public:

  Omu_IntEuler(); 	///< constructor
  ~Omu_IntEuler(); 	///< destructor

  /**
   * @name Implementation of predefined methods.
   * @see Omu_IntODE
   */

  //@{

  char *name() {return "Euler";}

  void ode_solve(double tstart, VECP y, const VECP u, double tend);

  //@}

 private:

  VECP		_f;
};  

#endif


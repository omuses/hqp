/*
 * Omu_IntRKsuite.h --
 *   -- integrate ODE over a stage using RKsuite
 *
 * rf, 11/13/96
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

#ifndef Omu_IntRKsuite_H
#define Omu_IntRKsuite_H

#include "Omu_IntODE.h"

/**
 * Solve ordinary differential equations using RKsuite.
 * @see http://www.netlib.org
 */
class Omu_IntRKsuite: public Omu_IntODE {

 public:

  Omu_IntRKsuite(); 	///< constructor
  ~Omu_IntRKsuite(); 	///< destructor

  /**
   * @name Implementation of predefined methods.
   * @see Omu_IntODE
   */

  //@{

  char *name() {return "RKsuite";}

  void ode_solve(double tstart, VECP y, const VECP u, double tend);

  //@}

  /** Callback routine for RKsuite for the calculation of system equations. */
  void F(double *T, double *Y, double *YP);

 protected:

  /**
   * Specifies which Runge-Kutta method is to be used (default: 2).
   * Possible choices are:
   *  - 1: use the (2,3) Runge-Kutta pair
   *  - 2: use the (4,5) Runge-Kutta pair
   *  - 3: use the (7,8) Runge-Kutta pair
   */
  int  _method;

 private:

  int  _neq;
  int  _npar;
  double _hnext;
  double _tlast;
  VECP _thres;
  VECP _work;
  VECP _u;
  VECP _yp;
  VEC _y_head;
  VEC _yp_head;
};  

#endif

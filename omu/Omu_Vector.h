/*
 * Omu_Vector.h --
 *   -- variable vector for Omuses problem setup
 *
 * rf, 2/3/97
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

#ifndef Omu_Vector_H
#define Omu_Vector_H

#include <math.h>
#include <Meschach.h>

#define Inf HUGE_VAL

//--------------------------------------------------------------------------
class Omu_Vector: public VECP {
  			// base class holds variable values
 public:
  VECP min;		// minimal permitted values (default: -Inf)
  VECP max;		// maximal permitted values (default: +Inf)
  VECP initial;		// initial values (default: 0.0)

  Omu_Vector();
  virtual ~Omu_Vector();

  virtual void alloc(int n, int n_expand = -1);
  	// -- allocate vectors of size n_expand 
        //    for min, max, initial, and this
	// -- (default n_expand: n)
	// -- derived classes may overload alloc() to restrict 
        //    the capability of allocations for specific vectors

 private:
	// dismiss copy constructor and operator=
  Omu_Vector(const Omu_Vector &);
  Omu_Vector &operator=(const Omu_Vector &);
};

#endif

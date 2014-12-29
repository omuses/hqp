/**
 * @file Omu_Vec.h
 *    Vector with automatic construction and destruction
 *
 * rf, 7/31/01
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

/** Avoid multiple inclusion */
#if !defined(Omu_Vec_H)
#define Omu_Vec_H

#include <Meschach.h>
#include "Omu.h"

/** Vector class with automatic construction / destruction */
class OMU_API Omu_Vec: public Mesch::VECP {
public:
  Omu_Vec() {_v = v_resize(v_get(1), 0);} ///< allocate empty vector
  virtual ~Omu_Vec() {v_free(_v);} 	  ///< destroy vector

  void resize(int dim) {v_resize(_v, dim);}///< resize vector

 protected:
  /// protect copy constructor and operator= as they should not be used
  //@{
  Omu_Vec(const Omu_Vec &v): VECP() {_v = v_copy(v, _v);}
  Omu_Vec &operator=(const Omu_Vec &v) {
    _v = v_copy(v, _v);
    return *this;
  }
  //@}
};


#endif

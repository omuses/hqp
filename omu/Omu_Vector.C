/*
 * Omu_Vector.C
 *   -- class definition
 *
 * rf, 2/3/97
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

#include "Omu_Vector.h"

//--------------------------------------------------------------------------
Omu_Vector::Omu_Vector()
{
  _v = v_resize(v_get(1), 0);
  min = v_resize(v_get(1), 0);
  max = v_resize(v_get(1), 0);
  initial = v_resize(v_get(1), 0);
}

//--------------------------------------------------------------------------
Omu_Vector::Omu_Vector(const Omu_Vector &cv)
  : VECP()
{
  _v = v_copy(cv._v, VNULL);
  min = v_copy(cv.min, VNULL);
  max = v_copy(cv.max, VNULL);
  initial = v_copy(cv.initial, VNULL);
}

//--------------------------------------------------------------------------
Omu_Vector &Omu_Vector::operator=(const Omu_Vector &cv)
{
  v_copy(cv._v, _v);
  v_copy(cv.min, min);
  v_copy(cv.max, max);
  v_copy(cv.initial, initial);

  return *this;
}

//--------------------------------------------------------------------------
Omu_Vector::~Omu_Vector()
{
  v_free(initial);
  v_free(max);
  v_free(min);
  v_free(_v);
}

//--------------------------------------------------------------------------
void Omu_Vector::alloc(int n, int n_expand)
{
  if (n_expand < n)
    n_expand = n;

  v_resize(_v, n_expand);
  v_resize(initial, n_expand);
  v_resize(min, n);
  v_resize(max, n);

  v_set(_v, 0.0);
  v_set(initial, 0.0);
  v_set(min, -Inf);
  v_set(max, Inf);
}


//==========================================================================

/*
 * Omu_Variables.C
 *   -- class definitions
 *
 * rf, 10/10/01
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

#include "Omu_Variables.h"
#include <Hqp.h>  // include this file to get Inf

//--------------------------------------------------------------------------
Omu_VariableVec::Omu_VariableVec()
{
  min = v_resize(v_get(1), 0);
  max = v_resize(v_get(1), 0);
  initial = v_resize(v_get(1), 0);
  integer = iv_resize(iv_get(1), 0);
  _n = 0;
}

//--------------------------------------------------------------------------
Omu_VariableVec::Omu_VariableVec(const Omu_VariableVec &cv)
  : Omu_Vec(cv)
{
  min = v_copy(cv.min, VNULL);
  max = v_copy(cv.max, VNULL);
  initial = v_copy(cv.initial, VNULL);
  integer = iv_copy(cv.integer, IVNULL);
}

//--------------------------------------------------------------------------
Omu_VariableVec &Omu_VariableVec::operator=(const Omu_VariableVec &cv)
{
  v_copy(cv._v, _v);
  v_copy(cv.min, min);
  v_copy(cv.max, max);
  v_copy(cv.initial, initial);
  iv_copy(cv.integer, integer);

  return *this;
}

//--------------------------------------------------------------------------
Omu_VariableVec::~Omu_VariableVec()
{
  iv_free(integer);
  v_free(initial);
  v_free(max);
  v_free(min);
}

//--------------------------------------------------------------------------
void Omu_VariableVec::alloc(int n, int n_expand)
{
  if (n_expand < n)
    n_expand = n;
  _n = n;

  v_resize(_v, n_expand);
  v_resize(min, n_expand);
  v_resize(max, n_expand);
  v_resize(initial, n_expand);
  iv_resize(integer, n_expand);

  v_set(_v, 0.0);
  v_set(min, -Inf);
  v_set(max, Inf);
  v_set(initial, 0.0);
  iv_set(integer, 0);
}

//==========================================================================

//--------------------------------------------------------------------------
Omu_OptVarVec::Omu_OptVarVec()
{
  weight1 = v_resize(v_get(1), 0);
  weight2 = v_resize(v_get(1), 0);
  ref = v_resize(v_get(1), 0);
  active = iv_resize(iv_get(1), 0);
}

//--------------------------------------------------------------------------
Omu_OptVarVec::~Omu_OptVarVec()
{
  v_free(ref);
  v_free(weight2);
  v_free(weight1);
  iv_free(active);
}

//--------------------------------------------------------------------------
void Omu_OptVarVec::resize(int n)
{
  Omu_VariableVec::alloc(n);
  v_resize(weight1, n);
  v_resize(weight2, n);
  v_resize(ref, n);
  iv_resize(active, n);

  v_set(weight1, 0.0);
  v_set(weight2, 0.0);
  v_set(ref, 0.0);
  iv_zero(active);
}

//==========================================================================

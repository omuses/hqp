/*
 * Omu_Vars.C
 *   -- class definitions
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

#include "Omu_Vars.h"

//--------------------------------------------------------------------------
Omu_Vars::Omu_Vars()
{
  c_alloc = false;
  c_expand = false;
}

//--------------------------------------------------------------------------
void Omu_Vars::alloc(int n, int n_expand)
{
  n_expand = max(n, n_expand);

  if (!c_alloc)
    error(E_OVERWRITE, "Omu_Vars::alloc");

  if (!c_expand && n_expand > n)
    error(E_RANGE, "Omu_Vars::alloc_expanded");

  v_resize(_v, n_expand);
  v_resize(initial, n_expand);
  v_resize(min, n_expand);
  v_resize(max, n_expand);

  v_set(_v, 0.0);
  v_set(initial, 0.0);
  v_set(min, -Inf);
  v_set(max, Inf);
}

const int Omu_States::Discrete = 1;
const int Omu_States::Algebraic = 2;

//--------------------------------------------------------------------------
Omu_States::Omu_States()
{
  nd = -1;	// must be determined
  na = 0;
  nv = 0;
  D = v_resize(v_get(1), 0);
  D_is_const = false;
  flags = iv_resize(iv_get(1), 0);
  sbw_u = -1;
  sbw_l = -1;
}

//--------------------------------------------------------------------------
Omu_States::~Omu_States()
{
  iv_free(flags);
  v_free(D);
}

//--------------------------------------------------------------------------
void Omu_States::alloc(int n, int n_expand)
{
  n_expand = max(n, n_expand);
  Omu_Vars::alloc(n, n_expand);
  nv = n_expand - n;
  v_resize(D, n_expand);
  v_zero(D);
  iv_resize(flags, n_expand);
  iv_zero(flags);
}


//==========================================================================

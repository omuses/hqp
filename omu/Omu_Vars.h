/*
 * Omu_Vars.h --
 *   -- Omu_Vector with capabilities
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

#ifndef Omu_Vars_H
#define Omu_Vars_H

#include "Omu_Vector.h"

//--------------------------------------------------------------------------
class Omu_Vars: public Omu_Vector {
 public:
  bool c_alloc;
  bool c_expand;

  Omu_Vars();

  void alloc(int n, int n_expand = -1);
};

//--------------------------------------------------------------------------
class Omu_States: public Omu_Vars {
 public:
  int 	nd;		// number of discrete-time states
  int	na;		// number of algebraic states
  int	nv;		// number of expansion variables
  VECP	D;		// diagonal of dF/dxp
  bool	D_is_const;	// F can be treated as explicit ODE
  int	sbw_u;		// upper semi-bandwidth of dF/dx + dF/dxp
  int	sbw_l;		// lower semi-bandwidth of dF/dx + dF/dxp

  /**
   * Define flagbits to characterize individual states.
   * This makes nd and nv obsolete. They should not be accessed
   * anymore as they will be removed in a later version.
   */
  IVECP flags;

  /**
   * Indicate a discrete-time state, which is not defined with F.
   */
  static const int Discrete;

  /**
   * Indicate an algebraic state, which is defined with F,
   * but whose time derivative does not appear.
   */
  static const int Algebraic;

  Omu_States();
  ~Omu_States();

  void alloc(int n, int n_expand = -1);
};

#endif

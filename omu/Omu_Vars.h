/*
 * Omu_Vars.h --
 *   -- extensions for managing independent variables
 *
 * rf, 2/3/97
 */

/*
    Copyright (C) 1997--2001  Ruediger Franke

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

#include "Omu_Variables.h"

/** Internally used vector of optimization variables. */
class Omu_VarVec: public Omu_VariableVec {
 public:
  bool c_alloc;
  bool c_expand;

  Omu_VarVec();

  void alloc(int n, int n_expand = -1);
};

/** Internally used vector of dynamic variables with additional structural
    information. */
class Omu_DynVarVec: public Omu_VarVec {
 public:
  int 	nd;		// number of discrete-time state variables
  int	na;		// number of algebraic state variables
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

  Omu_DynVarVec();
  ~Omu_DynVarVec();

  void alloc(int n, int n_expand = -1);
};

/** Depreciated name for Omu_DynVarVec. */
typedef Omu_DynVarVec Omu_States;

/** Internally used vector of state variables. */
class Omu_SVec: public Omu_StateVec {
public:
  void alloc(int dim, int nx, int nu);
  void realloc(int dim, int nx, int nu) {alloc(dim, nx, nu);}
};

#endif

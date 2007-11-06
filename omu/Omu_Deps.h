/**
 * @file Omu_Deps.h
 *    Extensions for managing dependent variables internally
 *
 * rf, 7/31/01
 */

/*
    Copyright (C) 1997--2007  Ruediger Franke

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

#ifndef Omu_Deps_H
#define Omu_Deps_H

#include "Omu_Dependents.h"


/** Internally used single dependent variable. */
class Omu_Dep: public Omu_Dependent {
public:
  Omu_Dep();

  /** Assign a new value */
  const Omu_Dep &operator = (const Real &value) {
    _value = value;
    return *this;
  }

  /** Allocate gradients for given numbers of independent variables */
  void size(int nx, int nu, int nxf);

  bool c_setup;  ///< indicate setup process (call to set_linear* allowed)

  /** @name Implementation of Omu_Dependent interface */
  //@{
  void set_linear(int wrt = Omu_Dependent::WRT_ALL, bool value = true) {
    if (!c_setup)
      m_error(E_OVERWRITE, "Omu_Dep::set_linear");
    if (value)
      _linear_flags |= wrt;
    else
      _linear_flags &= ~wrt;
  }
  bool is_linear(int wrt = Omu_Dependent::WRT_ALL) const {
    return (_linear_flags & wrt) == wrt;
  }
  //@}

  /** Analyze dependencies after setup process */
  void analyze_struct();

protected:
  int _linear_flags;

  // dismiss copy constructor and assignment operator
  Omu_Dep(const Omu_Dep &cv): Omu_Dependent(cv) {}
  Omu_Dep &operator=(const Omu_Dep &) {return *this;}
};


/** Internally used vector of dependent variables. */
class Omu_DepVec: public Omu_DependentVec {
public:
  Omu_DepVec();
  ~Omu_DepVec();

  /** Allocate and initialize memory */
  void size(int dim, int nx, int nu, int ndx, int nxf, int nq);

  /** Resize dimension without reinitializing memory.
      Argument dim must not be larger than allocated dim. */
  void adapt_size(int dim);

  bool c_setup; ///< indicate setup process (call to set_linear* allowed)

  /** @name Implementation of Omu_DependentVec interface */
  //@{
  void set_linear(int wrt = Omu_Dependent::WRT_ALL, bool value = true) {
    if (!c_setup)
      m_error(E_OVERWRITE, "Omu_DepVec::set_linear");

    // mark elements linear
    for (int i = 0; i < (int)_linear_flags->dim; i++)
      if (value)
	_linear_flags[i] |= wrt;
      else
	_linear_flags[i] &= ~wrt;

    // mark variables linear
    int val = value? 1: 0;
    int offs = 0;
    if ((wrt & Omu_Dependent::WRT_x))
      for (int j = 0; j < Jx->n; j++)
        _linear_vars[offs + j] = val;
    offs += Jx->n;
    if ((wrt & Omu_Dependent::WRT_u))
      for (int j = 0; j < Ju->n; j++)
        _linear_vars[offs + j] = val;
    offs += Ju->n;
    if ((wrt & Omu_Dependent::WRT_dx))
      for (int j = 0; j < Jdx->n; j++)
        _linear_vars[offs + j] = val;
    offs += Jdx->n;
    if ((wrt & Omu_Dependent::WRT_xf))
      for (int j = 0; j < Jxf->n; j++)
        _linear_vars[offs + j] = val;
    offs += Jxf->n;
    if ((wrt & Omu_Dependent::WRT_q))
      for (int j = 0; j < Jq->n; j++)
        _linear_vars[offs + j] = val;
  }
  bool is_linear(int wrt = Omu_Dependent::WRT_ALL) const {
    bool result = true;
    for (int i = 0; i < (int)_linear_flags->dim; i++)
      result &= (_linear_flags[i] & wrt) == wrt;
    return result;
  }

  void set_linear_element(int i, int wrt = Omu_Dependent::WRT_ALL,
			  bool value = true) {
    if (!c_setup)
      m_error(E_OVERWRITE, "Omu_DepVec::set_linear_element");
    if (value)
      _linear_flags[i] |= wrt;
    else
      _linear_flags[i] &= ~wrt;
  }
  bool is_linear_element(int i, int wrt = Omu_Dependent::WRT_ALL) const {
    return (_linear_flags[i] & wrt) == wrt;
  }

  void set_linear_variable(int vec, int j, bool value = true) {
    if (!c_setup)
      m_error(E_OVERWRITE, "Omu_DepVec::set_linear_variable");
    int offs = vec > Omu_Dependent::WRT_x? Jx->n: 0
      + vec > Omu_Dependent::WRT_u? Ju->n: 0
      + vec > Omu_Dependent::WRT_dx? Jdx->n: 0
      + vec > Omu_Dependent::WRT_xf? Jxf->n: 0
      + vec > Omu_Dependent::WRT_q? Jq->n: 0;
    if (value)
      _linear_vars[offs + j] = 1;
    else
      _linear_vars[offs + j] = 0;
  }

  bool is_linear_variable(int vec, int j) const {
    int offs = vec > Omu_Dependent::WRT_x? Jx->n: 0
      + vec > Omu_Dependent::WRT_u? Ju->n: 0
      + vec > Omu_Dependent::WRT_dx? Jdx->n: 0
      + vec > Omu_Dependent::WRT_xf? Jxf->n: 0
      + vec > Omu_Dependent::WRT_q? Jq->n: 0;
    return _linear_vars[offs + j] != 0;
  }
  //@}

  /** Analyze dependencies after setup process */
  void analyze_struct();

protected:
  IVECP _linear_flags;
  IVECP _linear_vars;

  // dismiss copy constructor and assignment operator
  Omu_DepVec(const Omu_DepVec &cv): Omu_DependentVec(cv) {}
  Omu_DepVec &operator=(const Omu_DepVec &) {return *this;}
};


#endif

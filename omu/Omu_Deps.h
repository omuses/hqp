/**
 * @file Omu_Deps.h
 *    Extensions for managing dependent variables internally
 *
 * rf, 7/31/01
 */

/*
    Copyright (C) 1997--2003  Ruediger Franke

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
  bool is_linear(int wrt = Omu_Dependent::WRT_ALL) {
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
  void size(int dim, int nx, int nu, int nxp, int nxf);

  /** Resize dimension without reinitializing memory.
      Argument dim must not be larger than allocated dim. */
  void adapt_size(int dim);

  bool c_setup; ///< indicate setup process (call to set_linear* allowed)

  /** @name Implementation of Omu_DependentVec interface */
  //@{
  void set_linear(int wrt = Omu_Dependent::WRT_ALL, bool value = true) {
    if (!c_setup)
      m_error(E_OVERWRITE, "Omu_DepVec::set_linear");
    for (int i = 0; i < (int)_linear_flags->dim; i++)
      if (value)
	_linear_flags[i] |= wrt;
      else
	_linear_flags[i] &= ~wrt;
  }
  bool is_linear(int wrt = Omu_Dependent::WRT_ALL) {
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
  bool is_linear_element(int i, int wrt = Omu_Dependent::WRT_ALL) {
    return (_linear_flags[i] & wrt) == wrt;
  }
  //@}

  /** Analyze dependencies after setup process */
  void analyze_struct();

protected:
  IVECP _linear_flags;

  // dismiss copy constructor and assignment operator
  Omu_DepVec(const Omu_DepVec &cv): Omu_DependentVec(cv) {}
  Omu_DepVec &operator=(const Omu_DepVec &) {return *this;}
};


#endif

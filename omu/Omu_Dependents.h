/**
 * @file Omu_Dependents.h 
 *   Dependent variables (single and vector, including Jacobians).
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

#ifndef Omu_Dependents_H
#define Omu_Dependents_H

#include <Meschach.h>
#include "Omu_Vec.h"

/** Vector extended with structural information for a gradient. */
class Omu_Gradient: public Mesch::VECP {
  friend class Omu_Dep;
public:
  /** Create an empty matrix. */
  Omu_Gradient();

  /** Free memory. */
  ~Omu_Gradient();

  // access methods for properties
  bool is_zero() {return _is_zero;}
  bool is_constant() {return _is_constant;}

protected:
  /** Allocate and initialize Gradient. */
  void alloc(int dim);

  /** Obtain properties for current matrix. */
  void analyze_struct(bool is_constant);

  bool _is_zero;
  bool _is_constant;
};


/** Single dependent variable. */
class Omu_Dependent {
public:
  // flags for indicating dependencies
  static const int WRT_x;
  static const int WRT_u;
  static const int WRT_xp;
  static const int WRT_xf;
  static const int WRT_ALL;

  Omu_Gradient gx;	// gradient w.r.t. x (i.e. initial states of period)
  Omu_Gradient gu;	// gradient w.r.t. u (i.e. control parameters of stage)
  Omu_Gradient gxf;	// gradient w.r.t. xf (i.e. final states of period)

  Omu_Dependent();
  virtual ~Omu_Dependent();

  const Omu_Dependent &operator = (const Real &value) {
    _value = value;
    return *this;
  }

  operator Real&() {
    return _value;
  }

  virtual void set_linear(int wrt = Omu_Dependent::WRT_ALL,
			  bool value = true) = 0;
  virtual bool is_linear(int wrt = Omu_Dependent::WRT_ALL) = 0;

  void set_required_g(bool value = true) {
    _required_g = value;
  }
  bool is_required_g() const {
    return _required_g;
  }

protected:
  Real _value;
  bool _required_g;
  // dismiss copy constructor and assignment operator
  Omu_Dependent(const Omu_Dependent &) {}
  Omu_Dependent &operator=(const Omu_Dependent &) {return *this;}
};


/** Matrix extended with structural information for a Jacobian. */
class Omu_Jacobian: public Mesch::MATP {
  friend class Omu_DepVec;
public:
  /** Create an empty matrix. */
  Omu_Jacobian();

  /** Free memory. */
  ~Omu_Jacobian();

  // access methods for properties
  bool is_zero() {return _is_zero;}
  bool is_ident() {return _is_ident;}
  bool is_constant() {return _is_constant;}
  int sbw() {return max(_sbw_lower, _sbw_upper);}
  int sbw_lower() {return _sbw_lower;}
  int sbw_upper() {return _sbw_upper;}

protected:
  /** Allocate and initialize Jacobian. */
  void alloc(int nrows, int ncols);

  /** Resize dimension without reinitializing memory.
      Argument nrows must not be larger than allocated nrows.*/
  void adapt_size(int nrows);

  /** Obtain properties for current matrix. */
  void analyze_struct(bool is_constant);

  bool _is_zero;
  bool _is_ident;
  bool _is_constant;
  int _sbw_lower;
  int _sbw_upper;
};


/** Vector of dependent variables. */
class Omu_DependentVec: public Omu_Vec {
public:
  Omu_Jacobian Jx;	// Jacobian w.r.t. x (initial states of sample period)
  Omu_Jacobian Ju;	// Jacobian w.r.t. u (control parameters of stage)
  Omu_Jacobian Jxp;	// Jacobian w.r.t. xp (time derivative of x)
  Omu_Jacobian Jxf;	// Jacobian w.r.t. xf (final states of sample period)

  Omu_DependentVec();

  virtual void set_linear(int wrt = Omu_Dependent::WRT_ALL,
			  bool value = true) = 0;
  virtual bool is_linear(int wrt = Omu_Dependent::WRT_ALL) = 0;

  virtual void set_linear_element(int i, int wrt = Omu_Dependent::WRT_ALL,
				  bool value = true) = 0;
  virtual bool is_linear_element(int i, int wrt = Omu_Dependent::WRT_ALL) = 0;

  void set_required_J(bool value = true) {
    _required_J = value;
  }
  bool is_required_J() const {
    return _required_J;
  }

protected:
  bool _required_J;

  // dismiss copy constructor and assignment operator
  Omu_DependentVec(const Omu_DependentVec &): Omu_Vec() {}
  Omu_DependentVec &operator=(const Omu_DependentVec &) {return *this;}
};


#endif

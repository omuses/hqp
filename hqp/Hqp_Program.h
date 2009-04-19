/**
 * @file Hqp_Program.h 
 *   Optimization program for sparse quadratic programming
 *
 * rf, 5/15/94
 */

/*
    Copyright (C) 1994--2009  Ruediger Franke

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

#ifndef Hqp_Program_H
#define Hqp_Program_H

#include "Hqp_impl.h"

/** Optimization program for sparse quadratic programming */
class Hqp_Program {

 public:

  /** @name flags for variables */
  //@{
  static const int IS_LOCAL; 	///< local variable of nonlinear group
  static const int IS_SLACK; 	///< slack variable
  //@}

  VECP		x;	///< optimization variables
  IVECP 	x_flags;///< classification of variables
  IVECP 	x_int; 	///< mark integer variables with values > 0

  SPMATP	Q;	///< criterion: 1/2 x'Qx + c'x -> min
  VECP		c; 	///< linear term of criterion

  SPMATP	A;	///< equality constraints: Ax + b = 0
  VECP		b; 	///< constant term of equality constraints

  SPMATP	C;	///< inequality constraints: Cx + d >= 0
  VECP		d; 	///< constant term of inequality constraints

  Hqp_Program(); 	///< constructor
  ~Hqp_Program(); 	///< destructor

  /** change size of optimization program */
  void resize(int n, int me, int m,
	      int el_n = 0, int el_me = 0, int el_m = 0);

  /** dump optimization program to FILE */
  void foutput(FILE *fp);
};  

#endif



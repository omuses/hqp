/*
 * Hqp_Program.h -- program for sparse quadratic programming
 *
 * rf, 5/15/94
 */

/*
    Copyright (C) 1994--2000  Ruediger Franke

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


class Hqp_Program {

 public:

  VECP		x;	// variables to optimize on

  SPMATP	Q;	// criterion: 1/2 x'Qx + c'x -> min
  VECP		c;

  SPMATP	A;	// equality constraints: Ax + b = 0
  VECP		b;

  SPMATP	C;	// inequality constraints: Cx + d >= 0
  VECP		d;

  Hqp_Program();
  ~Hqp_Program();

  void resize(int n, int me, int m,
	      int el_n = 0, int el_me = 0, int el_m = 0);
  void foutput(FILE *fp);
};  

#endif



/*
 * Hqp_DocpAdol.h --
 *   -- discrete time optimal control problems with automatic differentiation
 *      using operator overloading
 *
 *      Refences:
 *        E. Arnold: Hqp_Docp_Adol, personal communication
 * 
 *        A. Griewank, D. Juedes, and J. Utke:
 *           ADOL-C: A Package for the Automatic Differentiation of
 *           Algorithms Written in C/C++. 
 * rf, 12/19/95
 *
 * rf, 9/15/96
 *   -- replace update_ders() with update_stage()
 */

/*
    Copyright (C) 1994--1998  Ruediger Franke

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

#ifndef Hqp_DocpAdol_H
#define Hqp_DocpAdol_H

#include "Hqp_Docp.h"

#include <adolc/adolc.h>

//--------------------------------------------------------------------------
class Hqp_DocpAdol: public Hqp_Docp {

 public:

  Hqp_DocpAdol();
  ~Hqp_DocpAdol();

  void	horizon(int k0, int kmax);

  // interface for overloding by a derived class
  virtual void setup_constr(int k,
			    VECP xmin, VECP xmax,
			    VECP umin, VECP umax,
			    VECP cmin, VECP cmax) = 0;
  virtual void init_solution(int k,
			     VECP x, VECP u) = 0;
  virtual void update_avals(int k, 
			    const adoublev &x, const adoublev &u,
			    adoublev &f, adouble &f0, adoublev &c) = 0;

 private:

  int _k0;
  int _kmax;
  int _static_struct;
  int _hela;
  int _ad;

  // arrays passed to Adol-C and their dimensions
  int _am;	// number of dependent variables
  int _an;	// number of independent variables
  double **_X;
  double **_Y;
  double **_U;
  double **_Z;
  double ***_Z3;
  short  **_nz;
  void adalloc(int m, int n);
  void adfree();
  void adresize(int m, int n);

  // methods of Hqp_Docp, that are handled internally 
  // with the help of update_avals() and using Adol-C
  void setup_struct(int k,
		    MATP fx, MATP fu, IVECP f_lin,
		    MATP cx, MATP cu, IVECP c_lin,
		    MATP Lxx, MATP Luu, MATP Lxu);
  void update_vals(int k, const VECP x, const VECP u,
		   VECP f, Real &f0, VECP c);
  void update_stage(int k, const VECP x, const VECP u,
		    VECP f, Real &f0, VECP c,
		    MATP fx, MATP fu, VECP f0x, VECP f0u,
		    MATP cx, MATP cu,
		    const VECP rf, const VECP rc,
		    MATP Lxx, MATP Luu, MATP Lxu);
};  

#endif


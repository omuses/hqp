/*
 * Omu_IntOdeTs.h --
 *   -- class for integrating an Ode over a stage
 *   -- using Taylor series
 *
 * hl, 28/01/98
 */

/*
    Copyright (C) 1997--1999  Hartmut Linke

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

#ifndef Omu_IntOdeTs_H
#define Omu_IntOdeTs_H

#include "Omu_Integrator.h"

//--------------------------------------------------------------------------
class Omu_IntOdeTs: public Omu_Integrator {

 public:

  Omu_IntOdeTs(int );
  Omu_IntOdeTs();
  ~Omu_IntOdeTs();

  char *name() {return "OdeTs";}

  // interface routine from Omu_Integrator
  void solve(int kk, Real tstart, Real tend,
	     const Omu_States &x, const Omu_Vector &u,
	     Omu_Program *sys, VECP xt,
	     MATP Sx, MATP Su);

  private:

  bool          _adtaylor;
  bool          _multiple_record;

  int           _output;

  int		_n;	// number of continuous states
  int		_m;	// number of controls

  int           _max_deg;
  int           _max_deg0;

  double        _rho;
  double	_tau;

/*    double    _atol; */
/*    double  	_rtol; */

  VECP		_x;
  VECP		_u;

  VECP          _a;

  MATP          _Sxd;
  MATP          _Sxd0;
  MATP          _Sxx;
  MATP          _Sxx0;
  MATP          _Sxu;
  MATP          _Sxu0;
  MATP          _Sh;

  short         **_nz;         // sparsity pattern

  double        *_aa1;
  double        *_ab1;
  double        **_aa;
  double        **_ab;
  double        **_ai;
  double        ***_aA;
  double        ***_aB;

  void		realloc(int , int , int , int , int );
  void          free_adtaylor();

};

#endif





/*
 * Omu_IntIMP.h --
 *   -- integrate an ODE over a stage using the implicit midpoint rule
 *
 * E. Arnold   1999-04-12
 *             2000-05-12 _rtol, _atol -> Omu_Integrator 
 *             2000-05-30 step size control
 *             2001-08-16 prg_int_nsteps --> _stepsize
 *
 */

/*
    Copyright (C) 1999--2000  Eckhard Arnold

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

#ifndef Omu_IntIMP_H
#define Omu_IntIMP_H

#include "Omu_IntODE.h"

//--------------------------------------------------------------------------
class Omu_IntIMP: public Omu_IntODE {

 public:

  Omu_IntIMP();
  ~Omu_IntIMP();

  char *name() {return "IMP";}

  // interface routines
  void init_stage(int k,
		  const Omu_States &x, const Omu_Vector &u,
		  const Omu_DepVec &F, bool sa = false);

  void ode_solve(double tstart, VECP y, const VECP u, double tend);

 private:

  void resize();
  void jac(double t, VECP y);
  void step(double tstep, double dt, VECP y);

  VECP		_y;
  VECP		_u;
  VECP		_k1;
  VECP		_y0;
  VECP		_res;
  VECP		_fh;
  VECP		_yjac;
  VECP		_yjacp;
  MATP		_yy;
  MATP		_yyn;
  PERM          *_ppivot;
  VECP		_y1;
  VECP		_y2;

  int           _npar;
  int           _max_modnewtonsteps;
  int           _modnewtonsteps;
  int           _maxiters;
  double        _hinit;
  double        _dt;
  int           _IMP;

};  

#endif


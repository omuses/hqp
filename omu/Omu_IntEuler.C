/*
 * Omu_IntEuler.C --
 *   -- class implementation
 *
 * rf, 10/2/96
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

#include <assert.h>

#include <If_Class.h>

#include "Omu_IntEuler.h"

IF_CLASS_DEFINE("Euler", Omu_IntEuler, Omu_Integrator);

//--------------------------------------------------------------------------
Omu_IntEuler::Omu_IntEuler()
{
  _f = v_get(1);
}

//--------------------------------------------------------------------------
Omu_IntEuler::~Omu_IntEuler()
{
  v_free(_f);
}

//--------------------------------------------------------------------------
void Omu_IntEuler::ode_solve(Real tstart, VECP y, const VECP u, Real tend)
{
  Real dt;
  int i, nsteps;
  int neq = y->dim;

  v_resize(_f, neq);

  // obtain number of steps
  if (_stepsize > 0.0)
    nsteps = (int)ceil((tend - tstart) / _stepsize);
  else
    nsteps = 1;
  dt = (tend - tstart) / nsteps;

  // apply Euler's integration rule
  for (i = 0; i < nsteps; ++i) {
    syseq(tstart + i*dt, y, u, _f);
    v_mltadd(y, _f, dt, y);
  }
}


//========================================================================

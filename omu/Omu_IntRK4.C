/*
 * Omu_IntRK4.C --
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

#include <adutils.h>

#include <If_Int.h>
#include <If_Float.h>
#include <If_Class.h>

#include "Omu_IntRK4.h"

IF_CLASS_DEFINE("RK4", Omu_IntRK4, Omu_Integrator);

//--------------------------------------------------------------------------
Omu_IntRK4::Omu_IntRK4()
{
  _y1 = v_get(1);
  _yh = v_get(1);
  _fh = v_get(1);
}

//--------------------------------------------------------------------------
Omu_IntRK4::~Omu_IntRK4()
{
  v_free(_fh);
  v_free(_yh);
  v_free(_y1);
}

//--------------------------------------------------------------------------
void Omu_IntRK4::ode_solve(Real tstart, VECP y, const VECP u, Real tend)
{
  Real t, dt;
  int i, nsteps;
  int neq = y->dim;

  v_resize(_y1, neq);
  v_resize(_yh, neq);
  v_resize(_fh, neq);

  if (_stepsize > 0.0)
    nsteps = (int)ceil((tend - tstart) / _stepsize);
  else
    nsteps = 1;
  dt = (tend - tstart) / nsteps;

  for (i = 0; i < nsteps; ++i) {
    t = tstart + i*dt;

    syseq(t, y, u, _fh);
    v_mltadd(y, _fh, dt/6.0, _y1);

    v_mltadd(y, _fh, dt/2.0, _yh);
    syseq(t+dt/2.0, _yh, u, _fh);
    v_mltadd(_y1, _fh, dt/3.0, _y1);

    v_mltadd(y, _fh, dt/2.0, _yh);
    syseq(t+dt/2.0, _yh, u, _fh);
    v_mltadd(_y1, _fh, dt/3.0, _y1);

    v_mltadd(y, _fh, dt, _yh);
    syseq(t+dt, _yh, u, _fh);
    v_mltadd(_y1, _fh, dt/6.0, y);
  }
}


//========================================================================

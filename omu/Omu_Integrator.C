/*
 * Omu_Integrator.C --
 *   -- class definition
 *
 * rf, 10/2/96
 */

/*
    Copyright (C) 1997--2000  Ruediger Franke

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

#include "Omu_Integrator.h"

#include <If_Real.h>
#include <If_Int.h>
#include <If_Bool.h>

IF_BASE_DEFINE(Omu_Integrator);

//--------------------------------------------------------------------------
Omu_Integrator::Omu_Integrator()
{
  _kk = 0;
  _nxt = 0;
  _nd = 0;
  _nv = 0;
  _nx = 0;
  _nu = 0;
  _sa = false;
  _serr = false;
  _n = 0;
  _m = 0;
  _stepsize = 0.0;
  _rtol = 1e-8;
  _atol = 1e-8;
  _res_evals = 0;
  _sen_evals = 0;
  _jac_evals = 0;
  _ifList.append(new If_Bool("prg_int_serr", &_serr));
  _ifList.append(new If_Real("prg_int_stepsize", &_stepsize));
  _ifList.append(new If_Real("prg_int_rtol", &_rtol));
  _ifList.append(new If_Real("prg_int_atol", &_atol));
  _ifList.append(new If_Int("prg_int_res_evals", &_res_evals));
  _ifList.append(new If_Int("prg_int_sen_evals", &_sen_evals));
  _ifList.append(new If_Int("prg_int_jac_evals", &_jac_evals));
}

//--------------------------------------------------------------------------
Omu_Integrator::~Omu_Integrator()
{
}

//--------------------------------------------------------------------------
void Omu_Integrator::init_stage(int k,
				const Omu_States &x, const Omu_Vector &u,
				bool sa)
{
  _nxt = x->dim;
  _nd = x.nd;
  _nv = x.nv;
  _nx = _nxt - _nv;
  _nu = u->dim;
  _sa = sa;

  _n = _nxt - _nd;	// number of states for integration
  _m = _nd + _nu;	// number of parameters for integration
}

//--------------------------------------------------------------------------
void Omu_Integrator::init_sample(int kk, double tstart, double tend)
{
  _kk = kk;
}


//==========================================================================

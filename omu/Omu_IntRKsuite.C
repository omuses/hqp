/*
 * Omu_IntRKsuite.C --
 *   -- class implementation
 *
 * rf, 10/2/96
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

#include <adutils.h>

#include <If_Int.h>
#include <If_Class.h>

#include "Omu_IntRKsuite.h"

IF_CLASS_DEFINE("RKsuite", Omu_IntRKsuite, Omu_Integrator);

#if defined(HAVE_G2C_H) || defined(HAVE_F2C_H)
#  ifdef min
#    undef min
#  endif
#  ifdef max
#    undef max
#  endif
#  if defined(HAVE_G2C_H)
#    include <g2c.h>
#  else
#    include <f2c.h>
#  endif
#else
typedef long int integer;
typedef long int logical;
typedef double doublereal;
#endif
typedef doublereal dreal;

typedef void (f_t)(dreal *T, dreal *Y, dreal *YP);

extern "C" {
  void setup_(integer *NEQ, dreal *TSTART, dreal *YSTART, dreal *TEND,
	      dreal *TOL, dreal *THRES, integer *METHOD, char *TASK,
	      logical *ERRASS, dreal *HSTART, dreal *WORK, integer *LENWRK,
	      logical *MESAGE);
  void reset_(dreal *TENDNU);
  void ct_(f_t *F, dreal *TNOW, dreal *YNOW, dreal *YPNOW, dreal *WORK,
	   integer *CFLAG);
  void stat_(integer *TOTF, integer *STPCST, dreal *WASTE, integer *STPSOK,
	     dreal *HNEXT);
}

static Omu_IntRKsuite *theOmu_IntRKsuite;

//--------------------------------------------------------------------------
static void F(dreal *T, dreal *Y, dreal *YP)
{
  theOmu_IntRKsuite->F(T, Y, YP);
}

//--------------------------------------------------------------------------
Omu_IntRKsuite::Omu_IntRKsuite()
{
  theOmu_IntRKsuite = this;

  _y_head.max_dim = 0;
  _yp_head.max_dim = 0;
  _yp = v_get(1);
  _u = v_get(1);
  _work = v_get(1);
  _thres = v_get(1);
  _hnext = 0.0;
  _tlast = 0.0;
  _method = 2;

  _ifList.append(new If_Int("prg_int_method", &_method));
}

//--------------------------------------------------------------------------
Omu_IntRKsuite::~Omu_IntRKsuite()
{
  v_free(_thres);
  v_free(_work);
  v_free(_u);
  v_free(_yp);
}

//--------------------------------------------------------------------------
void Omu_IntRKsuite::ode_solve(double tstart, VECP y, const VECP u,
			       double tend)
{
  if (_sa)
    _neq = y->dim;
  else
    _neq = _n;
  _npar = u->dim;

  v_resize(_work, 32 * _neq);
  v_resize(_thres, _neq);
  v_resize(_u, _npar);
  v_resize(_yp, _neq);
  _y_head.dim = _neq;
  _yp_head.dim = _neq;

  v_set(_thres, _atol);

  // exclude sensitivities from error test
  if (!_serr) {
    for (int i = _n; i < (int)_thres->dim; i++)
      _thres[i] = 1e128;
  }

  v_zero(_yp);
  v_copy(u, _u);

  integer NEQ = _neq;

  // a hack to reinitialize _hnext in each simulation
  if (tstart < _tlast)
    _hnext = 0.0;
  _tlast = tstart;

  dreal TNOW = tstart;
  integer CFLAG;

  //integer METHOD = (_rtol > 5e-4)? 1: (_rtol > 5e-6)? 2: 3;
  integer METHOD = _method;
  char TASK = 'c';
  logical ERRASS = 0;
  dreal HSTART = _hnext;
  integer LENWRK = _work->dim;
  logical MESAGE = 0;

  setup_(&NEQ, &tstart, y->ve, &tend, &_rtol, _thres->ve,
	 &METHOD, &TASK, &ERRASS, &HSTART, _work->ve, &LENWRK, &MESAGE);

  while (TNOW < tend) {
    ct_(&::F, &TNOW, y->ve, _yp->ve, _work->ve, &CFLAG);
    if (CFLAG > 1)
      fprintf(stderr, "RKsuite message %d at time %g\n", CFLAG, TNOW);
    if (CFLAG > 4)
      m_error(E_UNKNOWN, "Omu_IntRKsuite::step");
  }

  integer TOTF, STPCST, STPSOK;
  dreal WASTE;
  stat_(&TOTF, &STPCST, &WASTE, &STPSOK, &HSTART);
  _hnext = HSTART;
}

//--------------------------------------------------------------------------
void Omu_IntRKsuite::F(double *T, double *Y, double *YP)
{
  _y_head.ve = Y;
  _yp_head.ve = YP;

  syseq(*T, &_y_head, _u, &_yp_head);
}


//========================================================================

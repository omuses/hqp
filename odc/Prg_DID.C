/*
 * Prg_DID.C -- class definition
 * (derived from Di example accessing Docp interface)
 *
 * rf, 05/05/01
 */

#include "Prg_DID.h"

#include <If_Bool.h>

IF_CLASS_DEFINE("DID", Prg_DID, Omu_Program);

//--------------------------------------------------------------------------
Prg_DID::Prg_DID()
{
  _K = 60;
  _nx = 2;
  _nu = 1;
  _with_cns = true;

  // interface elements
  _ifList.append(new If_Bool("prg_with_cns", &_with_cns));
}

//--------------------------------------------------------------------------
void Prg_DID::setup_stages(IVECP ks, VECP ts)
{
  stages_alloc(ks, ts, _K, 1, 0.0, 1.0);
}

//--------------------------------------------------------------------------
void Prg_DID::setup(int k,
		    Omu_Vector &x, Omu_Vector &u, Omu_Vector &c)
{
  x.alloc(_nx);
  if (k < _K)
    u.alloc(_nu);

  // initial values
  if (k == 0) {
    x.initial[0] = 1.0;
    x.initial[1] = 0.0;
  }
  if (k < _K)
    u.initial[0] = -2.0;

  // initial state constraints
  if (k == 0) {
    x.min[0] = x.max[0] = x.initial[0];
    x.min[1] = x.max[1] = x.initial[1];
  }
  // path constraint
  else if (k < _K) {
    x.max[1] = 0.01;
  }
  // final state constraints
  else {
    x.min[0] = x.max[0] = -1.0;
    x.min[1] = x.max[1] = 0.0;
  }

  // optional additional treatment for path constraint
  if (_with_cns) {
    if (k < _K) {
      c.alloc(1);
      c.max[0] = 0.01;
    }
  }
}

//--------------------------------------------------------------------------
void Prg_DID::update(int kk,
		     const adoublev &x, const adoublev &u,
		     adoublev &f, adouble &f0, adoublev &c)
{
  double dt = 1.0/_KK;

  // update constraints and objective for given x and u
  if (kk < _KK) {
    f[0] = x[0] + u[0]*dt;
    f[1] = x[1] + x[0]*dt + u[0]*0.5*dt*dt;

    f0 = u[0] * u[0] * dt;

    if (_with_cns) {
      c[0] = x[1] + 0.5*dt*x[0];
    }
  }
  else
    f0 = 0.0;
}


//==========================================================================

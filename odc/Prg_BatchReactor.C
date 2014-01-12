/*
 * Prg_BatchReactor.C -- class definition
 *
 * rf, 13/12/16
 */

#include "Prg_BatchReactor.h"

#include <If_Real.h>

IF_CLASS_DEFINE("BatchReactor", Prg_BatchReactor, Omu_Program);

//--------------------------------------------------------------------------
Prg_BatchReactor::Prg_BatchReactor()
{
  set_K(40); // can be modified via prg_K
  _kinf = 0.5;

  // interface elements
  _ifList.append(new If_Real("prg_kinf", &_kinf));
}

//--------------------------------------------------------------------------
void Prg_BatchReactor::setup_stages(IVECP ks, VECP ts)
{
  stages_alloc(ks, ts, K(), 1, 0.0, 1.0);
}

//--------------------------------------------------------------------------
void Prg_BatchReactor::setup(int k,
			     Omu_Vector &x, Omu_Vector &u, Omu_Vector &c)
{
  x.alloc(2);
  if (k < K())
    u.alloc(1);

  // initial control inputs and bounds
  if (k < K()) {
    u.initial[0] = 1.0;
    u.min[0] = 0.0;
    u.max[0] = 5.0;
  }
  // initial states and initial state constraints
  if (k == 0) {
    x.min[0] = x.max[0] = x.initial[0] = 1.0;
    x.min[1] = x.max[1] = x.initial[1] = 0.0;
  }
  // path constraint
  else {
    x.initial[0] = 0.5;
    x.initial[1] = 0.5;
    x.min[0] = x.min[1] = 0.0;
    x.max[1] = x.max[1] = 1.0;
  }
}

//--------------------------------------------------------------------------
void Prg_BatchReactor::update(int kk,
			      const adoublev &x, const adoublev &u,
			      adoublev &f, adouble &f0, adoublev &c)
{
  // update objective
  if (kk == KK())
    f0 = -x[1];
}

//--------------------------------------------------------------------------
void Prg_BatchReactor::continuous(int kk, double t, 
				  const adoublev &x, const adoublev &u,
				  const adoublev &dx, adoublev &F)
{
  F[0] = -(u[0] + _kinf*u[0]*u[0])*x[0] - dx[0];
  F[1] = u[0]*x[0] - dx[1];
}

//==========================================================================

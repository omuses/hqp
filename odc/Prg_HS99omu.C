/*
 * Prg_HS99omu.C -- class definition
 *
 * rf, 1/13/97
 */


#include "Prg_HS99omu.h"

#include <If_Class.h>

// propagate the class to the command interface
IF_CLASS_DEFINE("HS99omu", Prg_HS99omu, Omu_Program);

//--------------------------------------------------------------------------
static double a[] = {0.0, 50.0, 50.0, 75.0, 75.0, 75.0, 100.0, 100.0};
static double b = 32.0;
static double t[] = {0.0, 25.0, 50.0, 100.0, 150.0, 200.0, 290.0, 380.0};

//--------------------------------------------------------------------------
void Prg_HS99omu::setup_stages(IVECP ks, VECP ts)
{
  // allocate seven stages with one sample period per stage
  stages_alloc(ks, ts, 7, 1);

  // initialize communication time points
  for (int i = 0; i <= KK(); i++) {
    ts[i] = ::t[i];
  }
}

//--------------------------------------------------------------------------
void Prg_HS99omu::setup(int k,
			Omu_Vector &x, Omu_Vector &u, Omu_Vector &c)
{
  x.alloc(3);
  if (k < K())
    u.alloc(1);

  // bounds on variables and initial values
  if (k < K()) {
    u.min[0] = 0.0;
    u.max[0] = 1.58;
    u.initial[0] = 0.5;
  }

  // initial and final state constraints
  if (k == 0) {
    for (int i = 0; i < 3; i++)
      x.min[i] = x.max[i] = 0.0;
  }
  else if (k == K()) {
    x.min[1] = x.max[1] = 1e5;	// q(t_K)
    x.min[2] = x.max[2] = 1e3;	// s(t_K)
  }
}

//--------------------------------------------------------------------------
void Prg_HS99omu::update(int kk,
			 const adoublev &x, const adoublev &u,
			 adoublev &f, adouble &f0, adoublev &c)
{
  if (kk == KK()) {
    f0 = -x[0]*x[0];	// -r(t_K)^2
  }
}

//--------------------------------------------------------------------------
void Prg_HS99omu::continuous(int kk, double t,
			     const adoublev &x, const adoublev &u,
			     const adoublev &xp, adoublev &F)
{
  F[0] = ::a[kk+1] * cos(u[0]);
  F[1] = x[2];	// s(t)
  F[2] = ::a[kk+1] * sin(u[0]) - ::b;

  F -= xp;
}

//==========================================================================

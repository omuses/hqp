/*
 * Prg_TP383omu.C -- class definition
 *
 * rf, 1/13/97
 */

#include "Prg_TP383omu.h"

#include <If_Class.h>

// propagate the class to the command interface
IF_CLASS_DEFINE("TP383omu", Prg_TP383omu, Omu_Program);

//--------------------------------------------------------------------------
static double a[] = {
  12842.275, 634.25, 634.25, 634.125, 1268.0, 633.875, 633.75,
  1267.0, 760.05, 633.25, 1266.25, 632.875, 394.46, 940.838
};
static double c[] = {
  5.47934, 0.83234, 0.94749, 1.11082, 2.64824, 1.55868, 1.73215,
  3.90896, 2.74284, 2.60541, 5.96184, 3.29522, 1.83517, 2.81372
};

//--------------------------------------------------------------------------
void Prg_TP383omu::setup_stages(IVECP ks, VECP ts)
{
  // setup a problem with 14 stages and one sample period per stage
  // stages_alloc() initializes ks and ts and sets K = KK = 14
  stages_alloc(ks, ts, 14, 1);
}

//--------------------------------------------------------------------------
void Prg_TP383omu::setup(int k,
			 Omu_Vector &x, Omu_Vector &u, Omu_Vector &)
{
  x.alloc(1);
  if (k < K())
    u.alloc(1);

  // bounds on variables and initial values
  if (k < 5) {
    u.min[0] = 0.0;
    u.max[0] = 0.04;
    u.initial[0] = 0.01;
  }
  else if (k < K()) {
    u.min[0] = 0.0;
    u.max[0] = 0.03;
    u.initial[0] = 0.01;
  }

  // initial and final state constraints
  if (k == 0) {
    x.min[0] = x.max[0] = 0.0;	// s^0
    x.initial[0] = 0.0;
  }
  else if (k == K()) {
    x.min[0] = x.max[0] = 1.0;	// s^K
  }
}

//--------------------------------------------------------------------------
void Prg_TP383omu::update(int kk,
			  const adoublev &x, const adoublev &u,
			  adoublev &f, adouble &f0, adoublev &)
{
  if (kk < KK()) {
    f0 = ::a[kk] / u[0];
    f[0] = x[0] + ::c[kk] * u[0];
  }
}


//==========================================================================

/*
 * Prg_TP383.C -- class definition
 *
 * rf, 1/13/97
 */

#include "Prg_TP383.h"

#include <If_Class.h>

// propagate the class to the command interface
IF_CLASS_DEFINE("TP383", Prg_TP383, Omu_Program);

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
void Prg_TP383::setup(int k,
		      Omu_Vector &x, Omu_Vector &u, Omu_Vector &cns)
{
  // allocate optimization variables
  x.alloc(14);

  // bounds on variables and initial values
  for (int i = 0; i < 5; i++) {
    x.min[i] = 0.0;
    x.max[i] = 0.04;
    x.initial[i] = 0.01;
  }
  for (int i = 5; i < 14; i++) {
    x.min[i] = 0.0;
    x.max[i] = 0.03;
    x.initial[i] = 0.01;
  }

  // one equality constraint
  cns.alloc(1);
  cns.min[0] = cns.max[0] = 1.0;
}

//--------------------------------------------------------------------------
void Prg_TP383::update(int kk,
		       const adoublev &x, const adoublev &u,
		       adoublev &f, adouble &f0, adoublev &cns)
{
  for (int i = 0; i < 14; i++) {
    f0 += ::a[i] / x[i];
    cns[0] += ::c[i] * x[i];
  }
}


//==========================================================================

/*
 * Prg_HS99.C -- class definition
 *
 * rf, 1/13/97
 */

#include <math.h>

#include <If_Class.h>

#include "Prg_HS99.h"

// propagate the class to the command interface
IF_CLASS_DEFINE("HS99", Prg_HS99, Omu_Program);

//--------------------------------------------------------------------------
static double a[] = {0.0, 50.0, 50.0, 75.0, 75.0, 75.0, 100.0, 100.0};
static double b = 32.0;
static double t[] = {0.0, 25.0, 50.0, 100.0, 150.0, 200.0, 290.0, 380.0};

//--------------------------------------------------------------------------
void Prg_HS99::setup(int k,
		     Omu_Vector &x, Omu_Vector &u, Omu_Vector &c)
{
  // allocate optimization variables
  x.alloc(7);

  // bounds and initial values
  for (int i = 0; i < 7; i++) {
    x.min[i] = 0.0;
    x.max[i] = 1.58;
    x.initial[i] = 0.5;
  }

  // two equality constraints for final state constraints
  c.alloc(2);
  c.min[0] = c.max[0] = 1e5;
  c.min[1] = c.max[1] = 1e3;
}

//--------------------------------------------------------------------------
void Prg_HS99::update(int kk, 
		      const adoublev &x, const adoublev &u,
		      adoublev &f, adouble &f0, adoublev &c)
{
  adouble r, q, s, p;

  // initial states
  r = 0.0;
  q = 0.0;
  s = 0.0;

  // system integration
  for (int i = 1; i < 8; i++) {
    r = r + ::a[i] * cos(x[i-1]) * (::t[i] - ::t[i-1]);
    p = (::a[i] * sin(x[i-1]) - ::b) * (::t[i] - ::t[i-1]);
    q = q + (0.5 * p + s) * (::t[i] - ::t[i-1]);
    s = s + p;
  }

  f0 = -r*r;

  c[0] = q;
  c[1] = s;
}


//==========================================================================

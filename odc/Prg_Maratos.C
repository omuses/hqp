/*
 * Prg_Maratos.C -- class definition
 *
 * rf, 1/12/97
 */

#include <If_Class.h>

#include "Prg_Maratos.h"

IF_CLASS_DEFINE("Maratos", Prg_Maratos, Omu_Program);

//--------------------------------------------------------------------------
void Prg_Maratos::setup(int k,
			Omu_Vector &x, Omu_Vector &u, Omu_Vector &c)

{
  x.alloc(2);
  x.initial[0] = 0.8;
  x.initial[1] = 0.6;

  c.alloc(1);
  c.min[0] = c.max[0] = 0.0;
}

//--------------------------------------------------------------------------
void Prg_Maratos::update(int kk, 
			 const adoublev &x, const adoublev &u,
			 adoublev &f, adouble &f0, adoublev &c)
{
  adouble x1, x2;

  x1 = x[0];
  x2 = x[1];
  f0 = -x1 + 10.0*(x1*x1 + x2*x2 - 1.0);
  c[0] = x1*x1 + x2*x2 - 1.0;
}


//==========================================================================

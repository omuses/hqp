/*
 * Prg_Maratos.h -- 
 *   - Maratos problem
 *
 * rf, 1/12/97
 */

#ifndef Prg_Maratos_H
#define Prg_Maratos_H

#include <Omu_Program.h>

//--------------------------------------------------------------------------
class Prg_Maratos: public Omu_Program {

 protected:

  void setup(int k,
	     Omu_Vector &x, Omu_Vector &u, Omu_Vector &c);

  void update(int kk, 
	      const adoublev &x, const adoublev &u,
	      adoublev &f, adouble &f0, adoublev &c);

 public:

  char *name() {return "Maratos";}
};  

#endif


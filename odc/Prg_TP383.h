/*
 * Prg_TP383.h -- 
 *   - test example TP383
 *
 * rf, 1/13/97
 */

#ifndef Prg_TP383_H
#define Prg_TP383_H

#include <Omu_Program.h>

//--------------------------------------------------------------------------
class Prg_TP383: public Omu_Program {

 protected:

  void setup(int k,
	     Omu_Vector &x, Omu_Vector &u, Omu_Vector &c);

  void update(int kk,
	      const adoublev &x, const adoublev &u,
	      adoublev &f, adouble &f0, adoublev &c);

 public:

  char *name() {return "TP383";}
};  

#endif


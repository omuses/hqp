/*
 * Prg_TP383omu.h -- 
 *   - test example TP383, multistage formulation
 *
 * rf, 1/13/97
 */

#ifndef Prg_TP383omu_H
#define Prg_TP383omu_H

#include <Omu_Program.h>

//--------------------------------------------------------------------------
class Prg_TP383omu: public Omu_Program {

 protected:

  void setup_stages(IVECP ks, VECP ts);

  void setup(int k,
	     Omu_Vector &x, Omu_Vector &u, Omu_Vector &c);

  void update(int kk, 
	      const adoublev &x, const adoublev &u,
	      adoublev &f, adouble &f0, adoublev &c);

 public:

  char *name() {return "TP383omu";}
};  

#endif


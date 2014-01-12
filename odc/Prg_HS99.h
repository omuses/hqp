/*
 * Prg_HS99.h -- 
 *   - test example HS99
 *
 * rf, 1/13/97
 */

#ifndef Prg_HS99_H
#define Prg_HS99_H

#include <Omu_Program.h>

//--------------------------------------------------------------------------
class Prg_HS99: public Omu_Program {

 protected:

  void setup(int k,
	     Omu_Vector &x, Omu_Vector &u, Omu_Vector &c);

  void update(int kk, 
	      const adoublev &x, const adoublev &u,
	      adoublev &f, adouble &f0, adoublev &c);

 public:

  const char *name() {return "HS99";}
};  

#endif


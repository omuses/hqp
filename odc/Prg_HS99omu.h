/*
 * Prg_HS99omu.h -- 
 *   - test example HS99 exploiting Omuses
 *
 * rf, 1/13/97
 */

#ifndef Prg_HS99omu_H
#define Prg_HS99omu_H

#include <Omu_Program.h>

//--------------------------------------------------------------------------
class Prg_HS99omu: public Omu_Program {

 protected:

  void setup_stages(IVECP ks, VECP ts);

  void setup(int k,
	     Omu_Vector &x, Omu_Vector &u, Omu_Vector &c);

  void update(int kk,
	      const adoublev &x, const adoublev &u,
	      adoublev &f, adouble &f0, adoublev &c);

  void continuous(int kk, double t,
		  const adoublev &x, const adoublev &u,
		  const adoublev &xp, adoublev &F);

 public:

  char *name() {return "HS99omu";}
};  

#endif


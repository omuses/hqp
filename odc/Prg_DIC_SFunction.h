/*
 * Prg_DIC_SFunction.h -- 
 *   - double integrator example accessing continuous-time Simulink S-function
 *     (derived from Di example accessing Docp interface)
 *
 * rf, 05/06/01
 */

#ifndef Prg_DIC_SFunction_H
#define Prg_DIC_SFunction_H

#include <Omu_Program.h>

#include "simstruc.h"

//--------------------------------------------------------------------------
class Prg_DIC_SFunction: public Omu_Program {

 protected:

  SimStruct *_S;

  void setup_stages(IVECP ks, VECP ts);

  void setup(int k,
	     Omu_Vector &x, Omu_Vector &u, Omu_Vector &c);

  void update(int kk,
	      const adoublev &x, const adoublev &u,
	      adoublev &f, adouble &f0, adoublev &c);

  void continuous(int kk, double t,
		  const adoublev &x, const adoublev &u, const adoublev &xp,
		  adoublev &F);

  int _nx;	// number of states
  int _nu;	// number of control parameters
  int _mdl_ny;	// number of model outputs
  double _dt; 	// sample time
  bool _with_cns;// treat overshoot with additional constraint

 public:

  Prg_DIC_SFunction();
  ~Prg_DIC_SFunction();

  char *name() {return "DIC_SFunction";}
};  

#endif


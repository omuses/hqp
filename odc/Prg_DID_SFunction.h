/*
 * Prg_DID_SFunction.h -- 
 *   - double integrator example accessing discrete-time Simulink S-function
 *     (derived from double integrator example accessing Docp interface)
 *
 * rf, 05/05/01
 */

#ifndef Prg_DID_SFunction_H
#define Prg_DID_SFunction_H

#include <Omu_Program.h>

#include "simstruc.h"

//--------------------------------------------------------------------------
class Prg_DID_SFunction: public Omu_Program {

 protected:

  SimStruct *_S;

  void setup_stages(IVECP ks, VECP ts);

  void setup(int k,
	     Omu_Vector &x, Omu_Vector &u, Omu_Vector &c);

  void update(int kk,
	      const adoublev &x, const adoublev &u,
	      adoublev &f, adouble &f0, adoublev &c);

  int _nx;	// number of states
  int _nu;	// number of control parameters
  int _mdl_ny;	// number of model outputs
  double _dt; 	// sample time
  adouble _adt;	// active sample time parameter
  bool _with_cns;// treat overshoot with additional constraint

 public:

  Prg_DID_SFunction();
  ~Prg_DID_SFunction();

  char *name() {return "DID_SFunction";}
};  

#endif


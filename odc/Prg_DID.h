/*
 * Prg_DID.h -- 
 *   - double integrator example, discrete-time
 *     (derived from double integrator example accessing Docp interface)
 *
 * rf, 05/05/01
 */

#ifndef Prg_DID_H
#define Prg_DID_H

#include <Omu_Program.h>

//--------------------------------------------------------------------------
class Prg_DID: public Omu_Program {

 protected:

  void setup_stages(IVECP ks, VECP ts);

  void setup(int k,
	     Omu_Vector &x, Omu_Vector &u, Omu_Vector &c);

  void update(int kk,
	      const adoublev &x, const adoublev &u,
	      adoublev &f, adouble &f0, adoublev &c);

  int _nx;	// number of states
  int _nu;	// number of control parameters
  bool _with_cns;// treat overshoot with additional constraint

 public:

  Prg_DID();

  const char *name() {return "DID";}
};  

#endif


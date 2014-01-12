/*
 * Prg_BatchReactor.h -- 
 *   Non-isothermal batch reactor example with bounds.
 *   See: Biegler: Nonlinear Programming Concepts, Algorithms and
 *        Applications to Chemical Processes, MOS-SIAM 2010
 *   Given first order parallel reactions: A -> B, A -> C
 *     da/dt = -k10*exp(-E1/RT)*a(t) - k20*exp(-E2/RT)*a(t)
 *     db/dt =  k10*exp(-E1/RT)*a(t)
 *   Task: starting from a(0)=1, b(0)=0, find temperature profile such that
 *   the amount of product B is maximized at the final time; b(tf)->max.
 *
 * rf, 13/12/16
 */

#ifndef Prg_BatchReactor_H
#define Prg_BatchReactor_H

#include <Omu_Program.h>

//--------------------------------------------------------------------------
class Prg_BatchReactor: public Omu_Program {

 protected:

  double _kinf; // ratio of pre-exponential factors k20/k10

  void setup_stages(IVECP ks, VECP ts);

  void setup(int k,
	     Omu_Vector &x, Omu_Vector &u, Omu_Vector &c);

  void update(int kk,
	      const adoublev &x, const adoublev &u,
	      adoublev &f, adouble &f0, adoublev &c);

  void continuous(int kk, double t, 
                  const adoublev &x, const adoublev &u,
                  const adoublev &dx, adoublev &F);

 public:

  Prg_BatchReactor();

  const char *name() {return "BatchReactor";}
};  

#endif


/*
 * Prg_BatchReactor_bare.h -- 
 *   Non-isothermal batch reactor example with bounds.
 *   Alternative implementation using low-level interface without ADOL-C.
 *
 * rf, 13/12/16
 */

#ifndef Prg_BatchReactor_bare_H
#define Prg_BatchReactor_bare_H

#include <Omu_Program.h>

//--------------------------------------------------------------------------
class Prg_BatchReactor_bare: public Omu_Program {

 protected:

  double _kinf; // ratio of pre-exponential factors k20/k10

  void setup_stages(IVECP ks, VECP ts);

  void setup_struct(int k,
                    const Omu_VariableVec &x,
                    const Omu_VariableVec &u,
                    Omu_DependentVec &xt, Omu_DependentVec &F,
                    Omu_DependentVec &f,
                    Omu_Dependent &f0, Omu_DependentVec &c);
  void setup(int k,
	     Omu_Vector &x, Omu_Vector &u, Omu_Vector &c);

  void consistic(int kk, double t,
                 const Omu_StateVec &x, const Omu_Vec &u,
                 Omu_DependentVec &xt);

  void continuous(int kk, double t,
                  const Omu_StateVec &x, const Omu_Vec &u,
                  const Omu_StateVec &dx, Omu_DependentVec &F);

  void update(int kk,
              const Omu_StateVec &x, const Omu_Vec &u,
              const Omu_StateVec &xf,
              Omu_DependentVec &f, Omu_Dependent &f0,
              Omu_DependentVec &c);

 public:

  Prg_BatchReactor_bare();

  const char *name() {return "BatchReactor_bare";}
};  

#endif


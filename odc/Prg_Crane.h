/*
 * Prg_Crane.h -- 
 *   - crane according to the laboratory Opt3
 *
 * rf, 1/15/97
 */

#ifndef Prg_Crane_H
#define Prg_Crane_H

#include <Omu_Program.h>

//--------------------------------------------------------------------------
class Prg_Crane: public Omu_Program {

 protected:

  double _tf_guess;	// initial guess for the final time
  double _u_bound;	// bound on the control variable
  double _phi_bound;	// bound on the angle state
  bool 	 _bound_init;	// apply bound constraints to initial states

  void setup_stages(IVECP ks, VECP ts);

  void setup(int k,
	     Omu_Vector &x, Omu_Vector &u, Omu_Vector &c);

  void init_simulation(int k,
		       Omu_Vector &x, Omu_Vector &u);

  void update(int kk,
	      const adoublev &x, const adoublev &u,
	      adoublev &f, adouble &f0, adoublev &c);

  void continuous(int kk, double t,
		  const adoublev &x, const adoublev &u, const adoublev &xp,
		  adoublev &F);

  void model_eq(int kk, double t,
		const adoublev &x, const adoublev &u,
		adoublev &xp);

  int _nx;	// number of states
  int _nu;	// number of control parameters
  int offs;	// offset for the automatically generated equations

  // model inputs and parameters
  double Fscale;
  double g;
  double l;
  double md;
  double mdl;
  double ml;

 public:

  Prg_Crane();

  const char *name() {return "Crane";}
};  

#endif


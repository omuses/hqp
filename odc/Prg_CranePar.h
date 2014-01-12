/*
 * Prg_CranePar.h -- 
 *   - parameter/initial states estimation example
 *
 * rf, 1/15/97
 */

#ifndef Prg_CranePar_H
#define Prg_CranePar_H

#include <Omu_Program.h>
#include <If_Command.h>

//--------------------------------------------------------------------------
class Prg_CranePar: public Omu_Program {

 protected:

  VECP   _s_ref;
  double _maxdev;
  int	 _seed;
  bool	 _multistage;

  void setup_stages(IVECP ks, VECP ts);

  void setup(int k,
	     Omu_Vector &x, Omu_Vector &u, Omu_Vector &c);

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
  int offs;	// offset for the automatically generated equations

  // model inputs and parameters
  double Fscale;
  double g;
  double l;
  double md;
  double mdl;
  double ml;

 public:

  Prg_CranePar();
  ~Prg_CranePar();

  int disturb(IF_DEF_ARGS);		// method to be called trough Tcl

  const char *name() {return "CranePar";}
};  

#endif


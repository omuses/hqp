/*
 * Prg_CranePar.C -- class definition
 *
 * rf, 1/15/97
 */

#include <math.h>

#include "Prg_CranePar.h"

#include <If_RealVec.h>
#include <If_Real.h>
#include <If_Bool.h>
#include <If_Int.h>
#include <If_Method.h>

IF_CLASS_DEFINE("CranePar", Prg_CranePar, Omu_Program);

typedef If_Method<Prg_CranePar> If_Cmd;

//--------------------------------------------------------------------------
Prg_CranePar::Prg_CranePar()
{
  _nx = 5;	// 4 states plus the mass parameter
  offs = 1;	// for the container mass parameter
  _KK = 1;

  _maxdev = 0.05;
  _seed = 1234;
  _s_ref = v_get(_KK + 1);
  _multistage = true;

  _ifList.append(new If_RealVec("prg_s_ref", &_s_ref));
  _ifList.append(new If_Int("prg_seed", &_seed));
  _ifList.append(new If_Real("prg_maxdev", &_maxdev));
  _ifList.append(new If_Bool("prg_multistage", &_multistage));
  _ifList.append(new If_Cmd("prg_disturb", &Prg_CranePar::disturb, this));

  // default values for parameters
  Fscale = 1000.0;
  g = 9.81;
  l = 10.0;
  md = 1000.0;
  ml = 4000.0;

  // interface elements for unbound variables
  _ifList.append(new If_Real("prg_Fscale", &Fscale));
  _ifList.append(new If_Real("prg_g", &g));
  _ifList.append(new If_Real("prg_l", &l));
  _ifList.append(new If_Real("prg_md", &md));
  _ifList.append(new If_Real("prg_ml", &ml));
}

//--------------------------------------------------------------------------
Prg_CranePar::~Prg_CranePar()
{
  v_free(_s_ref);
}

//--------------------------------------------------------------------------
void Prg_CranePar::setup_stages(IVECP ks, VECP ts)
{
  if (_multistage) {
    _K = _KK;
    stages_alloc(ks, ts, _K, 1);
  }
  else {
    _K = 1;
    stages_alloc(ks, ts, 1, _KK);
  }

  // resize the vectors for measurement data
  v_resize(_s_ref, _KK + 1);
}

//--------------------------------------------------------------------------
void Prg_CranePar::setup(int k,
			 Omu_Vector &x, Omu_Vector &u, Omu_Vector &c)
{
  // allocate and initialize the free variables
  if (k == 0) {
    x.alloc(_nx);

    // bound parameters and implicit discrete part
    // (this calculation is redefined with adoubles in model_eq())
    mdl = md + ml;

    // initial guess for the mass
    x.initial[0] = ml / 1000.0;

    // initial states
    x.initial[offs+0] = 0.0;	// phi
    x.initial[offs+1] = 0.0;	// omega
    x.initial[offs+2] = 0.0;	// v
    x.initial[offs+3] = 25.0;	// s
  }

  // allocate additional state variables in the multistage case
  else if (_multistage) {
    x.alloc(_nx);
  }
}

//--------------------------------------------------------------------------
int Prg_CranePar::disturb(IF_CMD_ARGS)
{
  // disturb the recorded data once during the problem setup
  int i, i_end, r;
  i_end = _s_ref->dim;

  srand((unsigned int)_seed);
  for (i = 0; i < i_end; i++) {
    _s_ref[i] += _maxdev * ((double)rand() / (double)RAND_MAX * 2.0 - 1.0);
  }

  return IF_OK;
}

//--------------------------------------------------------------------------
void Prg_CranePar::update(int kk,
			  const adoublev &x, const adoublev &u, 
			  adoublev &f, adouble &f0, adoublev &c)
{
  adouble s;

  // discrete-time state equation for the mass parameter
  if (kk < _KK)
    f[0] = x[0];

  if (_multistage) {
    //
    // We have an optimization problem with one stage for each measurement
    // point (_K == _KK). This results in 5*(_KK+1) optimization variables.
    // The criterion is defined as function of the optimization variable
    // that describes the state s in each stage.
    // 
    s = x[offs+3];
    f0 = pow(s - _s_ref[kk], 2);
  }
  else {
    //
    // Alternatively we can formulate an optimization problem with one
    // stage (_K == 1) and without final states. This results in 5 
    // optimization variables. The initial state value of s is used
    // to formulate the criterion for the first measurement point. 
    // Afterwards the intermediate results for each sample period are used.
    //
    if (kk == 0) {
      s = x[offs+3];
      f0 = pow(s - _s_ref[kk], 2);
      s = f[offs+3];
      f0 += pow(s - _s_ref[kk+1], 2);
    }
    else if (kk < _KK) {
      s = f[offs+3];
      f0 = pow(s - _s_ref[kk+1], 2);
    }
  }
}

//--------------------------------------------------------------------------
void Prg_CranePar::continuous(int kk, double t,
			      const adoublev &x, const adoublev &u,
			      const adoublev &xp, adoublev &F)
{
  model_eq(kk, t, x, u, F);
  for (int i = offs; i < _nx; i++)
    F[i] -= xp[i];
}

//--------------------------------------------------------------------------
void Prg_CranePar::model_eq(int kk, double t,
			    const adoublev &x, const adoublev &u,
			    adoublev &xp)
{
  // dynamic model variables
  adouble den, omega, phi, s, sinphi, u_control, v;

  // redefine the mass parameter and its descendants as adouble
  adouble ml, mdl;

  ml = 1000.0 * x[0];

  // recalculate the bound parameters and implicit discrete part
  mdl = md + ml;

  // state assignments
  phi = x[offs+0];
  omega = x[offs+1];
  v = x[offs+2];
  s = x[offs+3];

  // constant control
  u_control = -1.0;

  // dynamic model equations
  //u_control = u[0];
  sinphi = sin(phi);
  den = md + ml*pow(sinphi, 2);
  xp[offs+0] = omega;
  xp[offs+1] = -(mdl*g*sinphi + 0.5*ml*l*pow(omega, 2)*sin(2*phi) + u_control*Fscale*cos(phi))/(l*den);
  xp[offs+2] = (0.5*ml*g*sin(2*phi) + ml*l*pow(omega, 2)*sinphi + u_control*Fscale)/den;
  xp[offs+3] = v;
}

//==========================================================================

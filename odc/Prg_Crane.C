/*
 * Prg_Crane.C -- class definition
 *
 * rf, 7/15/96
 */

#include <math.h>

#include "Prg_Crane.h"

#include <If_Bool.h>
#include <If_Real.h>

IF_CLASS_DEFINE("Crane", Prg_Crane, Omu_Program);

//--------------------------------------------------------------------------
Prg_Crane::Prg_Crane()
{
  set_K(50);
  _tf_guess = 15.0;
  _u_bound = 5.0;
  _phi_bound = 5.0 / 180.0 * 3.14159;
  _bound_init = false;

  _nx = 6;
  _nu = 1;
  offs = 1;	// for the final time parameter

  _ifList.append(new If_Real("prg_tf_guess", &_tf_guess));
  _ifList.append(new If_Real("prg_u_bound", &_u_bound));
  _ifList.append(new If_Real("prg_phi_bound", &_phi_bound));
  _ifList.append(new If_Bool("prg_bound_init", &_bound_init));

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
void Prg_Crane::setup_stages(IVECP ks, VECP ts)
{
  stages_alloc(ks, ts, K(), 1, 0.0, 1.0);
}

//--------------------------------------------------------------------------
void Prg_Crane::setup(int k,
		      Omu_Vector &x, Omu_Vector &u, Omu_Vector &c)
{
  double u_guess;

  x.alloc(_nx);
  if (k < K())
    u.alloc(_nu);

  if (k == 0) {
    // bound parameters and implicit discrete part
    mdl = md + ml;

    // initial state constraints
    x.min[offs+0] = x.max[offs+0] = 0.0;	// phi
    x.min[offs+1] = x.max[offs+1] = 0.0;	// omega
    x.min[offs+2] = x.max[offs+2] = 0.0;	// v
    x.min[offs+3] = x.max[offs+3] = 25.0;	// s

    // initial states
    x.initial[offs+0] = 0.0;	// phi
    x.initial[offs+1] = 0.0;	// omega
    x.initial[offs+2] = 0.0;	// v
    x.initial[offs+3] = 25.0;	// s
  }
  else if (k == K()) {
    // final state constraints
    x.min[offs+0] = x.max[offs+0] = 0.0;	// phi
    x.min[offs+1] = x.max[offs+1] = 0.0;	// omega
    x.min[offs+2] = x.max[offs+2] = 0.0;	// v
    x.min[offs+3] = x.max[offs+3] = 0.0;	// s
  }
  else {
    // bounds for phi
    x.min[offs+0] = -_phi_bound;
    x.max[offs+0] =  _phi_bound;

    // bounds for s
    x.min[offs+3] =  0.0;
    x.max[offs+3] = 25.0;
  }

  // lower bound on the final time
  x.min[0] = 1.0;

  // control bounds
  x.min[5] = -_u_bound;
  x.max[5] = _u_bound;

  //
  // init solution
  //

  x.initial[0] = _tf_guess;

  if (k < K()/2)
    u_guess = -100.0 * mdl / Fscale / (_tf_guess*_tf_guess);
  else
    u_guess = 100.0 * mdl / Fscale / (_tf_guess*_tf_guess);

  x.initial[5] = u_guess;

  if (k < K()) {
    if (k <= K()/2 && k+1 > K()/2)
      u.initial[0] = 2.0 * u_guess / (_tf_guess / (double)K());
    else
      u.initial[0] = 0.0;
  }
}

//--------------------------------------------------------------------------
void Prg_Crane::init_simulation(int k,
				Omu_Vector &x, Omu_Vector &u)
{
  int i;

  // set initial states for the first stage;
  // afterwards simulation results of the preceding stage are used
  if (k == 0) {
    for (i = 0; i < _nx; i++)
      x[i] = x.initial[i];
  }

  // set the initial control for each stage
  if (k < K())
    u[0] = u.initial[0];

  // apply bounds to simulated states
  if (_bound_init) {
    for (i = 0; i < _nx; i++) {
      x[i] = min(x[i], x.max[i]);
      x[i] = max(x[i], x.min[i]);
    }
  }
}

//--------------------------------------------------------------------------
void Prg_Crane::update(int kk,
		       const adoublev &x, const adoublev &u,
		       adoublev &f, adouble &f0, adoublev &c)
{
  if (kk < KK())
    f[0] = x[0];	// constant final time
  else
    f0 = x[0];		// optimization criterion
}

//--------------------------------------------------------------------------
void Prg_Crane::continuous(int kk, double t,
			   const adoublev &x, const adoublev &u,  
			   const adoublev &xp, adoublev &F)
{
  adouble tscale;

  tscale = x[0];
  model_eq(kk, value(tscale) * t, x, u, F);
  for (int i = offs; i < _nx; i++) {
    F[i] = tscale * F[i] - xp[i];
  }
}

//--------------------------------------------------------------------------
void Prg_Crane::model_eq(int kk, double t,
			 const adoublev &x, const adoublev &u,
			 adoublev &xp)
{
  // dynamic model variables
  adouble den, omega, phi, s, sinphi, u_control, v;

  // state assignments
  phi = x[offs+0];
  omega = x[offs+1];
  v = x[offs+2];
  s = x[offs+3];

  // piecewise linear control
  xp[5] = u[0];
  u_control = x[5];

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

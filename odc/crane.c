    // bound parameters and implicit discrete part
    mdl = md + ml;

  // dynamic model equations
  u_control = u[0];
  sinphi = sin(phi);
  den = md + ml*pow(sinphi, 2);
  xp[offs+0] = omega;
  xp[offs+1] = -(mdl*g*sinphi + 0.5*ml*l*pow(omega, 2)*sin(2*phi) + u_control*Fscale*cos(phi))/(l*den);
  xp[offs+2] = (0.5*ml*g*sin(2*phi) + ml*l*pow(omega, 2)*sinphi + u_control*Fscale)/den;
  xp[offs+3] = v;

  // state assignments
  phi = x[offs+0];
  omega = x[offs+1];
  v = x[offs+2];
  s = x[offs+3];

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

  // model inputs and parameters
  double Fscale;
  double g;
  double l;
  double md;
  double mdl;
  double ml;

  // dynamic model variables
  adouble den, omega, phi, s, sinphi, u_control, v;

  # interface elements for unresolved variables

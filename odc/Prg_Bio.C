//-----------------------------------------------------------------------------
// Prg_Bio.C -- class definition
//   - Hqp/Omuses example: bio-technological process
//
// Adapted from:
//    @PhdThesis{Pfaff:91,
//       author = {Pfaff, M.},
//       title =  {Entwurf optimaler {S}teuerstrategien ausgew{\"a}hlter 
//                 biotechnologischer {P}rozesse anhand aggregierter 
//                 kinetischer {M}odelle},
//       school = {Technische Hochschule Ilmenau},
//       year =   1991}
//
// E. Arnold 10/20/96
//           2003-02-16 odc example
//
//-----------------------------------------------------------------------------
// A bio-technological fermentation process is considered.
// The continuous fed-batch process is described by a second
// order nonlinear state equation. State variables are the amount (mass) of
// product and the added substrate, resp.. The only control variable
// is the substrate inflow rate. The production term Pi(.) on 
// the right hand side of the first state equation depends on the product 
// concentration Cp and substrate concentration Cs, which both are 
// (nonlinear) functions of the state variables.
// The aim is to minimize over a fixed time horizon a functional which 
// takes into account the cost of the added substrate and the profit 
// obtained from the product that depends on the product mass and product
// concentration at final time.
//
// State equation:     x1' = Pi(Cp,Cs,t) with Cp = Cp(x1,x2), Cs = Cs(x1,x2)
//                     x2' = u
// Initial state:      x(0) = [ p0 0 ]'   (fixed)
// Terminal state:     x(T) free, T fixed
// Cost functional:    J = c1*x1(T) + c2*x2(T) + c3 ---> min
//                     Constants c1, c2, c3 result from real prices of
//                     product and substrate.
// Control bound:      0 <= u <= umax
//-----------------------------------------------------------------------------

#include <math.h>

#include <If_Real.h>
#include <If_Int.h>
#include <If_Bool.h>
#include <If_RealVec.h>

#include "Prg_Bio.h"

IF_CLASS_DEFINE("Bio", Prg_Bio, Omu_Program);


//--------------------------------------------------------------------------
Prg_Bio::Prg_Bio()
{
    // number of stages
    set_K(51);
    
    // initialize local variables
    _fscale = 1.0;   // cost function scaling parameter
    _y = VNULL;
    _controller = 1; // use the nonlinear state controller for generation of
                     // the initial control trajectory (0/1)
    _uinit = 0.01;   // initial control level without the nonlinear controller
    
    // model parameters
    pimax = 0.16; // various kinetic and stochiometric parameters    
    ks = 1;        
    kis = 160;    
    kip = 75;      
    kd = 0.006;
    yps = 0.55;     
    kappa = 600;   
    cdos = 750;   
    kp = 0.08;     
    kap = 0.1;
    kos = 0.02;     
    cs0 = 5.0;    // initial value substrate concentration 
    v0 = 5;       // initial value volume 
    p0 = 0*v0;    // initial value product mass
    x0 = 30*v0;   // bio mass (constant within optimization horizon)
    Fsmin = 0;    // minimum substrate flow rate  
    Fsmax = 0.1;  // maximum substrate flow rate   
    tf = 10;      // optimization horizon

    // interface elements configurable via Tcl
    _ifList.append(new If_Real("prg_fscale", &_fscale));
    _ifList.append(new If_Real("prg_tf", &tf));
    _ifList.append(new If_Real("prg_cs0", &cs0));
    _ifList.append(new If_Bool("prg_controller", &_controller));
    _ifList.append(new If_Real("prg_uinit", &_uinit));
    _ifList.append(new If_RealVec("prg_y", &_y));
}

//--------------------------------------------------------------------------
Prg_Bio::~Prg_Bio()
{
    V_FREE(_y);
}

//--------------------------------------------------------------------------
void Prg_Bio::setup_stages(IVECP ks, VECP ts)
{
    stages_alloc(ks, ts, K(), 1, 0.0, tf);

    // allocate output variables
    _y = v_resize(_y, 2*(K()+1));    
}

//--------------------------------------------------------------------------
void Prg_Bio::setup(int k,
		    Omu_Vector &x, Omu_Vector &u, Omu_Vector &c)
{
    // number of state variables
    x.alloc(2);

    if ( k == 0 ) {
	// fixed initial state
	x.min[0] = x.max[0] = x.initial[0] = p0;
	x.min[1] = x.max[1] = x.initial[1] = 0.0;
    } else {
	// state bounds
	x.min[0] = x.min[1] = 0.0;
    } 

    if ( k < _K ) {
	// number of control variables
	u.alloc(1);
	// control bounds
	u.min[0] = Fsmin;
	u.max[0] = Fsmax;
	// default initial control trajectory 
	u.initial[0] = _uinit;
    }
}

//--------------------------------------------------------------------------
void Prg_Bio::init_simulation(int k,
			      Omu_Vector &x, Omu_Vector &u)
{
    int i;
    double dt = tf/(double) K();

    // set initial states for the first stage;
    // afterwards simulation results of the preceding stage are used
    if ( k == 0 ) {
	for ( i = 0; i < (int) x->dim; i++ )
	    x[i] = x.initial[i];
    }

    // set the initial control for each stage
    if ( k < _K ) {
	if ( _controller )
	    u[0] = Controller((const VECP) x, 0.0+k*dt, dt);
	else
	    u[0] = _uinit;
    }
}

//--------------------------------------------------------------------------
void  Prg_Bio::update(int kk, 
		      const adoublev &x, const adoublev &u, 
		      adoublev &f, adouble &f0, adoublev &c)
{
    double cs, v, s, cp;

    v = v0 +(value(x[0])-p0)/kappa+value(x[1]);
    s = cs0*v0-(value(x[0])-p0)/yps+cdos*value(x[1]);
    cs = s/v;            // substrate concentartion
    if ( cs < 0.0 ) 
	cs = 0.0;
    cp = value(x[0])/v;  // product concentration
    if ( cp < 0.0 ) 
	cp = 0.0;
    // output 
    _y[2*kk] = cs;
    _y[2*kk+1] = cp;

    // cost function
    if ( kk == KK() ) {
	f0 = -((kp+kap/kappa)*x[0]-(kos*cdos+kap)*x[1]-kap*v0 + kap/kappa*p0);
	f0 *= _fscale;
    }
}

//--------------------------------------------------------------------------
void Prg_Bio::continuous(int kk, double t, 
			 const adoublev &x, const adoublev &u,
			 const adoublev &xp, adoublev &F)
{
    adouble v, s, cs, cp, Pi;

    v = v0 +(x[0]-p0)/kappa+x[1];
    s = cs0*v0-(x[0]-p0)/yps+cdos*x[1];
//    cs = fmax(0.0, s/v);    // substrate concentration
    cs = (s/v>0.0) ? s/v : (adouble) 0.0;
//    cp = fmax(x[0]/v, 0.0); // product concentration
    cp = (x[0]/v>0.0) ? x[0]/v : (adouble) 0.0;
    Pi = x0*pimax*exp(-kd*t-cp/kip)*cs/(ks+cs+cs*cs/kis);

    // state equations
    F[0] = Pi - xp[0]; 
    F[1] = u[0]-xp[1];
}

//--------------------------------------------------------------------------
// nonlinear state controller for generation of initial control trajectory
double Prg_Bio::Controller(const VECP x, const double t, const double dt)
{
    double cs_ref, Fs, kr, v, s, cs, cp, Pi;

    cs_ref = sqrt(ks*kis); // substrate reference value
    kr = 0.2-1;            // controller gain

    v = fabs(v0 +(x[0]-p0)/kappa+x[1]);
    s = fabs(cs0*v0-(x[0]-p0)/yps+cdos*x[1]);
    cs = s/v;
    cp = fabs(x[0]/v);
    Pi = x0*pimax*exp(-kd*t-cp/kip)*cs/(ks+cs+cs*cs/kis);
  
    Fs = (kr/dt*(cs-cs_ref)*v*v+Pi*(s/kappa+v/yps))/(cdos*v-s);
    Fs = (Fs>Fsmax) ? Fsmax : ( (Fs<Fsmin) ? Fsmin : Fs );

    return Fs;
}

//==========================================================================

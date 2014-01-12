/*
 * Prg_DIC_SFunction.C -- class definition
 * (derived from double integrator example accessing Docp interface)
 *
 * rf, 05/06/01
 */

#include "Prg_DIC_SFunction.h"

#include <assert.h>
#include <If_Bool.h>

IF_CLASS_DEFINE("DIC_SFunction", Prg_DIC_SFunction, Omu_Program);

//--------------------------------------------------------------------------
Prg_DIC_SFunction::Prg_DIC_SFunction()
{
  set_K(60);
  _nx = 2;
  _nu = 1;
  _mdl_ny = 2;
  _dt = 1.0/K();

 _with_cns = true;

  _SS = NULL;

  // interface elements
  _ifList.append(new If_Bool("prg_with_cns", &_with_cns));
}

//--------------------------------------------------------------------------
Prg_DIC_SFunction::~Prg_DIC_SFunction()
{
  if (_SS) {
    mdlTerminate(_SS);
    Hxi_SimStruct_destroy(_SS);
  }
}

//--------------------------------------------------------------------------
void Prg_DIC_SFunction::setup_stages(IVECP ks, VECP ts)
{
  // (re)create a new SimStruct
  if (_SS) {
    mdlTerminate(_SS);
    Hxi_SimStruct_destroy(_SS);
  }
  _SS = Hxi_SimStruct_create(""); // no path argument for inlined S-function

  // initialize model sizes
  mdlInitializeSizes(_SS);
  if (ssGetErrorStatus(_SS)) {
    cerr << "Error from mdlInitializeSizes: "
	 << ssGetErrorStatus(_SS) << "\n";
    exit(-1);
  }

  // check for number of parameters
  if (ssGetNumSFcnParams(_SS) != ssGetSFcnParamsCount(_SS)) {
    cerr << "S-function parameter count mismatch: expected "
	 << ssGetNumSFcnParams(_SS) << ", provided "
	 << ssGetSFcnParamsCount(_SS) << "!\n";
    exit(-1);
  }

  // check model sizes
  assert(_nx == ssGetNumContStates(_SS));
  assert(ssGetNumInputPorts(_SS) == 1);
  assert(_nu == ssGetInputPortWidth(_SS, 0));
  assert(ssGetNumOutputPorts(_SS) == 1);
  assert(_mdl_ny == ssGetOutputPortWidth(_SS, 0));
  // ... further checks may follow

  // initialize sample times of model
  mdlInitializeSampleTimes(_SS);
  if (ssGetErrorStatus(_SS)) {
    cerr << "Error from mdlInitializeSampleTimes: "
	 << ssGetErrorStatus(_SS) << "\n";
    exit(-1);
  }

  // check sample times of model
  assert(ssGetNumSampleTimes(_SS) == 1);
  assert(value(ssGetSampleTime(_SS, 0)) == CONTINUOUS_SAMPLE_TIME);

  // allocate stages for optimization
  stages_alloc(ks, ts, K(), 1, 0.0, 1.0);
}

//--------------------------------------------------------------------------
void Prg_DIC_SFunction::setup(int k,
			      Omu_Vector &x, Omu_Vector &u, Omu_Vector &c)
{
  x.alloc(_nx);
  if (k < K())
    u.alloc(_nu);

  // initial values
  if (k == 0) {
    // get initial states from model
    mdlInitializeConditions(_SS);
    if (ssGetErrorStatus(_SS)) {
      cerr << "Error from mdlInitializeConditions: "
	   << ssGetErrorStatus(_SS) << "\n";
      ssSetErrorStatus(_SS, NULL);
    }
    real_T *mdl_x = ssGetContStates(_SS);
    x.initial[0] = value(mdl_x[0]);
    x.initial[1] = value(mdl_x[1]);
  }
  if (k < K())
    u.initial[0] = -2.0;

  // initial state constraints
  if (k == 0) {
    x.min[0] = x.max[0] = x.initial[0];
    x.min[1] = x.max[1] = x.initial[1];
  }
  // path constraint
  else if (k < K()) {
    x.max[1] = 0.01;
  }
  // final state constraints
  else {
    x.min[0] = x.max[0] = -1.0;
    x.min[1] = x.max[1] = 0.0;
  }

  // optional additional treatment for path constraint
  if (_with_cns) {
    if (k < K()) {
      c.alloc(1);
      c.max[0] = 0.01;
    }
  }
}

//--------------------------------------------------------------------------
void Prg_DIC_SFunction::update(int kk,
			       const adoublev &x, const adoublev &u,
			       adoublev &f, adouble &f0, adoublev &c)
{
  int i;

  // update constraints and objective for given x and u
  if (kk < KK()) {
    // get pointers to model variables
    real_T *mdl_x = ssGetContStates(_SS);
    real_T **mdl_u = (real_T **)ssGetInputPortRealSignalPtrs(_SS, 0);
    real_T *mdl_y = ssGetOutputPortRealSignal(_SS, 0);

    // initialize model states and inputs
    for (i = 0; i < _nx; ++i)
      mdl_x[i] = x[i];
    for (i = 0; i < _nu; ++i)
      *mdl_u[i] = u[i];

    // obtain model outputs for current states and inputs
    mdlOutputs(_SS, 0);
    if (ssGetErrorStatus(_SS)) {
      cerr << "Error from mdlOutputs: " << ssGetErrorStatus(_SS) << "\n";
      ssSetErrorStatus(_SS, NULL);
    }

    // objective
    f0 = u[0] * u[0] * _dt;

    // constraints
    if (_with_cns) {
      c[0] = mdl_y[1] + 0.5*_dt*mdl_y[0];
    }
  }
  else
    f0 = 0.0;
}

//--------------------------------------------------------------------------
void Prg_DIC_SFunction::continuous(int kk, double t,
				   const adoublev &x, const adoublev &u,  
				   const adoublev &dx, adoublev &F)
{
  int i;

  // get pointers to model variables
  real_T *mdl_x = ssGetContStates(_SS);
  real_T **mdl_u = (real_T **)ssGetInputPortRealSignalPtrs(_SS, 0);
  real_T *mdl_dx = ssGetdX(_SS);

  // set simulation time
  ssSetT(_SS, t);

  // initialize model states and inputs
  for (i = 0; i < _nx; ++i)
    mdl_x[i] = x[i];
  for (i = 0; i < _nu; ++i)
    *mdl_u[i] = u[i];

  // evaluate continuous model equations
  mdlDerivatives(_SS);
  if (ssGetErrorStatus(_SS)) {
    cerr << "Error from mdlDerivatives: " << ssGetErrorStatus(_SS) << "\n";
    ssSetErrorStatus(_SS, NULL);
  }

  // change to residual form
  for (i = 0; i < _nx; ++i)
    F[i] = mdl_dx[i] - dx[i];
}

//==========================================================================

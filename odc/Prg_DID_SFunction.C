/*
 * Prg_DID_SFunction.C -- class definition
 * (derived from double integrator example accessing Docp interface)
 *
 * rf, 05/05/01
 */

#include "Prg_DID_SFunction.h"

#include <assert.h>
#include <If_Bool.h>
#include <If_Real.h>

// include model definition
#include "sfun_did.c"

IF_CLASS_DEFINE("DID_SFunction", Prg_DID_SFunction, Omu_Program);

//--------------------------------------------------------------------------
Prg_DID_SFunction::Prg_DID_SFunction()
{
  set_K(60);
  _nx = 2;
  _nu = 1;
  _mdl_ny = 2;
  _dt = 1.0/K();

  _mx_dt = mxCreateDoubleMatrix(1, 1, mxREAL);

  _with_cns = true;

  _S = NULL;

  // interface elements
  _ifList.append(new If_Bool("prg_with_cns", &_with_cns));
  _ifList.append(new If_Real("prg_dt", &_dt));
}

//--------------------------------------------------------------------------
Prg_DID_SFunction::~Prg_DID_SFunction()
{
  if (_S) {
    mdlTerminate(_S);
    Hxi_SimStruct_destroy(_S);
  }
  mxDestroyArray(_mx_dt);
}

//--------------------------------------------------------------------------
void Prg_DID_SFunction::setup_stages(IVECP ks, VECP ts)
{
  // (re)create a new SimStruct
  if (_S) {
    mdlTerminate(_S);
    Hxi_SimStruct_destroy(_S);
  }
  _S = Hxi_SimStruct_create();

  // initialize model parameters
  mxGetPr(_mx_dt)[0] = _dt;

  // pass model parameters to SimStruct
  ssSetSFcnParamsCount(_S, 1);
  ssSetSFcnParam(_S, 0, _mx_dt);

  // initialize model sizes
  mdlInitializeSizes(_S);

  // check for initialization errors
  if (ssGetErrorStatus(_S)) {
    cerr << "S-function initialization error: "
	 << ssGetErrorStatus(_S) << "\n";
    exit(-1);
  }
  if (ssGetNumSFcnParams(_S) != ssGetSFcnParamsCount(_S)) {
    cerr << "S-function parameter count mismatch: expected "
	 << ssGetNumSFcnParams(_S) << ", provided "
	 << ssGetSFcnParamsCount(_S) << "!\n";
    exit(-1);
  }

  // check model sizes
  assert(_nx == ssGetNumDiscStates(_S));
  assert(ssGetNumInputPorts(_S) == 1);
  assert(_nu == ssGetInputPortWidth(_S, 0));
  assert(ssGetNumOutputPorts(_S) == 1);
  assert(_mdl_ny == ssGetOutputPortWidth(_S, 0));
  // ... further checks may follow

  // initialize sample times of model
  mdlInitializeSampleTimes(_S);

  // check sample times of model
  assert(ssGetNumSampleTimes(_S) == 1);
  assert(value(ssGetSampleTime(_S, 0)) > 0.0);

  // allocate stages for optimization
  stages_alloc(ks, ts, K(), 1, 0.0, K()*value(ssGetSampleTime(_S, 0)));
}

//--------------------------------------------------------------------------
void Prg_DID_SFunction::setup(int k,
			      Omu_Vector &x, Omu_Vector &u, Omu_Vector &c)
{
  x.alloc(_nx);
  if (k < K())
    u.alloc(_nu);

  // initial values
  if (k == 0) {
    // get initial states from model
    mdlInitializeConditions(_S);
    real_T *mdl_x = ssGetDiscStates(_S);
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
void Prg_DID_SFunction::update(int kk,
			       const adoublev &x, const adoublev &u,
			       adoublev &f, adouble &f0, adoublev &c)
{
  int i;
  real_T dt = ssGetSampleTime(_S, 0);

  // update constraints and objective for given x and u
  if (kk < KK()) {
    real_T *mdl_x = ssGetDiscStates(_S);
    real_T **mdl_u = (real_T **)ssGetInputPortRealSignalPtrs(_S, 0);
    real_T *mdl_y = ssGetOutputPortRealSignal(_S, 0);

    // initialize model states and inputs
    for (i = 0; i < _nx; ++i)
      mdl_x[i] = x[i];
    for (i = 0; i < _nu; ++i)
      *mdl_u[i] = u[i];

    // obtain model outputs for current states and inputs
    mdlOutputs(_S, 0);

    // objective
    f0 = u[0] * u[0] * dt;

    // constraints
    if (_with_cns) {
      c[0] = mdl_y[1] + 0.5*dt*mdl_y[0];
    }

    // evaluate model equations
    mdlUpdate(_S, 0);

    // get new states
    for (i = 0; i < _nx; ++i)
      f[i] = mdl_x[i];
  }
  else
    f0 = 0.0;
}


//==========================================================================

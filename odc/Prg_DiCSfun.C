/*
 * Prg_DiCSfun.C -- class definition
 * (derived from Di example accessing Docp interface)
 *
 * rf, 05/06/01
 */

#include "Prg_DiCSfun.h"

#include <assert.h>
#include <If_Bool.h>

// include model definition
#include "sfun_dic.c"

/** Define if local memory should be used for adoubles.
    (required to avoid overflow of max number of live active variables) */
#define PRG_WITH_LOCAL_MEMORY 1

IF_CLASS_DEFINE("DiCSfun", Prg_DiCSfun, Omu_Program);

//--------------------------------------------------------------------------
Prg_DiCSfun::Prg_DiCSfun()
{
  _K = 60;
  _nx = 2;
  _nu = 1;
  _mdl_ny = 2;
  _dt = 1.0/_K;

 _with_cns = true;

  // create SimStruct for communication with S-function
  _S = new SimStruct;

  // interface elements
  _ifList.append(new If_Bool("prg_with_cns", &_with_cns));
}

//--------------------------------------------------------------------------
Prg_DiCSfun::~Prg_DiCSfun()
{
  delete _S;
}

//--------------------------------------------------------------------------
void Prg_DiCSfun::setup_stages(IVECP ks, VECP ts)
{
  // initialize model sizes
  mdlInitializeSizes(_S);

  // check model sizes
  assert(_nx == ssGetNumContStates(_S));
  assert(ssGetNumInputPorts(_S) == 1);
  assert(_nu == ssGetInputPortWidth(_S, 0));
  assert(ssGetNumOutputPorts(_S) == 1);
  assert(_mdl_ny == ssGetOutputPortWidth(_S, 0));
  // ... further checks may follow

  // initialize sample times of model
  mdlInitializeSampleTimes(_S);

  // check sample times of model
  assert(ssGetNumSampleTimes(_S) == 1);
  assert(value(ssGetSampleTime(_S, 0)) == CONTINUOUS_SAMPLE_TIME);

  // allocate stages for optimization
  stages_alloc(ks, ts, _K, 1, 0.0, 1.0);
}

//--------------------------------------------------------------------------
void Prg_DiCSfun::setup(int k,
			Omu_Vector &x, Omu_Vector &u, Omu_Vector &c)
{
  x.alloc(_nx);
  if (k < _K)
    u.alloc(_nu);

  // initial values
  if (k == 0) {
    // get initial states from model
#if defined(PRG_WITH_LOCAL_MEMORY)
    // allocate local memory
    adouble work[_nx];
    _S->set_xc_ext(work);
#endif
    mdlInitializeConditions(_S);
    real_T *mdl_x = ssGetContStates(_S);
    x.initial[0] = value(mdl_x[0]);
    x.initial[1] = value(mdl_x[1]);
#if defined(PRG_WITH_LOCAL_MEMORY)
    // reset as local memory is freed
    _S->set_xc_ext(NULL);
#endif
  }
  if (k < _K)
    u.initial[0] = -2.0;

  // initial state constraints
  if (k == 0) {
    x.min[0] = x.max[0] = x.initial[0];
    x.min[1] = x.max[1] = x.initial[1];
  }
  // path constraint
  else if (k < _K) {
    x.max[1] = 0.01;
  }
  // final state constraints
  else {
    x.min[0] = x.max[0] = -1.0;
    x.min[1] = x.max[1] = 0.0;
  }

  // optional additional treatment for path constraint
  if (_with_cns) {
    if (k < _K) {
      c.alloc(1);
      c.max[0] = 0.01;
    }
  }
}

//--------------------------------------------------------------------------
void Prg_DiCSfun::update(int kk,
			 const adoublev &x, const adoublev &u,
			 adoublev &f, adouble &f0, adoublev &c)
{
#if defined(PRG_WITH_LOCAL_MEMORY)
  // allocate local memory
  adouble work[_nx + _nu + _mdl_ny];
  _S->set_xc_ext(work);
  _S->set_u_ext(work + _nx);
  _S->set_y_ext(work + _nx + _nu);
#endif

  int i;

  // update constraints and objective for given x and u
  if (kk < _KK) {
    // get pointers to model variables
    real_T *mdl_x = ssGetContStates(_S);
    real_T *mdl_u = ssGetInputPortRealSignal(_S, 0);
    real_T *mdl_y = ssGetOutputPortRealSignal(_S, 0);

    // initialize model states and inputs
    for (i = 0; i < _nx; ++i)
      mdl_x[i] = x[i];
    for (i = 0; i < _nu; ++i)
      mdl_u[i] = u[i];

    // obtain model outputs for current states and inputs
    mdlOutputs(_S, 0);

    // objective
    f0 = u[0] * u[0] * _dt;

    // constraints
    if (_with_cns) {
      c[0] = mdl_y[1] + 0.5*_dt*mdl_y[0];
    }
  }
  else
    f0 = 0.0;

#if defined(PRG_WITH_LOCAL_MEMORY)
  // reset as local memory is freed
  _S->set_xc_ext(NULL);
  _S->set_u_ext(NULL);
  _S->set_y_ext(NULL);
#endif
}

//--------------------------------------------------------------------------
void Prg_DiCSfun::continuous(int kk, double t,
			     const adoublev &x, const adoublev &u,  
			     const adoublev &dx, adoublev &F)
{
#if defined(PRG_WITH_LOCAL_MEMORY)
  // allocate local memory
  adouble work[_nx + _nu + _nx];
  _S->set_xc_ext(work);
  _S->set_u_ext(work + _nx);
  _S->set_dxc_ext(work + _nx + _nu);
#endif

  int i;

  // get pointers to model variables
  real_T *mdl_x = ssGetContStates(_S);
  real_T *mdl_u = ssGetInputPortRealSignal(_S, 0);
  real_T *mdl_dx = ssGetdX(_S);

  // set simulation time
  ssSetT(_S, t);

  // initialize model states and inputs
  for (i = 0; i < _nx; ++i)
    mdl_x[i] = x[i];
  for (i = 0; i < _nu; ++i)
    mdl_u[i] = u[i];

  // evaluate continuous model equations
  mdlDerivatives(_S);

  // change to residual form
  for (i = 0; i < _nx; ++i)
    F[i] = mdl_dx[i] - dx[i];

#if defined(PRG_WITH_LOCAL_MEMORY)
  // reset as local memory is freed
  _S->set_xd_ext(NULL);
  _S->set_u_ext(NULL);
  _S->set_y_ext(NULL);
#endif
}

//==========================================================================

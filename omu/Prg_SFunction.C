/*
 * Prg_SFunction.C -- class definition
 *
 */

/*
    Copyright (C) 1997--2005  Ruediger Franke

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Library General Public
    License as published by the Free Software Foundation; 
    version 2 of the License.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Library General Public License for more details.

    You should have received a copy of the GNU Library General Public
    License along with this library (file COPYING.LIB);
    if not, write to the Free Software Foundation, Inc.,
    59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include "Prg_SFunction.h"

#include <stdlib.h>

#include <Hxi_mx_parse.h>

#include <If_RealVec.h>
#include <If_String.h>

// redefine assert to throw an error instead of aborting
#undef assert
#define assert(expr) if (!(expr)) m_error(E_INTERN, "assert(" #expr ")");

#define GET_SET_CB(vartype, name) \
  #name, \
  IF_GET_CB(vartype, Prg_SFunction, name), \
  IF_SET_CB(vartype, Prg_SFunction, set_##name)

// Call an S-function method and check for errors.
// Throw E_FORMAT as errors occuring during initialization done here
// are generally due to wrong configuration data (e.g. bad mdl_args).
#define SMETHOD_CALL(method, S) \
  ssSetErrorStatus(S, NULL); \
  method(S); \
  if (ssGetErrorStatus(S)) { \
    fprintf(stderr, "Error from " #method ": %s\n", \
	    ssGetErrorStatus(S)); \
    m_error(E_FORMAT, ssGetErrorStatus(S)); \
  }

#define SMETHOD_CALL2(method, S, tid) \
  ssSetErrorStatus(S, NULL); \
  method(S, tid); \
  if (ssGetErrorStatus(S)) { \
    fprintf(stderr, "Error from " #method ": %s\n", \
	    ssGetErrorStatus(S)); \
    m_error(E_FORMAT, ssGetErrorStatus(S)); \
  }

//--------------------------------------------------------------------------
Prg_SFunction::Prg_SFunction()
{
  _mdl_name = strdup("SFunction");
  _mdl_path = strdup("");

  // initialize _mdl_args and _mx_args
  _mdl_args = strdup("");
  _mdl_nargs = 0;
  _mx_args = NULL;

  _mdl_needs_setup = true;

  _SS = NULL;

  _mdl_np = 0;
  _mdl_nx = 0;
  _mdl_nu = 0;
  _mdl_ny = 0;

  _mdl_p = v_get(_mdl_np);
  _mdl_x0 = v_get(_mdl_nx);

  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, mdl_p)));
  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, mdl_x0)));
  _ifList.append(new If_String(GET_SET_CB(const char *, mdl_name)));
  _ifList.append(new If_String(GET_SET_CB(const char *, mdl_path)));
  _ifList.append(new If_String(GET_SET_CB(const char *, mdl_args)));
}

//--------------------------------------------------------------------------
Prg_SFunction::~Prg_SFunction()
{
  int i;
  if (_SS) {
    mdlTerminate(_SS);
    Hxi_SimStruct_destroy(_SS);
  }
  for (i = 0; i < _mdl_nargs; i++)
    mxDestroyArray(_mx_args[i]);
  delete [] _mx_args;

  v_free(_mdl_x0);
  v_free(_mdl_p);
  free(_mdl_args);
  free(_mdl_path);
  free(_mdl_name);
}

//--------------------------------------------------------------------------
void Prg_SFunction::set_mdl_name(const char *str)
{
  free(_mdl_name);
  _mdl_name = strdup(str);
  _mdl_needs_setup = true;
}

//--------------------------------------------------------------------------
void Prg_SFunction::set_mdl_path(const char *str)
{
  free(_mdl_path);
  _mdl_path = strdup(str);
  _mdl_needs_setup = true;
}

//--------------------------------------------------------------------------
void Prg_SFunction::set_mdl_args(const char *arg_str)
{
  const char *str, *str1;
  mxArray **args;
  int i, nargs = 0;

  // parse args
  str = arg_str;
  str1 = Hxi::mx_count_columns(str, nargs);
  if (*str1 != '\0') {
    // did not arrive at the end of the string
    m_error(E_FORMAT, "Prg_SFunction::set_mdl_args that "
	    "failed to parse S-function args");
  }
  args = new mxArray* [nargs];
  for (i = 0; i < nargs; i++) {
    args[i] = Hxi::mx_parse_argument(_SS, str);
    str = Hxi::mx_forward_argument(str);
    if (*str == ',')
      str = Hxi::mx_forward_whitespaces(++str);  // skip arg delimiter
  }

  // take over successfully parsed args
  free(_mdl_args);
  for (i = 0; i < _mdl_nargs; i++)
    mxDestroyArray(_mx_args[i]);
  delete [] _mx_args;
  _mx_args = args;
  _mdl_nargs = nargs;
  _mdl_args = strdup(arg_str);

  _mdl_needs_setup = true;
}

//--------------------------------------------------------------------------
void Prg_SFunction::read_mx_args(VECP p)
{
  mxArray *arg;
  int i, j, idx, nel;

  for (idx = 0, j = 0; j < _mdl_nargs; j++) {
    arg = _mx_args[j];
    if (mxIsDouble(arg)) {
      nel = mxGetNumberOfElements(arg);
      for (i = 0; i < nel; i++, idx++)
	p[idx] = mxGetPr(arg)[i];
    }
  }
  assert(idx == _mdl_np); // S-function parameters must not have changed
}

//--------------------------------------------------------------------------
void Prg_SFunction::write_mx_args(VECP p)
{
  mxArray *arg;
  int i, j, idx, nel;

  for (idx = 0, j = 0; j < _mdl_nargs; j++) {
    arg = _mx_args[j];
    if (mxIsDouble(arg)) {
      nel = mxGetNumberOfElements(arg);
      for (i = 0; i < nel; i++, idx++)
	mxGetPr(arg)[i] = p[idx];
    }
  }
  assert(idx == _mdl_np); // S-function parameters must not have changed
}

//--------------------------------------------------------------------------
void Prg_SFunction::setup_model()
{
  int i;

  // setup S-function
  if (_SS) {
    mdlTerminate(_SS);
    Hxi_SimStruct_destroy(_SS);
  }

  _SS = Hxi_SimStruct_create(_mdl_path[0] != '\0'? _mdl_path: _mdl_name);
  if (ssGetErrorStatus(_SS)) {
    fprintf(stderr, "Error creating SimStruct: %s\n", ssGetErrorStatus(_SS));
    m_error(E_FORMAT, ssGetErrorStatus(_SS));
  }

  // initialize model name
  ssSetModelName(_SS, _mdl_name);

  // initialize S-function parameters
  ssSetSFcnParamsCount(_SS, _mdl_nargs);
  for (i = 0; i < _mdl_nargs; i++)
    ssSetSFcnParam(_SS, i, _mx_args[i]);

  // initialize solver
  // (note: variable step size is indicated as the model must allow
  //  simulation time stepping back)
  ssSetVariableStepSolver(_SS, 1);
  // (note: preselect major time steps requiring complete model evaluation --
  //  minor time steps would require support for events and zero crossings)
  ssSetMinorTimeStep(_SS, 0);

  // initialize model
  SMETHOD_CALL(mdlInitializeSizes, _SS);
  if (ssGetNumSFcnParams(_SS) != ssGetSFcnParamsCount(_SS)) {
    fprintf(stderr, "Parameter count mismatch: expected: %d, provided: %d\n",
	    ssGetNumSFcnParams(_SS), ssGetSFcnParamsCount(_SS));
    m_error(E_FORMAT, "Prg_SFunction::setup_model: parameter count mismatch");
  }

  // obtain model sizes
  _mdl_nx = ssGetNumContStates(_SS);
  // count inputs of ports as long as they are contiguous in memory
  // (this limitation is because later on we will only access port 0)
  _mdl_nu = 0;
  for (i = 0; i < ssGetNumInputPorts(_SS); i++) {
    if (*ssGetInputPortRealSignalPtrs(_SS, i)
	== *ssGetInputPortRealSignalPtrs(_SS, 0) + _mdl_nu)
      _mdl_nu += ssGetInputPortWidth(_SS, i);
    else {
      m_warning(WARN_UNKNOWN,
		"Prg_SFunction::setup_model: ignoring non-contiguous inputs");
      break;
    }
  }
  // count outputs of ports as long as they are contiguous in memory
  // (this limitation is because later on we will only access port 0)
  _mdl_ny = 0;
  for (i = 0; i < ssGetNumOutputPorts(_SS); i++) {
    if (ssGetOutputPortRealSignal(_SS, i)
	== ssGetOutputPortRealSignal(_SS, 0) + _mdl_ny)
      _mdl_ny += ssGetOutputPortWidth(_SS, i);
    else {
      m_warning(WARN_UNKNOWN,
		"Prg_SFunction::setup_model: ignoring non-contiguous outputs");
      break;
    }
  }

  // set simulation time
  ssSetT(_SS, _t0);

  // initialize and check sample times of model
  mdlInitializeSampleTimes(_SS);
  assert(ssGetNumSampleTimes(_SS) > 0);
  assert(value(ssGetSampleTime(_SS, 0)) == CONTINUOUS_SAMPLE_TIME);
  if (ssGetNumSampleTimes(_SS) > 1)
    m_warning(WARN_UNKNOWN,
	      "Prg_SFunction::setup_model: ignoring multiple sample times");

  // start using S-function
  if (ssGetmdlStart(_SS) != NULL)
    mdlStart(_SS);

  // get parameters
  // determine number of parameters
  mxArray *arg;
  _mdl_np = 0;
  for (int j = 0; j < _mdl_nargs; j++) {
    arg = _mx_args[j];
    // only consider parameters in double format for accessing via mxGetPr()
    if (mxIsDouble(arg))
      _mdl_np += mxGetNumberOfElements(arg);
  }
  v_resize(_mdl_p, _mdl_np);
  read_mx_args(_mdl_p);

  // get initial states
  // (states shall be initialized with mdlInitializeConditions;
  //  mdlOutputs is called according to the docu afterwards)
  if (ssGetmdlInitializeConditions(_SS) != NULL) {
    SMETHOD_CALL(mdlInitializeConditions, _SS);
  }
  SMETHOD_CALL2(mdlOutputs, _SS, 0);
  real_T *mdl_x = ssGetContStates(_SS);
  v_resize(_mdl_x0, _mdl_nx);
  for (i = 0; i < _mdl_nx; i++)
    _mdl_x0[i] = mdl_x[i];

  _mdl_needs_setup = false;
}


//==========================================================================

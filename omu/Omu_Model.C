/*
 * Omu_Model.C -- class definition
 *
 */

/*
    Copyright (C) 1997--2014  Ruediger Franke

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

#include "Omu_Model.h"

#include <stdlib.h>

#include <Hxi_mx_parse.h>

#include <If_RealVec.h>
#include <If_String.h>

// redefine assert to throw an error instead of aborting
#undef assert
#define assert(expr) if (!(expr)) m_error(E_INTERN, "assert(" #expr ")");

#define GET_SET_CB(vartype, name) \
  #name, \
  IF_GET_CB(vartype, Omu_Model, name), \
  IF_SET_CB(vartype, Omu_Model, set_##name)

// Call an S-function method and check for errors.
// Throw E_FORMAT as errors occuring during initialization done here
// are generally due to wrong configuration data (e.g. bad mdl_args).
#define SMETHOD_CALL(method, S) { \
  ssSetErrorStatus(S, NULL); \
  method(S); \
  if (ssGetErrorStatus(S)) { \
    fprintf(stderr, "Error from " #method ": %s\n", \
	    ssGetErrorStatus(S)); \
    m_error(E_FORMAT, ssGetErrorStatus(S)); \
  } \
}
#define SMETHOD_CALL2(method, S, tid) { \
  ssSetErrorStatus(S, NULL); \
  method(S, tid); \
  if (ssGetErrorStatus(S)) { \
    fprintf(stderr, "Error from " #method ": %s\n", \
	    ssGetErrorStatus(S)); \
    m_error(E_FORMAT, ssGetErrorStatus(S)); \
  } \
}

//--------------------------------------------------------------------------
Omu_Model::Omu_Model()
{
  _mdl_name = strdup("SFunction");
  _mdl_path = strdup("");
  _mdl_is_fmu = false;

  // initialize _mdl_args and _mx_args
  _mdl_args = strdup("");
  _mdl_nargs = 0;
  _mx_args = NULL;

  _mdl_needs_setup = true;

  _SS = NULL;

  _mdl_np = 0;
  _mdl_nd = 0;
  _mdl_nx = 0;
  _mdl_nu = 0;
  _mdl_ny = 0;

  _mdl_p = v_get(_mdl_np);
  _mdl_x_start = v_get(_mdl_nx);
  _mdl_u_start = v_get(_mdl_nu);
  _mdl_x_nominal = v_get(_mdl_nx);
  _mdl_u_nominal = v_get(_mdl_nu);
  _mdl_y_nominal = v_get(_mdl_ny);

  _ifList_model.append(new If_RealVec(GET_SET_CB(const VECP, mdl_p)));
  _ifList_model.append(new If_RealVec(GET_SET_CB(const VECP, mdl_x_start)));
  _ifList_model.append(new If_RealVec(GET_SET_CB(const VECP, mdl_u_start)));
  _ifList_model.append(new If_RealVec(GET_SET_CB(const VECP, mdl_x_nominal)));
  _ifList_model.append(new If_RealVec(GET_SET_CB(const VECP, mdl_u_nominal)));
  _ifList_model.append(new If_RealVec(GET_SET_CB(const VECP, mdl_y_nominal)));
  _ifList_model.append(new If_String(GET_SET_CB(const char *, mdl_name)));
  _ifList_model.append(new If_String(GET_SET_CB(const char *, mdl_path)));
  _ifList_model.append(new If_String(GET_SET_CB(const char *, mdl_args)));
}

//--------------------------------------------------------------------------
Omu_Model::~Omu_Model()
{
  int i;
  if (_SS) {
    mdlTerminate(_SS);
    Hxi_SimStruct_destroy(_SS);
  }
  for (i = 0; i < _mdl_nargs; i++)
    mxDestroyArray(_mx_args[i]);
  delete [] _mx_args;

  v_free(_mdl_y_nominal);
  v_free(_mdl_u_nominal);
  v_free(_mdl_x_nominal);
  v_free(_mdl_u_start);
  v_free(_mdl_x_start);
  v_free(_mdl_p);
  free(_mdl_args);
  free(_mdl_path);
  free(_mdl_name);
}

//--------------------------------------------------------------------------
void Omu_Model::set_mdl_name(const char *str)
{
  free(_mdl_name);
  _mdl_name = strdup(str);
  _mdl_needs_setup = true;
}

//--------------------------------------------------------------------------
void Omu_Model::set_mdl_path(const char *str)
{
  free(_mdl_path);
  _mdl_path = strdup(str);
  _mdl_needs_setup = true;
  int len = strlen(_mdl_path);
  _mdl_is_fmu = (len > 4 && strcmp(_mdl_path + len - 4, ".fmu") == 0);
}

//--------------------------------------------------------------------------
void Omu_Model::set_mdl_args(const char *arg_str)
{
  const char *str, *str1;
  mxArray **args;
  int i, nargs = 0;

  // parse args
  str = arg_str;
  str1 = Hxi::mx_count_columns(str, nargs);
  if (*str1 != '\0') {
    // did not arrive at the end of the string
    m_error(E_FORMAT, "Omu_Model::set_mdl_args that "
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
  for (i = 0; i < _mdl_nargs; i++)
    mxDestroyArray(_mx_args[i]);
  delete [] _mx_args;
  _mx_args = args;
  _mdl_nargs = nargs;
  char *mdl_args = strdup(arg_str);
  free(_mdl_args);
  _mdl_args = mdl_args; // Note: intermediate mdl_args in case arg_str == _mdl_args

  _mdl_needs_setup = true;
}

//--------------------------------------------------------------------------
void Omu_Model::read_mx_args(VECP p)
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
void Omu_Model::write_mx_args(VECP p)
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
bool Omu_Model::setContinuousTask(bool val)
{
  bool hasContinuousSampleTime = false;
  for (int i = 0; i < ssGetNumSampleTimes(_SS); i++) {
    if (ssGetSampleTime(_SS, i) == CONTINUOUS_SAMPLE_TIME) {
      ssGetSampleHitPtr(_SS)[ssGetSampleTimeTaskID(_SS, i)] = val;
      hasContinuousSampleTime = val;
    }
  }
  return hasContinuousSampleTime;
}

//--------------------------------------------------------------------------
bool Omu_Model::setSampleHit(bool val)
{
  bool hasDiscreteSampleTime = false;
  for (int i = 0; i < ssGetNumSampleTimes(_SS); i++) {
    if (ssGetSampleTime(_SS, i) != CONTINUOUS_SAMPLE_TIME) {
      ssGetSampleHitPtr(_SS)[ssGetSampleTimeTaskID(_SS, i)] = val;
      hasDiscreteSampleTime = val;
    }
  }
  return hasDiscreteSampleTime;
}

//--------------------------------------------------------------------------
void Omu_Model::setup_model(double t0)
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

  // initialize model name (FMU name is initialized by model from path below)
  if (!_mdl_is_fmu)
    ssSetModelName(_SS, _mdl_name);

  // (re-)parse model arguments to
  // adapt mxArray to SimStruct of loaded S-function (MEX or Hxi)
  set_mdl_args(_mdl_args);

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
  // exploit model description available for FMU
  if (_mdl_is_fmu) {
    // obtain model name
    free(_mdl_name);
    _mdl_name = strdup(ssGetModelName(_SS));

    // obtain parameters from model description if no mdl_args given
    if (_mdl_nargs == 0 && ssGetNumSFcnParams(_SS) > 0) {
      if (Tcl_VarEval(theInterp, "::set {::fmu::", _mdl_name,
                      "::parameterStartValues}", NULL) != TCL_OK) {
        m_error(E_INTERN, "can't get parameter values");
      }
      set_mdl_args(Tcl_GetStringResult(theInterp));
      Tcl_ResetResult(theInterp);
      ssSetSFcnParamsCount(_SS, _mdl_nargs);
      for (i = 0; i < _mdl_nargs; i++)
        ssSetSFcnParam(_SS, i, _mx_args[i]);
    }
  }
  // check numbers of expected and given model parameters
  if (ssGetNumSFcnParams(_SS) != ssGetSFcnParamsCount(_SS)) {
    fprintf(stderr, "Parameter count mismatch: expected: %d, provided: %d\n",
	    ssGetNumSFcnParams(_SS), ssGetSFcnParamsCount(_SS));
    m_error(E_FORMAT, "Omu_Model::setup_model: parameter count mismatch");
  }

  // obtain model sizes
  _mdl_nd = ssGetNumDiscStates(_SS);
  _mdl_nx = _mdl_nd + ssGetNumContStates(_SS);
  // count inputs of ports as long as they are contiguous in memory
  // (this limitation is because later on we will only access port 0)
  _mdl_nu = 0;
  for (i = 0; i < ssGetNumInputPorts(_SS); i++) {
    if (ssGetInputPortWidth(_SS, i) < 1)
      continue;
    if (*ssGetInputPortRealSignalPtrs(_SS, i)
	== *ssGetInputPortRealSignalPtrs(_SS, 0) + _mdl_nu)
      _mdl_nu += ssGetInputPortWidth(_SS, i);
    else {
      m_warning(WARN_UNKNOWN,
		"Omu_Model::setup_model: ignoring non-contiguous inputs");
      break;
    }
  }
  // count outputs of ports as long as they are contiguous in memory
  // (this limitation is because later on we will only access port 0)
  _mdl_ny = 0;
  for (i = 0; i < ssGetNumOutputPorts(_SS); i++) {
    if (ssGetOutputPortWidth(_SS, i) < 1)
      continue;
    if (ssGetOutputPortRealSignal(_SS, i)
	== ssGetOutputPortRealSignal(_SS, 0) + _mdl_ny)
      _mdl_ny += ssGetOutputPortWidth(_SS, i);
    else {
      m_warning(WARN_UNKNOWN,
		"Omu_Model::setup_model: ignoring non-contiguous outputs");
      break;
    }
  }

  // set simulation time
  ssSetT(_SS, t0);
  _t0_setup_model = t0;

  // initialize and check sample times of model
  mdlInitializeSampleTimes(_SS);
  if (ssGetNumSampleTimes(_SS) < 1)
    m_warning(WARN_UNKNOWN,
	      "Omu_Model::setup_model: no sample times initialized");

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

  v_resize(_mdl_x_start, _mdl_nx);
  if (_mdl_is_fmu) {
    // get start values for states from model description
    if (Tcl_VarEval(theInterp, "mdl_x_start ${::fmu::", _mdl_name,
		    "::stateStartValues}", NULL) != TCL_OK)
      m_error(E_INTERN, "can't obtain start values for states of FMU");
  }
  else {
    // get initial states that shall be initialized with mdlInitializeConditions
    // Note: initialization might also be postponed by the model to the next
    //       mdlOutputs/mdlUpdate call, e.g. if initial states depend on inputs,
    //       however, don't call these methods here as inputs are not known.
    if (ssGetmdlInitializeConditions(_SS) != NULL) {
      SMETHOD_CALL(mdlInitializeConditions, _SS);
    }
    real_T *mdl_xd = ssGetDiscStates(_SS);
    real_T *mdl_xc = ssGetContStates(_SS);
    for (i = 0; i < _mdl_nd; i++)
      _mdl_x_start[i] = mdl_xd[i];
    for (i = _mdl_nd; i < _mdl_nx; i++)
      _mdl_x_start[i] = mdl_xc[i - _mdl_nd];
  }

  v_resize(_mdl_u_start, _mdl_nu);
  if (_mdl_is_fmu) {
    // get start values for inputs from model description
    if (Tcl_VarEval(theInterp, "mdl_u_start ${::fmu::", _mdl_name,
		    "::inputStartValues}", NULL) != TCL_OK)
      m_error(E_INTERN, "can't obtain start values for inputs of FMU");
  }
  else {
    v_set(_mdl_u_start, 0.0);
  }

  // setup nominal values
  v_resize(_mdl_x_nominal, _mdl_nx);
  v_resize(_mdl_u_nominal, _mdl_nu);
  v_resize(_mdl_y_nominal, _mdl_ny);
  if (_mdl_is_fmu) {
    // take over default values from model description
    if(Tcl_VarEval(theInterp, "mdl_u_nominal ${::fmu::", _mdl_name,
                   "::inputNominalValues}", NULL) != TCL_OK ||
       Tcl_VarEval(theInterp, "mdl_x_nominal ${::fmu::", _mdl_name,
                   "::stateNominalValues}", NULL) != TCL_OK ||
       Tcl_VarEval(theInterp, "mdl_y_nominal ${::fmu::", _mdl_name,
                   "::outputNominalValues}", NULL) != TCL_OK)
      m_error(E_INTERN, "can't obtain nominal values");
  }
  else {
    // ordinary S-function
    v_set(_mdl_u_nominal, 1.0);
    v_set(_mdl_x_nominal, 1.0);
    v_set(_mdl_y_nominal, 1.0);
  }

  _mdl_needs_setup = false;
}


//==========================================================================

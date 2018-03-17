/*
 * Omu_Model.C -- class definition
 *
 */

/*
    Copyright (C) 1997--2017  Ruediger Franke

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

#include <Hqp_omp.h>
#include <Hxi_mx_parse.h>

#include <If_Bool.h>
#include <If_Int.h>
#include <If_IntVec.h>
#include <If_RealVec.h>
#include <If_String.h>

// redefine assert to throw an error instead of aborting
#undef assert
#define assert(expr) if (!(expr)) m_error(E_INTERN, "assert(" #expr ")");

#define GET_CB(vartype, name) \
  #name, \
  IF_GET_CB(vartype, Omu_Model, name)

#define GET_SET_CB(vartype, name) \
  #name, \
  IF_GET_CB(vartype, Omu_Model, name), \
  IF_SET_CB(vartype, Omu_Model, set_##name)

//--------------------------------------------------------------------------
Omu_Model::Omu_Model(int ncpu)
{
  _mdl_name = strdup("SFunction");
  _mdl_path = strdup("");
  _mdl_is_fmu = false;

  // initialize _mdl_args and _mx_args
  _mdl_args = strdup("");
  _mdl_nargs = 0;
  _mx_args = NULL;

  _mdl_ncpu = ncpu;
  _mdl_logging = If_LogError;
  _mdl_jac = true;
  _mdl_needs_setup = true;
  _mdl_needs_init = iv_get(_mdl_ncpu);
  iv_set(_mdl_needs_init, 0);

  _SS = new SimStruct* [_mdl_ncpu];
  for (int tn = 0; tn < _mdl_ncpu; tn++)
    _SS[tn] = NULL;

  _mdl_np_total = 0;
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

  _mdl_jac_y_active = IVNULL;
  _mdl_jac_u_active = IVNULL;
  _mdl_jac_jc = IVNULL;
  _mdl_jac_ir = IVNULL;
  _mdl_jac_deps = IVNULL;

  _ifList_model.append(new If_Int(GET_CB(int, mdl_ncpu)));
  _ifList_model.append(new If_Int(GET_SET_CB(int, mdl_logging)));
  _ifList_model.append(new If_Bool(GET_SET_CB(bool, mdl_jac)));
  _ifList_model.append(new If_RealVec(GET_SET_CB(const VECP, mdl_p)));
  _ifList_model.append(new If_RealVec(GET_SET_CB(const VECP, mdl_x_start)));
  _ifList_model.append(new If_RealVec(GET_SET_CB(const VECP, mdl_u_start)));
  _ifList_model.append(new If_RealVec(GET_SET_CB(const VECP, mdl_x_nominal)));
  _ifList_model.append(new If_RealVec(GET_SET_CB(const VECP, mdl_u_nominal)));
  _ifList_model.append(new If_RealVec(GET_SET_CB(const VECP, mdl_y_nominal)));
  _ifList_model.append(new If_String(GET_SET_CB(const char *, mdl_name)));
  _ifList_model.append(new If_String(GET_SET_CB(const char *, mdl_path)));
  _ifList_model.append(new If_String(GET_SET_CB(const char *, mdl_args)));
  _ifList_model.append(new If_IntVec(GET_SET_CB(const IVECP, mdl_jac_deps)));
}

//--------------------------------------------------------------------------
Omu_Model::~Omu_Model()
{
  int i, tn;

  for (tn = 0; tn < _mdl_ncpu; tn++) {
    if (_SS[tn]) {
      SMETHOD_CALL_HOLD(mdlTerminate, _SS[tn]);
      Hxi_SimStruct_destroy(_SS[tn]);
    }
  }
  delete [] _SS;

  for (i = 0; i < _mdl_nargs; i++)
    mxDestroyArray(_mx_args[i]);
  delete [] _mx_args;

  iv_free(_mdl_needs_init);
  iv_free(_mdl_jac_deps);
  iv_free(_mdl_jac_ir);
  iv_free(_mdl_jac_jc);
  iv_free(_mdl_jac_u_active);
  iv_free(_mdl_jac_y_active);
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
  SimStruct *S = _SS[0];

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
    args[i] = Hxi::mx_parse_argument(S, str);
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
bool Omu_Model::setContinuousTask(SimStruct *S, bool val)
{
  bool hasContinuousSampleTime = false;
  for (int i = 0; i < ssGetNumSampleTimes(S); i++) {
    if (ssGetSampleTime(S, i) == CONTINUOUS_SAMPLE_TIME) {
      ssGetSampleHitPtr(S)[ssGetSampleTimeTaskID(S, i)] = val;
      hasContinuousSampleTime = val;
    }
  }
  return hasContinuousSampleTime;
}

//--------------------------------------------------------------------------
bool Omu_Model::setSampleHit(SimStruct *S, bool val)
{
  bool hasDiscreteSampleTime = false;
  for (int i = 0; i < ssGetNumSampleTimes(S); i++) {
    if (ssGetSampleTime(S, i) != CONTINUOUS_SAMPLE_TIME) {
      ssGetSampleHitPtr(S)[ssGetSampleTimeTaskID(S, i)] = val;
      hasDiscreteSampleTime = val;
    }
  }
  return hasDiscreteSampleTime;
}

//--------------------------------------------------------------------------
void Omu_Model::setSampleTime(SimStruct *S, double val)
{
  for (int i = 0; i < ssGetNumSampleTimes(S); i++) {
    if (ssGetSampleTime(S, i) != CONTINUOUS_SAMPLE_TIME) {
      ssSetSampleTime(S, i, val);
    }
  }
}

//--------------------------------------------------------------------------
void Omu_Model::setup_model(double t0)
{
  int i, tn;

  for (tn = 0; tn < _mdl_ncpu; tn++) {
    // setup S-function
    if (_SS[tn]) {
      SMETHOD_CALL(mdlTerminate, _SS[tn]);
      Hxi_SimStruct_destroy(_SS[tn]);
    }

    _SS[tn] = Hxi_SimStruct_create(_mdl_path[0] != '\0'? _mdl_path: _mdl_name);
    if (ssGetErrorStatus(_SS[tn])) {
      fprintf(stderr, "Error creating SimStruct: %s\n", ssGetErrorStatus(_SS[tn]));
      m_error(E_FORMAT, ssGetErrorStatus(_SS[tn]));
    }

    // initialize model name (FMU name is initialized by model from path below)
    if (!_mdl_is_fmu)
      ssSetModelName(_SS[tn], _mdl_name);

    // (re-)parse model arguments to
    // adapt mxArray to SimStruct of loaded S-function (MEX or Hxi)
    if (tn == 0)
      set_mdl_args(_mdl_args);

    // initialize S-function parameters
    ssSetSFcnParamsCount(_SS[tn], _mdl_nargs);
    for (i = 0; i < _mdl_nargs; i++)
      ssSetSFcnParam(_SS[tn], i, _mx_args[i]);

    // initialize solver
    // (note: variable step size is indicated as the model must allow
    //  simulation time stepping back)
    ssSetVariableStepSolver(_SS[tn], 1);
    // (note: preselect major time steps requiring complete model evaluation --
    //  minor time steps would require support for events and zero crossings)
    ssSetMinorTimeStep(_SS[tn], 0);

    // initialize model
    if (_mdl_is_fmu)
      ssSetOptions(_SS[tn], tn);
    SMETHOD_CALL(mdlInitializeSizes, _SS[tn]);
    // exploit model description available for FMU
    if (_mdl_is_fmu) {
      // obtain model name
      if (tn == 0) {
        free(_mdl_name);
        _mdl_name = strdup(ssGetModelName(_SS[tn]));
      }

      // obtain parameters from model description if no mdl_args given
      if (_mdl_nargs == 0 && ssGetNumSFcnParams(_SS[tn]) > 0) {
        if (tn == 0) {
          if (Tcl_VarEval(theInterp, "::set {::fmu::", _mdl_name,
                          "::parameterStartValues}", NULL) != TCL_OK) {
            m_error(E_INTERN, "can't get parameter values");
          }
          set_mdl_args(Tcl_GetStringResult(theInterp));
          Tcl_ResetResult(theInterp);
        }
        ssSetSFcnParamsCount(_SS[tn], _mdl_nargs);
        for (i = 0; i < _mdl_nargs; i++)
          ssSetSFcnParam(_SS[tn], i, _mx_args[i]);
      }
    }
    // check numbers of expected and given model parameters
    if (ssGetNumSFcnParams(_SS[tn]) != ssGetSFcnParamsCount(_SS[tn])) {
      fprintf(stderr, "Parameter count mismatch: expected: %d, provided: %d\n",
              ssGetNumSFcnParams(_SS[tn]), ssGetSFcnParamsCount(_SS[tn]));
      m_error(E_FORMAT, "Omu_Model::setup_model: parameter count mismatch");
    }

    // obtain model sizes
    if (tn == 0) {
      _mdl_nd = ssGetNumDiscStates(_SS[tn]);
      _mdl_nx = _mdl_nd + ssGetNumContStates(_SS[tn]);
      // count inputs of ports as long as they are contiguous in memory
      // (this limitation is because later on we will only access port 0)
      _mdl_nu = 0;
      for (i = 0; i < ssGetNumInputPorts(_SS[tn]); i++) {
        if (ssGetInputPortWidth(_SS[tn], i) < 1)
          continue;
        if (*ssGetInputPortRealSignalPtrs(_SS[tn], i)
            == *ssGetInputPortRealSignalPtrs(_SS[tn], 0) + _mdl_nu)
          _mdl_nu += ssGetInputPortWidth(_SS[tn], i);
        else {
          m_warning(WARN_UNKNOWN,
                    "Omu_Model::setup_model: ignoring non-contiguous inputs");
          break;
        }
      }
      // count outputs of ports as long as they are contiguous in memory
      // (this limitation is because later on we will only access port 0)
      _mdl_ny = 0;
      for (i = 0; i < ssGetNumOutputPorts(_SS[tn]); i++) {
        if (ssGetOutputPortWidth(_SS[tn], i) < 1)
          continue;
        if (ssGetOutputPortRealSignal(_SS[tn], i)
            == ssGetOutputPortRealSignal(_SS[tn], 0) + _mdl_ny)
          _mdl_ny += ssGetOutputPortWidth(_SS[tn], i);
        else {
          m_warning(WARN_UNKNOWN,
                    "Omu_Model::setup_model: ignoring non-contiguous outputs");
          break;
        }
      }

      // setup active variables for Jacobian
      _mdl_jac_y_active = iv_resize(_mdl_jac_y_active, _mdl_ny);
      _mdl_jac_u_active = iv_resize(_mdl_jac_u_active, _mdl_nu);
      iv_zero(_mdl_jac_y_active);
      iv_zero(_mdl_jac_u_active);
    }

    // set simulation time
    ssSetT(_SS[tn], t0);
    _t0_setup_model = t0;

    // initialize and check sample times of model
    SMETHOD_CALL(mdlInitializeSampleTimes, _SS[tn]);
    if (ssGetNumSampleTimes(_SS[tn]) < 1)
      m_warning(WARN_UNKNOWN,
                "Omu_Model::setup_model: no sample times initialized");

    // start using S-function
    if (ssGetmdlStart(_SS[tn]) != NULL)
      mdlStart(_SS[tn]);
  }

  // get parameters
  // determine number of parameters
  mxArray *arg;
  _mdl_np = 0;
  for (i = 0; i < _mdl_nargs; i++) {
    arg = _mx_args[i];
    // only consider parameters in double format for accessing via mxGetPr()
    if (mxIsDouble(arg)) {
      _mdl_np_total += mxGetNumberOfElements(arg);
      _mdl_np += mxGetNumberOfElements(arg);
    }
    else
      _mdl_np_total += 1; // string parameter
  }
  v_resize(_mdl_p, _mdl_np);
  read_mx_args(_mdl_p);

  v_resize(_mdl_x_start, _mdl_nx);
  if (_mdl_is_fmu) {
    // get start values for states from model description
    if (Tcl_VarEval(theInterp, "mdl_x_start ${::fmu::", _mdl_name,
		    "::stateStartValues}", NULL) != TCL_OK)
      m_error(E_INTERN, "can't obtain start values for states of FMU");
    // replace undefined start values for states with zero because
    // they might stay undefined if they are initialized from parameters
    for (i = 0; i < _mdl_nx; i++)
      if (!is_finite(_mdl_x_start[i]))
	_mdl_x_start[i] = 0.0;
  }
  else {
    // get initial states that shall be initialized with mdlInitializeConditions
    // Note: initialization might also be postponed by the model to the next
    //       mdlOutputs/mdlUpdate call, e.g. if initial states depend on inputs,
    //       however, don't call these methods here as inputs are not known.
    if (ssGetmdlInitializeConditions(_SS[0]) != NULL) {
      SMETHOD_CALL(mdlInitializeConditions, _SS[0]);
    }
    real_T *mdl_xd = ssGetDiscStates(_SS[0]);
    real_T *mdl_xc = ssGetContStates(_SS[0]);
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
    if (Tcl_VarEval(theInterp, "mdl_u_nominal ${::fmu::", _mdl_name,
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
  iv_set(_mdl_needs_init, 1);
}

//-------------------------------------------------------------------------
void Omu_Model::setup_jac()
{
  int offs, idx;
  int ir_offset = 0;
  SimStruct *S = _SS[0];
  int m = _mdl_nx + _mdl_ny; // number of rows
  int n = _mdl_nx + _mdl_nu; // number of cols

  // allocate storage for Jacobian structure
  _mdl_jac_jc = iv_resize(_mdl_jac_jc, ssGetJacobianNzMax(S));
  _mdl_jac_ir = iv_resize(_mdl_jac_ir, m + 1);

  // setup sparse structure of Jacobian for active variables
  _mdl_jac_ir[ir_offset] = 0;
  offs = _mdl_np_total;
  for (idx = 0; idx < _mdl_nd; idx++, ir_offset++) {
    setup_jac_row("discreteState", offs + idx, _mdl_jac_ir[ir_offset], ir_offset);
  }
  offs = _mdl_np_total + _mdl_nx + _mdl_nd;
  for (; idx < _mdl_nx; idx++, ir_offset++) {
    setup_jac_row("derivative", offs + idx, _mdl_jac_ir[ir_offset], ir_offset);
  }
  offs = _mdl_np_total + _mdl_nx + _mdl_nx + _mdl_nu;
  for (idx = 0; idx < _mdl_ny; idx++, ir_offset++) {
    if (_mdl_jac_y_active[idx])
      setup_jac_row("output", offs + idx, _mdl_jac_ir[ir_offset], ir_offset);
    else
      // empty row
      _mdl_jac_ir[ir_offset + 1] = _mdl_jac_ir[ir_offset];
  }

  //
  // setup sparse structure for FMU using Jacobian format of S-function
  //
  int_T *ir = ssGetJacobianIr(S);
  int_T *jc = ssGetJacobianJc(S);
  int i, j, jdx;

  // count row elements per column
  memset(jc, 0, (n + 1)*sizeof(int_T));
  for (i = 0; i < m; i++) {
    for (jdx = _mdl_jac_ir[i]; jdx < _mdl_jac_ir[i + 1]; jdx++) {
      j = _mdl_jac_jc[jdx] - _mdl_np_total;
      if (j >= 2*_mdl_nx)
        j -= _mdl_nx; // input
      else if (j >= _mdl_nd)
        j -= _mdl_nd; // continuous state or previous
      jc[j + 1] ++;
    }
  }
  // accumulate row start indices
  for (j = 0; j < n; j++) {
    jc[j+1] += jc[j];
  }
  // fill in row indices
  for (j = n; j > 0; j--) {
    jc[j] -= jc[j] - jc[j-1];
  }
  for (i = 0; i < m; i++) {
    for (jdx = _mdl_jac_ir[i]; jdx < _mdl_jac_ir[i + 1]; jdx++) {
      j = _mdl_jac_jc[jdx] - _mdl_np_total;
      if (j >= 2*_mdl_nx)
        j -= _mdl_nx; // input
      else if (j >= _mdl_nd)
        j -= _mdl_nd; // continuous state or previous
      ir[jc[j + 1]] = i;
      jc[j + 1] ++;
    }
  }

  // copy Jacobian structure to all instances
  for (int tn = 1; tn < _mdl_ncpu; tn++) {
    memcpy(ssGetJacobianIr(_SS[tn]), ir, ssGetJacobianNzMax(S) * sizeof(int_T));
    memcpy(ssGetJacobianJc(_SS[tn]), jc, (n + 1) * sizeof(int_T));
  }
}

//-------------------------------------------------------------------------
void Omu_Model::setup_jac_row(const char *category, int idx,
                              int jc_offset, int ir_offset)
{
  // obtain column indices from model description
  Tcl_Obj *idxObj = Tcl_NewIntObj(idx);
  Tcl_IncrRefCount(idxObj);
  if (Tcl_VarEval(If_Interp(), "mdl_jac_deps ${::fmu::", _mdl_name, "::",
                  category, "Dependencies(", Tcl_GetString(idxObj), ")}",
                  NULL) != TCL_OK)
    m_error(E_INTERN, "can't get model structure");
  Tcl_DecrRefCount(idxObj);

  // take over dependencies from previous/derivatives and active inputs
  int offs = _mdl_np_total + 2*_mdl_nx;
  for (int i = 0; i < _mdl_jac_deps->dim; i++) {
    int j = _mdl_jac_deps[i];
    if (offs <= j && j < offs + _mdl_nu && _mdl_jac_u_active[j - offs] == 0)
      continue;
    _mdl_jac_jc[jc_offset++] = j;
  }
  _mdl_jac_ir[ir_offset + 1] = jc_offset;
}

//==========================================================================

/*
 * Prg_SFunction.C -- class definition
 *
 */

/*
    Copyright (C) 1997--2003  Ruediger Franke

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

#include <Hxi_MEX_SFunction.h>
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

  _S = NULL;

  _mdl_nx = 0;
  _mdl_nu = 0;
  _mdl_ny = 0;

  _mdl_x0 = v_get(_mdl_nx);

  _ifList.append(new If_RealVec(GET_SET_CB(const VECP, mdl_x0)));
  _ifList.append(new If_String(GET_SET_CB(const char *, mdl_name)));
  _ifList.append(new If_String(GET_SET_CB(const char *, mdl_path)));
  _ifList.append(new If_String(GET_SET_CB(const char *, mdl_args)));
}

//--------------------------------------------------------------------------
Prg_SFunction::~Prg_SFunction()
{
  int i;
  if (_S) {
    mdlTerminate(_S);
    Hxi_SimStruct_destroy(_S);
  }
  for (i = 0; i < _mdl_nargs; i++)
    mxDestroyArray(_mx_args[i]);
  delete [] _mx_args;

  v_free(_mdl_x0);
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
    args[i] = Hxi::mx_parse_argument(str);
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
void Prg_SFunction::set_mdl_x0(const VECP value)
{
  v_copy_elements(value, _mdl_x0);
}

//--------------------------------------------------------------------------
void Prg_SFunction::setup_model()
{
  int i;

  // setup S-function
  if (_S) {
    mdlTerminate(_S);
    Hxi_SimStruct_destroy(_S);
  }

  _S = Hxi_SimStruct_create();

  // initialize model name and path
  ssSetModelName(_S, _mdl_name);
  if (_mdl_path[0] != '\0')
    ssSetPath(_S, _mdl_path);
  else
    ssSetPath(_S, _mdl_name);

  // initialize S-function parameters
  ssSetSFcnParamsCount(_S, _mdl_nargs);
  for (i = 0; i < _mdl_nargs; i++)
    ssSetSFcnParam(_S, i, _mx_args[i]);

  // initialize solver
  // (note: variable step size is indicated as the model must allow
  //  simulation time stepping back)
  ssSetVariableStepSolver(_S, 1);

  // initialize model
  mdlInitializeSizes(_S);
  if (ssGetErrorStatus(_S)) {
    fprintf(stderr, "Error from mdlInitializeSizes: %s\n",
	    ssGetErrorStatus(_S));
    ssSetErrorStatus(_S, NULL);
    m_error(E_FORMAT, "mdlInitializeSizes");
  }
  if (ssGetNumSFcnParams(_S) != ssGetSFcnParamsCount(_S)) {
    fprintf(stderr, "Parameter count mismatch: expected: %d, provided: %d\n",
	    ssGetNumSFcnParams(_S), ssGetSFcnParamsCount(_S));
    m_error(E_FORMAT, "Prg_SFunction::setup_model: parameter count mismatch");
  }

  // obtain model sizes
  _mdl_nx = ssGetNumContStates(_S);
  // count inputs of ports as long as they are contiguous in memory
  // (this limitation is because later on we will only access port 0)
  _mdl_nu = 0;
  for (i = 0; i < ssGetNumInputPorts(_S); i++) {
    if (*ssGetInputPortRealSignalPtrs(_S, i)
	== *ssGetInputPortRealSignalPtrs(_S, 0) + _mdl_nu)
      _mdl_nu += ssGetInputPortWidth(_S, i);
    else {
      m_warning(WARN_UNKNOWN,
		"Prg_SFunction::setup_model: ignoring non-contiguous inputs");
      break;
    }
  }
  // count outputs of ports as long as they are contiguous in memory
  // (this limitation is because later on we will only access port 0)
  _mdl_ny = 0;
  for (i = 0; i < ssGetNumOutputPorts(_S); i++) {
    if (ssGetOutputPortRealSignal(_S, i)
	== ssGetOutputPortRealSignal(_S, 0) + _mdl_ny)
      _mdl_ny += ssGetOutputPortWidth(_S, i);
    else {
      m_warning(WARN_UNKNOWN,
		"Prg_SFunction::setup_model: ignoring non-contiguous outputs");
      break;
    }
  }

  // initialize and check sample times of model
  mdlInitializeSampleTimes(_S);
  assert(ssGetNumSampleTimes(_S) > 0);
  assert(value(ssGetSampleTime(_S, 0)) == CONTINUOUS_SAMPLE_TIME);
  if (ssGetNumSampleTimes(_S) > 1)
    m_warning(WARN_UNKNOWN,
	      "Prg_SFunction::setup_model: ignoring multiple sample times");

  // start using S-function
  if (ssGetmdlStart(_S) != NULL)
    mdlStart(_S);

  // get initial states
  if (ssGetmdlInitializeConditions(_S) != NULL) {
    mdlInitializeConditions(_S);
    if (ssGetErrorStatus(_S)) {
      fprintf(stderr, "Error from mdlInitializeConditions: %s\n",
	      ssGetErrorStatus(_S));
      ssSetErrorStatus(_S, NULL);
      m_error(E_RANGE, "mdlInitializeConditions");
    }
  }
  real_T *mdl_x = ssGetContStates(_S);
  v_resize(_mdl_x0, _mdl_nx);
  for (i = 0; i < _mdl_nx; i++)
    _mdl_x0[i] = mdl_x[i];

  _mdl_needs_setup = false;
}


//==========================================================================

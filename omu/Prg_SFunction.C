/*
 * Prg_SFunction.C -- class definition
 *
 */

/*
    Copyright (C) 1997--2002  Ruediger Franke

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
#include <If_Method.h>

// redefine assert to throw an error instead of aborting
#undef assert
#define assert(expr) if (!(expr)) m_error(E_INTERN, "assert(" #expr ")");

typedef If_Method<Prg_SFunction> If_Cmd;

//--------------------------------------------------------------------------
Prg_SFunction::Prg_SFunction()
{
  _mdl_name = strdup("SFunction");
  _mdl_path = strdup("");

  // initialize _mdl_args and _mx_args
  _mdl_args = strdup("");
  _mdl_nargs = 0;
  _mx_args = NULL;

  _S = NULL;

  _mdl_nx = 0;
  _mdl_nu = 0;
  _mdl_ny = 0;

  _mdl_x0 = v_get(_mdl_nx);

  _ifList.append(new If_RealVec("mdl_x0", &_mdl_x0));
  _ifList.append(new If_Cmd("mdl_name", &Prg_SFunction::mdl_name, this));
  _ifList.append(new If_Cmd("mdl_path", &Prg_SFunction::mdl_path, this));
  _ifList.append(new If_Cmd("mdl_args", &Prg_SFunction::mdl_args, this));
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
int Prg_SFunction::mdl_name(int argc, char *argv[], char **result)
{
  if (argc == 1) {
    *result = _mdl_name;
  }
  else if (argc == 2) {
    free(_mdl_name);
    _mdl_name = strdup(argv[1]);
  }
  else {
    *result = "wrong # args, should be: mdl_name [new value]";
    return IF_ERROR;
  }
  return IF_OK;
}

//--------------------------------------------------------------------------
int Prg_SFunction::mdl_path(int argc, char *argv[], char **result)
{
  if (argc == 1) {
    *result = _mdl_path;
  }
  else if (argc == 2) {
    free(_mdl_path);
    _mdl_path = strdup(argv[1]);
  }
  else {
    *result = "wrong # args, should be: mdl_path [new value]";
    return IF_ERROR;
  }
  return IF_OK;
}

//--------------------------------------------------------------------------
int Prg_SFunction::mdl_args(int argc, char *argv[], char **result)
{
  if (argc == 1) {
    *result = _mdl_args;
  }
  else if (argc == 2) {
    const char *str, *str1;
    mxArray **args;
    int i, nargs;

    // parse args
    str = argv[1];
    str1 = Hxi::mx_count_columns(str, nargs);
    if (*str1 != '\0') {
      // we did not arrive at the end of the string
      *result = "failed to parse S-function args";
      return IF_ERROR;
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
    _mdl_args = strdup(argv[1]);
  }
  else {
    *result = "wrong # args, should be: mdl_args [new value]";
    return IF_ERROR;
  }
  return IF_OK;
}

//--------------------------------------------------------------------------
void Prg_SFunction::setup_sfun()
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
    m_error(E_FORMAT, "mdlInitializeSizes");
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
		"Prg_SFunction::setup_sfun: ignoring non-contiguous inputs");
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
		"Prg_SFunction::setup_sfun: ignoring non-contiguous outputs");
      break;
    }
  }

  // initialize and check sample times of model
  mdlInitializeSampleTimes(_S);
  assert(ssGetNumSampleTimes(_S) == 1);
  assert(value(ssGetSampleTime(_S, 0)) == CONTINUOUS_SAMPLE_TIME);

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
}


//==========================================================================

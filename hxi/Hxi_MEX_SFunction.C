/*
 * Hxi_MEX_SFunction.C:
 *   implementation of C MEX S-function interface for Hqp
 *
 * rf, 07/14/2001
 */

/*
    Copyright (C) 1994--2002  Ruediger Franke

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

#define Hxi_MEX_SFunction_C
#include "Hxi_MEX_SFunction.h"

#if defined(_MSC_VER) || defined(__MINGW32__)
#include <windows.h>
#else
#include <dlfcn.h>
#endif

// indices and size for MEX call argument
const int_T Hxi_NRHS = 4;	// number of right hand side arguments
const int_T Hxi_RHS_T = 0;	// index of time
const int_T Hxi_RHS_X = 1;	// index of states
const int_T Hxi_RHS_U = 2;	// index of inputs
const int_T Hxi_RHS_FLAG = 3;	// index of flag determining method to call

// prototype for MEX function
typedef void
(mexFunction_t)(int nlhs, mxArray *plhs[], int nrhs, mxArray *prhs[]);

//-------------------------------------------------------------------
// callback for allocating input ports
//-------------------------------------------------------------------
static int_T setNumInputPorts(void *arg, int_T nuPorts)
{
  SimStruct *S = (SimStruct *)arg;

  if (nuPorts < 0) {
    return 0;
  }

  _ssSetNumInputPorts(S, nuPorts);
  _ssSetSfcnUsesNumPorts(S, 1);

  if (nuPorts > 0) {
    struct _ssPortInputs *inputPorts = (struct _ssPortInputs*)
      mxCalloc(nuPorts, sizeof(struct _ssPortInputs));
    ssSetPortInfoForInputs(S, inputPorts);
  }

  return 1;
}

//-------------------------------------------------------------------
// callback for allocating output ports
//-------------------------------------------------------------------
static int_T setNumOutputPorts(void *arg, int_T nyPorts)
{
  SimStruct *S = (SimStruct *)arg;

  if (nyPorts < 0) {
    return 0;
  }

  _ssSetNumOutputPorts(S, nyPorts);
  _ssSetSfcnUsesNumPorts(S, 1);

  if (nyPorts > 0) {
    struct _ssPortOutputs *outputPorts = (struct _ssPortOutputs*)
      mxCalloc(nyPorts, sizeof(struct _ssPortOutputs));
    ssSetPortInfoForOutputs(S, outputPorts);
  }

  return 1;
}

//-------------------------------------------------------------------
// callback for setting input port dimensions
//-------------------------------------------------------------------
static int_T setInputPortDimensionInfo(SimStruct *S, int_T port, 
				       const DimsInfo_T *dimsInfo)
{
  // just propagate the with as if ssSetInputPortWidth had been called
  ssSetInputPortWidth(S, port, dimsInfo->width);
  return 1;
}

//-------------------------------------------------------------------
// callback for setting output port dimensions
//-------------------------------------------------------------------
static int_T setOutputPortDimensionInfo(SimStruct *S, int_T port, 
					const DimsInfo_T *dimsInfo)
{
  // just propagate the with as if ssSetOutputPortWidth had been called
  ssSetOutputPortWidth(S, port, dimsInfo->width);
  return 1;
}

//-------------------------------------------------------------------
// callback for allocating data work records
//-------------------------------------------------------------------
static int_T setNumDWork(void *arg, int_T nDWork)
{
  SimStruct *S = (SimStruct *)arg;

  if (nDWork < 0) {
    return 0;
  }

  if (nDWork > 0) {
    // allocate nDWork data work records
    struct _ssDWorkRecord *dWork = (struct _ssDWorkRecord *)
      mxCalloc(nDWork, sizeof(struct _ssDWorkRecord));
    
    ssSetSFcnDWork(S, dWork);
  }
  _ssSetNumDWork(S, nDWork);

  return 1;
}

//-------------------------------------------------------------------
SimStruct *Hxi_SimStruct_create()
{
  // allocate a SimStruct
  SimStruct *S = (SimStruct *)mxCalloc(1, sizeof(SimStruct));

  // initialize the SimStruct as root model
  ssSetModelName(S, "SimStruct");
  ssSetPath(S, "SimStruct");
  ssSetRootSS(S, S);

  // allocate a model info struct
  struct _ssMdlInfo *mdlInfo = (struct _ssMdlInfo *)
    mxCalloc(1, sizeof(struct _ssMdlInfo));
  ssSetMdlInfoPtr(S, mdlInfo);

  // allocate model methods 2 struct
  struct _ssSFcnModelMethods2 *modelMethods2 = (struct _ssSFcnModelMethods2 *)
    mxCalloc(1, sizeof(struct _ssSFcnModelMethods2));
  ssSetModelMethods2(S, modelMethods2);

  // allocate sparse header for Jacobian
  SparseHeader *jacobian = (SparseHeader *)
    mxCalloc(1, sizeof(SparseHeader));
  ssSetJacobianHeader(S, jacobian);

  // register callback methods
  ssSetRegNumInputPortsFcn(S, setNumInputPorts);
  ssSetRegNumInputPortsFcnArg(S, S);
  ssSetRegNumOutputPortsFcn(S, setNumOutputPorts);
  ssSetRegNumOutputPortsFcnArg(S, S);
  ssSetRegInputPortDimensionInfoFcn(S, setInputPortDimensionInfo);
  ssSetRegOutputPortDimensionInfoFcn(S, setOutputPortDimensionInfo);
  ssSetNumDWorkFcn(S, setNumDWork);

  return S;
}

//-------------------------------------------------------------------
void Hxi_SimStruct_destroy(SimStruct *S)
{
  if (S == NULL)
    return;

  // obtain handle and release MEX S-function
#if defined(_MSC_VER) || defined(__MINGW32__)
  HMODULE handle = GetModuleHandle(ssGetPath(S));
  if (handle)
    FreeLibrary(handle);
#else
  void *handle = dlopen(ssGetPath(S), RTLD_LAZY);
  if (handle)
    dlclose(handle);
#endif

  // free memory
  mxFree(ssGetTPtr(S));
  mxFree(ssGetIWork(S));
  mxFree(ssGetPWork(S));
  mxFree(ssGetRWork(S));
  mxFree(ssGetSFcnDWork(S));
  mxFree(ssGetPortInfoForInputs(S));
  mxFree(ssGetPortInfoForOutputs(S));
  mxFree(ssGetSFcnParamsPtr(S));
  mxFree(ssGetJacobianHeader(S));
  mxFree(ssGetModelMethods2(S));
  mxFree(ssGetMdlInfoPtr(S));
  mxFree(S);
}

//-------------------------------------------------------------------
void mdlInitializeSizes(SimStruct *S)
{
  int i, ip;

  // get handle to MEX S-function
#if defined(_MSC_VER) || defined(__MINGW32__)
  HMODULE handle = LoadLibrary(ssGetPath(S));
  if (!handle) {
    ssSetErrorStatus(S, "LoadLibrary failed");
    return;
  }
#else
  void *handle = dlopen(ssGetPath(S), RTLD_LAZY);
  if (!handle) {
    ssSetErrorStatus(S, dlerror());
    return;
  }
#endif

  // get pointer to entry point mexFunction
  mexFunction_t *mexFunction_p;
#if defined(_MSC_VER) || defined(__MINGW32__)
  mexFunction_p = (mexFunction_t *)GetProcAddress(handle, "mexFunction");
  if (!mexFunction_p) {
    ssSetErrorStatus(S, "GetProcAddress failed for mexFunction");
    return;
  }
#else
  mexFunction_p = (mexFunction_t *)dlsym(handle, "mexFunction");
  if (!mexFunction_p) {
    ssSetErrorStatus(S, dlerror());
    return;
  }
#endif
  //
  // initialize SimStruct via MEX calling interface
  //   - get pointers to S-function methods
  //   - get model sizes and check parameters
  //     via call to mdlInitializeSizes
  //

  // prepare MEX calling arguments
  mxArray *plhs[1];
  plhs[0] = NULL;

  mxArray *prhs[Hxi_NRHS];
  for (i = 0; i < Hxi_NRHS; i++)
    prhs[i] = NULL;

  // pass pointer to SimStruct in argument Hxi_RHS_X
  int m = sizeof(SimStruct *)/sizeof(int_T) + 1;
  prhs[Hxi_RHS_X] = mxCreateDoubleMatrix(m, 1, mxREAL);
  real_T *pr = (real_T*)mxGetPr(prhs[Hxi_RHS_X]);
  int_T *intS = (int_T *)&S;
  for (i = 0; i < m-1; i++)
    pr[i] = (real_T)intS[i];
  // indicate level 2 S-function to get pointers to S-function methods
  pr[m-1] = SIMSTRUCT_VERSION_LEVEL2;

  // pass flag 0 for initialization as argument Hxi_RHS_FLAG
  prhs[Hxi_RHS_FLAG] = mxCreateDoubleMatrix(1, 1, mxREAL);
  pr = (real_T*)mxGetPr(prhs[Hxi_RHS_FLAG]);
  pr[0] = 0;

  // call mexFunction for initialization
  (*mexFunction_p)(1, plhs, Hxi_NRHS, prhs);

  // free call arguments
  mxDestroyArray(plhs[0]);
  for (i = 0; i < Hxi_NRHS; i++)
    mxDestroyArray(prhs[i]);

  // return if S-function's mdlInitializeSizes failed
  if (ssGetErrorStatus(S) != NULL
      || ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S))
    return;

  //
  // allocate memory for SimStruct
  //

  int nuPorts = ssGetNumInputPorts(S);
  int nyPorts = ssGetNumOutputPorts(S);
  int nx = ssGetNumContStates(S);
  int nxd = ssGetNumDiscStates(S);
  int nrw = ssGetNumRWork(S);
  int npw = ssGetNumPWork(S);
  int niw = ssGetNumIWork(S);
  int nst = ssGetNumSampleTimes(S);
  int ndw = ssGetNumDWork(S);

  // count inputs and initialize ports
  int nu = 0;
  if (nuPorts > 0 && &ssGetInputPortWidth(S,0) == NULL) {
    // this seems to be a level 1 S-function;
    // we allocate NumInputs ports of width one
    // (sensible as NumInputPorts is stored at the same place as NumInputs)
    nu = nuPorts;
    setNumInputPorts(S, nuPorts);
    for (ip = 0; ip < nuPorts; ip++)
      ssSetInputPortWidth(S, ip, 1);
  }
  else {
    // level 2: count number of inputs of all ports,
    for (ip = 0; ip < ssGetNumInputPorts(S); ip++) {
      // set dynamically sized ports to width one
      if (ssGetInputPortWidth(S, ip) == DYNAMICALLY_SIZED)
	ssSetInputPortWidth(S, ip, 1);
      nu += ssGetInputPortWidth(S, ip);
    }
  }
    
  // count outputs and ports
  int ny = 0;
  if (nyPorts > 0 && &ssGetOutputPortWidth(S,0) == NULL) {
    // this seems to be a level 1 S-function;
    // we allocate NumOutputs ports of width one
    // (sensible as NumOutputPorts is stored at the same place as NumOutputs)
    ny = nyPorts;
    setNumOutputPorts(S, nyPorts);
    for (ip = 0; ip < nyPorts; ip++)
      ssSetOutputPortWidth(S, ip, 1);
  }
  else {
    // level 2: count number of outputs of all ports,
    for (ip = 0; ip < ssGetNumOutputPorts(S); ip++) {
      // set dynamically sized ports to width one
      if (ssGetOutputPortWidth(S, ip) == DYNAMICALLY_SIZED)
	ssSetOutputPortWidth(S, ip, 1);
      ny += ssGetOutputPortWidth(S, ip);
    }
  }

  // count number of dwork array elements
  int ndwel = 0;
  for (i = 0; i < ndw; i++)
    ndwel += ssGetDWorkWidth(S, i);

  // initialize Jacobian and obtain number of elements
  int njacm = nx + nxd + ny;
  int njacn = nx + nxd + nu;
  int njacel = ssGetJacobianNzMax(S);
  if (njacel < 0) {
    // full Jacobian
    njacel = njacm * njacn;
  }
  int njacjc = 0;
  if (njacel > 0) {
    njacjc = njacn + 1;
    ssSetJacobianM(S, njacm);
    ssSetJacobianN(S, njacn);
  }

  // allocate memory for real_T
  mxFree(ssGetRWork(S));
  real_T *rptr = (real_T *)mxCalloc(sizeof(real_T),
				    nrw + nu + ny + nxd + 2*nx + ndwel
				    + njacel);
  ssSetRWork(S, rptr);
  rptr += nrw;

  // allocate memory for pointers
  mxFree(ssGetPWork(S));
  void **pptr = (void **)mxCalloc(sizeof(void *), npw + nu);
  ssSetPWork(S, pptr);
  pptr += npw;

  // setup pointers to inputs
  for (i = 0; i < nu; i++)
    pptr[i] = &rptr[i];

  // setup inputs
  ssSetU(S, rptr);  // level 1 S-function
  for (ip = 0; ip < nuPorts; ip++) {
    if (ssGetInputPortRequiredContiguous(S, ip))
      // setup input signals
      ssSetInputPortSignal(S, ip, rptr);
    else
      // setup pointer to input signals
      ssSetInputPortSignalPtrs(S, ip, (InputPtrsType)pptr);
    rptr += ssGetInputPortWidth(S, ip);
    pptr += ssGetInputPortWidth(S, ip);
  }

  // setup outputs
  ssSetY(S, rptr);  // level 1 S-function
  for (ip = 0; ip < nyPorts; ip++) {
    // setup pointers to outputs in SimStruct
    ssSetOutputPortSignal(S, ip, rptr);
    rptr += ssGetOutputPortWidth(S, ip);
  }

  // setup states
  ssSetDiscStates(S, rptr);
  rptr += nxd;
  ssSetContStates(S, rptr);
  rptr += nx;
  ssSetdX(S, rptr);
  rptr += nx;

  // setup dwork arrays
  for (i = 0; i < ndw; i++) {
    ssSetDWork(S, i, rptr);
    rptr += ssGetDWorkWidth(S, i);
  }

  // setup Jacobian elements
  ssSetJacobianPr(S, rptr);
  rptr += njacel;

  // allocate memory for int_T
  mxFree(ssGetIWork(S));
  int_T *iptr = (int_T *)mxCalloc(sizeof(int_T), niw + 2*nst
				  + njacel + njacjc);
  ssSetIWork(S, iptr);
  iptr += niw;

  // allocate and setup times
  ssSetSampleHitPtr(S, iptr);
  iptr[0] = 1; // required for ssIsContinuousTask(0) == true
  iptr += nst;
  ssSetSampleTimeTaskIDPtr(S, iptr);
  for (i = 0; i < nst; i++)
    *iptr++ = i; // task id and sample time index are identical

  mxFree(ssGetTPtr(S));
  time_T *tptr = (time_T *)mxCalloc(sizeof(time_T), 3*nst);
  ssSetTPtr(S, tptr);
  tptr += nst;
  ssSetSampleTimePtr(S, tptr);
  tptr += nst;
  ssSetOffsetTimePtr(S, tptr);

  // setup Jacobian
  ssSetJacobianIr(S, iptr);
  iptr += njacel;
  ssSetJacobianJc(S, iptr);
  iptr += njacjc;
}


//===================================================================

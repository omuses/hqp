/*
 * Hxi_MEX_SFunction.h:
 *   declarations of C MEX S-function interface for Hqp
 *
 * rf, 07/14/2001
 */

/*
    Copyright (C) 1994--2001  Ruediger Franke

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

#if !defined(Hxi_MEX_SFunction_H)
#define Hxi_MEX_SFunction_H

#include <assert.h>

#if !defined(MATLAB_MEX_FILE)
// MATLAB_MEX_FILE must be defined prior to inclusion of simstruc.h
#define MATLAB_MEX_FILE 1
#endif

#include <simstruc.h>

/** Create SimStruct for S-function. */
SimStruct *Hxi_SimStruct_create();

/** Realease S-function. */
void Hxi_SimStruct_destroy(SimStruct *S);

/** Initialize sizes of data vectors. This function loads the shared object 
    found under ssGetPath(S), calls the S-function method mdlInitializeSizes
    and allocates memory required for a level 2 S-function.
    The models parameters must be initialized prior to calling this
    function. */
void mdlInitializeSizes(SimStruct *S);

#if !defined(Hxi_MEX_SFunction_C)

// redefine ssSetSFcnParamsCount to allocate memory for mxArray pointers
#undef ssSetSFcnParamsCount
#define ssSetSFcnParamsCount(S,n) \
  _ssSetSFcnParamsCount(S,n); \
  mxFree(ssGetSFcnParamsPtr(S)); \
  ssSetSFcnParamsPtr(S, (const mxArray **)mxCalloc(sizeof(mxArray *), n))

// disable macros that are critical for memory management
#undef ssSetRWork
#define ssSetRWork(S, rwork) ssSetRWork_cannot_be_used_with_Hxi_MEX_SFunction

#undef ssSetPWork
#define ssSetPWork(S, pwork) ssSetPWork_cannot_be_used_with_Hxi_MEX_SFunction

#undef ssSetIWork
#define ssSetIWork(S, iwork) ssSetIWork_cannot_be_used_with_Hxi_MEX_SFunction

#undef ssSetTPtr
#define ssSetTPtr(S, tptr) ssSetTPtr_cannot_be_used_with_Hxi_MEX_SFunction

#endif

/** Optional: Allocate local ressources for simulation. */
inline void mdlStart(SimStruct *S)
{
  assert(ssGetmdlStart(S) != NULL);
  sfcnStart(S);
}

/** Initialize sample times. */
inline void mdlInitializeSampleTimes(SimStruct *S)
{
  sfcnInitializeSampleTimes(S);
}

/** Optional: Compute initial conditions. */
inline void mdlInitializeConditions(SimStruct *S)
{
  assert(ssGetmdlInitializeConditions(S) != NULL);
  if (ssGetVersion(S) == SIMSTRUCT_VERSION_LEVEL1)
    sfcnInitializeConditionsLevel1(ssGetX(S), S);
  else
    sfcnInitializeConditions(S);
}

/** Compute model outputs. */
inline void mdlOutputs(SimStruct *S, int_T tid)
{
  if (ssGetVersion(S) == SIMSTRUCT_VERSION_LEVEL1)
    sfcnOutputsLevel1(ssGetY(S), ssGetX(S), ssGetU(S), S, tid);
  else
    sfcnOutputs(S, tid);
}

/** Optional: Update discrete-time states. */
inline void mdlUpdate(SimStruct *S, int_T tid)
{
  assert(ssGetmdlUpdate(S) != NULL);
  if (ssGetVersion(S) == SIMSTRUCT_VERSION_LEVEL1)
    sfcnUpdateLevel1(ssGetX(S), ssGetU(S), S, tid);
  else
    sfcnUpdate(S, tid);
}

/** Optional: Compute derivatives for continuous-time states. */
inline void mdlDerivatives(SimStruct *S)
{
  assert(ssGetmdlDerivatives(S) != NULL);
  if (ssGetVersion(S) == SIMSTRUCT_VERSION_LEVEL1)
    sfcnDerivativesLevel1(ssGetdX(S), ssGetX(S), ssGetU(S), S, 0);
  else
    sfcnDerivatives(S);
}

/** Release resources allocated for simulation. */
inline void mdlTerminate(SimStruct *S)
{
  sfcnTerminate(S);
}

#endif // !defined(Hxi_MEX_SFunction_H)

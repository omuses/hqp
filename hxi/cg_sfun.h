/**
 * @file cg_sfun.h
 * Native S-function code generation include file for Hqp.
 * Currently macros for checking optional S-function methods are defined.
 * This allows to write the same code for an inlined S-function and for 
 * an external MEX S-function.
 *
 * (Simulink is a registered trademark of The MathWorks, Inc.)
 *
 * rf, 05/05/2001
 */

/*
    Copyright (C) 1994--2005  Ruediger Franke

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

#define HXI_STRINGIFY(token) #token

#if !defined(S_FUNCTION_LEVEL) || S_FUNCTION_LEVEL < 2
#error "Hxi requires S_FUNCTION_LEVEL 2"
#endif

/**
 * Entry function into Hxi S-function. It initializes a
 * %SimStruct with S-function name, version, and with
 * pointers to S-function methods.
 */
#if defined(HXI_INLINE_S_FUNCTION)
static
#else
#  if defined(__cplusplus)
   extern "C"
#  endif
#  if defined(_MSC_VER)
   __declspec(dllexport)
#  endif
#endif
void Hxi_SimStruct_init(SimStruct *S) {
  // S-function name and level
  ssSetModelName(S, HXI_STRINGIFY(S_FUNCTION_NAME));
  ssSetVersion(S, S_FUNCTION_LEVEL);

  // required S-function methods
  ssSetmdlInitializeSizes(S, &mdlInitializeSizes);
  ssSetmdlInitializeSampleTimes(S, &mdlInitializeSampleTimes);
  ssSetmdlOutputs(S, &mdlOutputs);
  ssSetmdlTerminate(S, &mdlTerminate);

  // optional S-function methods
#if defined(MDL_CHECK_PARAMETERS)
  ssSetmdlStart(S, &mdlCheckParameters);
#endif

#if defined(MDL_START)
  ssSetmdlStart(S, &mdlStart);
#endif

#if defined(MDL_INITIALIZE_CONDITIONS)
  ssSetmdlInitializeConditions(S, &mdlInitializeConditions);
#endif

#if defined(MDL_UPDATE)
  ssSetmdlUpdate(S, &mdlUpdate);
#endif

#if defined(MDL_DERIVATIVES)
  ssSetmdlDerivatives(S, &mdlDerivatives);
#endif

#if defined(MDL_JACOBIAN)
  ssSetmdlJacobian(S, &mdlJacobian);
#endif
}

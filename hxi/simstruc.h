/**
 * @file simstruc.h
 *   S-function SimStruct include file for treating an S-function either
 *   as external MEX object or inlined with Hqp using the Hxi interface.
 *
 * Intended usage:
 *  An S-function including this file can either be compiled to a 
 *  stand-alone MEX S-function or it can be inlined with an Hqp 
 *  optimization problem.
 *  - Define the macro HXI_MEX_S_FUNCTION if the S-function should be
 *  compiled to a MATLAB MEX S-function. Then the original
 *  %simulink/%include/%simstruc.h found in the path will be included.
 *  - The alternative file Hxi_SimStruct.h will be included to compile
 *  the S-function with Hqp.
 *
 * (MATLAB and Simulink are registered trademarks of The MathWorks, Inc.)
 *
 * rf, 05/05/2001
 */

#if !defined(HXI_SIMSTRUC_H)
#define HXI_SIMSTRUC_H

#if defined(HXI_MEX_S_FUNCTION)
/* The S-function is being compiled as MATLAB MEX file;
 * we include the simstruc.h found under simulink/include. */

/* Note: define MathWorks_h to get in simulink/include/simstruc.h the 
 * use SS_SL_INTERNAL. This enables macros for setting up the SimStruct 
 * passed to model methods (introduced with Matlab 6.5). */
#define MathWorks_h // Matlab 6.5
#define SL_INTERNAL // Matlab 7.0
#include <simulink/include/simstruc.h>
#undef MathWorks_h

# if !defined(ssGetNumInputPorts)
# define ssGetNumInputPorts(S) 	(S)->sizes.in.numInputPorts
# endif

# if !defined(ssGetNumOutputPorts)
# define ssGetNumOutputPorts(S) (S)->sizes.out.numOutputPorts
# endif

# if !defined(ssSetMinorTimeStep)
# define ssSetMinorTimeStep(S, flag) \
  ssSetSimTimeStep(S, (flag)? MINOR_TIME_STEP: MAJOR_TIME_STEP)
# endif

# if defined(ssSetInputPortVectorDimension)
# undef ssSetInputPortVectorDimension
# endif
# define ssSetInputPortVectorDimension 	ssSetInputPortWidth

# if defined(ssSetOutputPortVectorDimension)
# undef ssSetOutputPortVectorDimension
# endif
# define ssSetOutputPortVectorDimension ssSetOutputPortWidth

// guess MATLAB version (times 10 to avoid decimal point)
#if defined(SS_OPTION_ASYNCHRONOUS_INTERRUPT)
#define MATLAB_VERSION 65
#elif defined(SS_OPTION_ASYNCHRONOUS_CUSTOM)
#define MATLAB_VERSION 61
#else
#define MATLAB_VERSION 60
#endif

#else
// We are compiling the S-function for Hqp.

/**
 * Hqp eXternal Interfaces: native implementations of types that are 
 * defined by external packages.
 * For example an S-function can be accessed via a MEX interface using
 * the MEX %SimStruct or it can be inlined with Hqp using Hxi::SimStruct.
 */

#include "Hxi_SimStruct.h"

#endif /* defined(HXI_MEX_S_FUNCTION) */

#endif /* !defined(HXI_SIMSTRUC_H) */

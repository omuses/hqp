/**
 * @file simstruc.h
 *   S-function SimStruct include file for treating an S-function either
 *   as external MEX object or inlined with Hqp using the Hxi interface.
 *
 * Intended usage:
 *  An S-function including this file can either be compiled to a 
 *  stand-alone MEX S-function or it can be inlined with an Hqp 
 *  optimization problem.
 *  - Define the macro MATLAB_MEX_FILE (as when invoking the MATLAB(R) 
 *  command mex) if the S-function should be compiled to a MEX S-function.
 *  Then the original %simulink/%include/%simstruc.h found in the path
 *  will be included.
 *  - The alternative file Hxi_SimStruct.h will be included to compile
 *  the S-function with Hqp.
 *
 * (MATLAB and Simulink are registered trademarks of The MathWorks, Inc.)
 *
 * rf, 05/05/2001
 */

#if !defined(HXI_SIMSTRUC_H)
#define HXI_SIMSTRUC_H

#if defined(MATLAB_MEX_FILE)
// The S-function is being compiled as Matlab MEX file;
// we include the simstruc.h found under simulink/include.

// Note: define MathWorks_h to get in simulink/include/simstruc.h the 
// use SS_SL_INTERNAL. This enables macros for setting up the SimStruct 
// passed to model methods (introduced with Matlab 6.5).
#undef MATLAB_MEX_FILE
#define MathWorks_h
#include <simulink/include/simstruc.h>
#undef MathWorks_h
#define MATLAB_MEX_FILE

# if !defined(ssGetNumInputPorts)
# define ssGetNumInputPorts(S) _ssGetNumInputPorts(S)
# endif

# if !defined(ssGetNumOutputPorts)
# define ssGetNumOutputPorts(S) _ssGetNumOutputPorts(S)
# endif

// guess MATLAB version (times 10 to avoid decimal point)
#if defined(SS_OPTION_ASYNCHRONOUS_INTERRUPT)
#define MATLAB_VERSION 65
#elif defined(SS_OPTION_ASYNCHRONOUS_CUSTOM)
#define MATLAB_VERSION 61
#else
#define MATLAB_VERSION 60
#endif

/** ADOL-C's value function for double. It may be needed to compile
    code written for real_T adouble with real_T double. */
inline double value(double a) {return a;}

#else
// We are compiling the S-function for Hqp.
// ADOL-C is used for automatic differentiation per default.

/**
 * Hqp eXternal Interfaces: native implementations of types that are 
 * defined by external packages.
 * For example an S-function can be accessed via a MEX interface using
 * the MEX %SimStruct or it can be inlined with Hqp using Hxi::SimStruct.
 */
namespace Hxi {};

#include "Hxi_SimStruct.h"

// make Hxi implementations visible in global scope
using namespace Hxi;

#endif // defined(MATLAB_MEX_FILE)

#endif // !defined(HXI_SIMSTRUC_H)

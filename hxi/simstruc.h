/**
 * @file simstruc.h
 *   S-function %SimStruct include file for compiling either for Simulink(R) 
 *   or with HQP.
 *
 * Intended usage:
 *  An S-function including this file can either be compiled to a 
 *  stand-alone MEX S-function or it can be inlined with an HQP 
 *  optimization problem.
 *  - Define the macro MATLAB_MEX_FILE (as when invoking the MATLAB(R) 
 *  command mex) if the S-function should be compiled to a MEX S-function.
 *  Then the original simulink/include/simstruc.h found in the path
 *  will be included.
 *  - The alternative file Hxi_SimStruct.h will be included to compile
 *  the S-function with HQP.
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

#include <simulink/include/simstruc.h>

/** ADOL-C's value function for double. It may be needed to compile
    code written for real_T adouble with real_T double. */
inline double value(double a) {return a;}

#else
// We are compiling the S-function for HQP.
// ADOL-C is used for automatic differentiation per default.

#include "Hxi_SimStruct.h"

#endif // defined(MATLAB_MEX_FILE)

#endif // !defined(HXI_SIMSTRUC_H)

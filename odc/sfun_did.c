/*
 *  sfun_did.c: S-function for a discrete-time double integrator
 *
 *  rf, 05/05/2001
 *
 *  derived from:
 *
 *  File    : dsfunc.c
 *  Abstract:
 *
 *      Example C-file S-function for defining a discrete system.  
 *      For more details about S-functions, see simulink/src/sfuntmpl_doc.c.
 * 
 *  Copyright 1990-2000 The MathWorks, Inc.
 *  $Revision: 1.1 $
 */

#define S_FUNCTION_NAME sfun_did
#define S_FUNCTION_LEVEL 2

#include "simstruc.h"

#define U(element) (*uPtrs[element])  /* Pointer to Input Port0 */

/*====================*
 * S-function methods *
 *====================*/

/* Function: mdlInitializeSizes ===============================================
 * Abstract:
 *    The sizes information is used by Simulink to determine the S-function
 *    block's characteristics (number of inputs, outputs, states, etc.).
 */
static void mdlInitializeSizes(SimStruct *S)
{
    ssSetNumSFcnParams(S, 0);  /* Number of expected parameters */
    if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)) {
        return; /* Parameter mismatch will be reported by Simulink */
    }

    ssSetNumContStates(S, 0);
    ssSetNumDiscStates(S, 2);

    if (!ssSetNumInputPorts(S, 1)) return;
    ssSetInputPortWidth(S, 0, 1);
    ssSetInputPortDirectFeedThrough(S, 0, 0);

    if (!ssSetNumOutputPorts(S, 1)) return;
    ssSetOutputPortWidth(S, 0, 2);

    ssSetNumSampleTimes(S, 1);
    ssSetNumRWork(S, 0);
    ssSetNumIWork(S, 0);
    ssSetNumPWork(S, 0);
    ssSetNumModes(S, 0);
    ssSetNumNonsampledZCs(S, 0);

    /* Take care when specifying exception free code - see sfuntmpl_doc.c */
    ssSetOptions(S, SS_OPTION_EXCEPTION_FREE_CODE);
}



/* Function: mdlInitializeSampleTimes =========================================
 * Abstract:
 *    Specifiy that we inherit our sample time from the driving block.
 */
static void mdlInitializeSampleTimes(SimStruct *S)
{
    ssSetSampleTime(S, 0, 1.0/60);
    ssSetOffsetTime(S, 0, 0.0);
}



#define MDL_INITIALIZE_CONDITIONS
/* Function: mdlInitializeConditions ========================================
 * Abstract:
 *    Initialize discrete states.
 */
static void mdlInitializeConditions(SimStruct *S)
{
    real_T *x0 = ssGetRealDiscStates(S);

    x0[0] = 1.0;
    x0[1] = 0.0;
}



/* Function: mdlOutputs =======================================================
 * Abstract:
 *      y = f(x)
 */
static void mdlOutputs(SimStruct *S, int_T tid)
{
    real_T            *y    = ssGetOutputPortRealSignal(S,0);
    real_T            *x    = ssGetRealDiscStates(S);
 
    UNUSED_ARG(tid); /* not used in single tasking mode */

    /* y = f(x) */
    y[0] = x[0];
    y[1] = x[1];
}



#define MDL_UPDATE
/* Function: mdlUpdate ======================================================
 * Abstract:
 *      x(k+1) = f(x(k),u(k))
 */
static void mdlUpdate(SimStruct *S, int_T tid)
{
    real_T            dt       = ssGetSampleTime(S, 0);
    real_T            tempX[2] = {0.0, 0.0};
    real_T            *x       = ssGetRealDiscStates(S);
    InputRealPtrsType uPtrs    = ssGetInputPortRealSignalPtrs(S,0);

    UNUSED_ARG(tid); /* not used in single tasking mode */

    /* x(k+1) = f(x(k),u(k)) */
    tempX[0] = x[0] + U(0)*dt;
    tempX[1] = x[1] + x[0]*dt + U(0)*0.5*dt*dt;
 
    x[0] = tempX[0];
    x[1] = tempX[1];
}



/* Function: mdlTerminate =====================================================
 * Abstract:
 *    No termination needed, but we are required to have this routine.
 */
static void mdlTerminate(SimStruct *S)
{
    UNUSED_ARG(S); /* unused input argument */
}

#ifdef  MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "simulink.c"      /* MEX-file interface mechanism */
#else
#include "cg_sfun.h"       /* Code generation registration function */
#endif


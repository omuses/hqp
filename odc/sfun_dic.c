/*
 *  sfun_dic.c: S-function for a continuous-time double integrator
 *
 *  rf, 05/06/2001
 *
 *  derived from:
 *
 *  File    : csfunc.c
 *  Abstract:
 *
 *      Example C-file S-function for defining a continuous system.  
 *      For more details about S-functions, see simulink/src/sfuntmpl_doc.c.
 * 
 *  Copyright 1990-2000 The MathWorks, Inc.
 *  Revision: 1.7 
 */

#define S_FUNCTION_NAME sfun_dic
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

    ssSetNumContStates(S, 2);
    ssSetNumDiscStates(S, 0);

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

    /* Setup Jacobian if not used with Hxi::SimStruct as
       Hxi::SimStruct does not support mdlJacobian, but works with ADOL-C. */
#if !defined(Hxi_SimStruct_H)
    ssSetJacobianNzMax(S, 4);
#endif

    /* Take care when specifying exception free code - see sfuntmpl_doc.c */
    ssSetOptions(S, SS_OPTION_EXCEPTION_FREE_CODE);
}



/* Function: mdlInitializeSampleTimes =========================================
 * Abstract:
 *    Specifiy that we have a continuous sample time.
 */
static void mdlInitializeSampleTimes(SimStruct *S)
{
    ssSetSampleTime(S, 0, CONTINUOUS_SAMPLE_TIME);
    ssSetOffsetTime(S, 0, 0.0);
}



#define MDL_INITIALIZE_CONDITIONS
/* Function: mdlInitializeConditions ========================================
 * Abstract:
 *    Initialize continuous states.
 */
static void mdlInitializeConditions(SimStruct *S)
{
    real_T *x0 = ssGetContStates(S);

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
    real_T            *x    = ssGetContStates(S);
 
    UNUSED_ARG(tid); /* not used in single tasking mode */

    /* y = f(x) */
    y[0] = x[0];
    y[1] = x[1];
}



#define MDL_DERIVATIVES
/* Function: mdlDerivatives =================================================
 * Abstract:
 *      dx = f(x,u)
 */
static void mdlDerivatives(SimStruct *S)
{
    real_T            *dx   = ssGetdX(S);
    real_T            *x    = ssGetContStates(S);
    InputRealPtrsType uPtrs = ssGetInputPortRealSignalPtrs(S,0);

    /* dx = f(x,u) */
    dx[0] = U(0);
    dx[1] = x[0];
}



#if !defined(Hxi_SimStruct_H)
#define MDL_JACOBIAN
/* Function: mdlJacobian =================================================
 * Abstract:
 *    J = d(dxc,xd,y)/d(xc,xd,u)
 *	=	x1	x2	u
 *	dx1	0	0	1
 *	dx2	1	0	0
 *	y1	1	0	0
 *	y2	0	1	0
 */
static void mdlJacobian(SimStruct *S)
{
  real_T *pr; /* Jacobian elements */
  int_T *ir;  /* row indices */
  int_T *jc;  /* start index for each column */
  int_T j;    /* column number */
  int_T idx;  /* index into data vectors pr and ir */

  pr = ssGetJacobianPr(S);
  ir = ssGetJacobianIr(S);
  jc = ssGetJacobianJc(S);
  j = 0;
  idx = 0;
  /* first column */
  jc[j++] = idx;
  ir[idx] = 1; pr[idx] = 1.0; idx++;
  ir[idx] = 2; pr[idx] = 1.0; idx++;
  /* second column */
  jc[j++] = idx;
  ir[idx] = 3; pr[idx] = 1.0; idx++;
  /* third column */
  jc[j++] = idx;
  ir[idx] = 0; pr[idx] = 1.0; idx++;
  /* end marker */
  jc[j] = idx;
}
#endif



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

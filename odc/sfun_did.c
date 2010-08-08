/*
 *  sfun_did.c: S-function for a discrete-time double integrator
 *
 *  rf, 05/05/2001
 *
 *  derived from:
 *
 *  File    : dsfunc.c (treatment of parameters: sfun_zc_sat.c)
 *  Abstract:
 *
 *      Example C-file S-function for defining a discrete system.  
 *      For more details about S-functions, see simulink/src/sfuntmpl_doc.c.
 * 
 *  Copyright 1990-2000 The MathWorks, Inc.
 *  Revision: 1.10 
 */

#define S_FUNCTION_NAME sfun_did
#define S_FUNCTION_LEVEL 2

#include "simstruc.h"

/*========================*
 * General Defines/macros *
 *========================*/

/* index to dt parameter */
#define I_PAR_DT                  0

/* total number of block parameters */
#define N_PAR                     1

/*
 *  Make access to data more readable.
 */
#define P_PAR_DT    ssGetSFcnParam(S, I_PAR_DT)

#define U(element) (*uPtrs[element])  /* Pointer to Input Port0 */

/*====================*
 * S-function methods *
 *====================*/

#define     MDL_CHECK_PARAMETERS
#if defined(MDL_CHECK_PARAMETERS)

  /* Function: mdlCheckParameters =============================================
   * Abstract:
   *   Check that parameter choices are allowable.
   */
  static void mdlCheckParameters(SimStruct *S)
  {
      int_T      i;
      int_T      numDT;

      /*
       * check parameter basics
       */
      for ( i = 0; i < N_PAR; i++ ) {
          if ( mxIsEmpty(    ssGetSFcnParam(S,i) ) ||
               mxIsSparse(   ssGetSFcnParam(S,i) ) ||
               mxIsComplex(  ssGetSFcnParam(S,i) ) ||
               !mxIsNumeric( ssGetSFcnParam(S,i) ) ) {
              ssSetErrorStatus(S, "Parameters must be real vectors.");
              return;
          }
      }

      /*
       * Check sizes of parameters.
       */
      numDT = mxGetNumberOfElements(P_PAR_DT);

      if ( numDT != 1 ) {
          ssSetErrorStatus(S, "Parameter dt must have size one.");
	  return;
      }
  }
#endif /* MDL_CHECK_PARAMETERS */



/* Function: mdlInitializeSizes ===============================================
 * Abstract:
 *    The sizes information is used by Simulink to determine the S-function
 *    block's characteristics (number of inputs, outputs, states, etc.).
 */
static void mdlInitializeSizes(SimStruct *S)
{
    ssSetNumSFcnParams(S, N_PAR);  /* Number of expected parameters */

    if (ssGetNumSFcnParams(S) == ssGetSFcnParamsCount(S)) {
        mdlCheckParameters(S);
        if (ssGetErrorStatus(S) != NULL) {
            return;
        }
    } else {
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

    /* Setup Jacobian */
    ssSetJacobianNzMax(S, 7);

    /* Take care when specifying exception free code - see sfuntmpl_doc.c */
    ssSetOptions(S, SS_OPTION_EXCEPTION_FREE_CODE);
}



/* Function: mdlInitializeSampleTimes =========================================
 * Abstract:
 *    Specifiy that we inherit our sample time from the driving block.
 */
static void mdlInitializeSampleTimes(SimStruct *S)
{
    const real_T *dt = mxGetPr(P_PAR_DT);

    ssSetSampleTime(S, 0, *dt);
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


#define MDL_JACOBIAN
/* Function: mdlJacobian =================================================
 * Abstract:
 *    J = d(dxc,xd,y)/d(xc,xd,u)
 *	=	x1	x2	u
 *	x1	1	0	dt
 *	x2	dt	1	0.5*dt*dt
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
  real_T dt = ssGetSampleTime(S, 0);

  pr = ssGetJacobianPr(S);
  ir = ssGetJacobianIr(S);
  jc = ssGetJacobianJc(S);
  j = 0;
  idx = 0;
  /* first column */
  jc[j++] = idx;
  ir[idx] = 0; pr[idx] = 1.0; idx++;
  ir[idx] = 1; pr[idx] = dt; idx++;
  ir[idx] = 2; pr[idx] = 1.0; idx++;
  /* second column */
  jc[j++] = idx;
  ir[idx] = 1; pr[idx] = 1.0; idx++;
  ir[idx] = 3; pr[idx] = 1.0; idx++;
  /* third column */
  jc[j++] = idx;
  ir[idx] = 0; pr[idx] = dt; idx++;
  ir[idx] = 1; pr[idx] = 0.5*dt*dt; idx++;
  /* end marker */
  jc[j] = idx;
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


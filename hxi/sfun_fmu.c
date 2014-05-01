/**
 * @file sfun_fmu.c
 *   Wrapper S-function for a Functional Model Unit (FMU)
 *   for model exchange according to FMI 2.0.
 *   This wrapper works together with fmi.tcl for the treatment of zip files
 *   and xml model descriptions.
 *   @see https://www.fmi-standard.org
 *
 * rf, 02/22/2014
 *
 */

/*
    Copyright (C) 1994--2014  Ruediger Franke

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

#define S_FUNCTION_NAME sfun_fmu
#define S_FUNCTION_LEVEL 2

#undef WITH_LOGGING

#include <string.h>
#include <stdarg.h>

/**
 * @name Dynamic loading
 * @{
 */

#if defined(_MSC_VER) || defined(__MINGW32__)
# include <windows.h>
# define DLHANDLE HMODULE
# define DLOPEN(filename) LoadLibrary(filename)
# define DLSYM GetProcAddress
# define DLCLOSE FreeLibrary
  static const char *DLERROR() {
    return "LoadLibrary failed for FMU";
  }
#else
# include <dlfcn.h>
# define DLHANDLE void*
# define DLOPEN(filename) dlopen(filename, RTLD_LAZY)
# define DLSYM dlsym
# define DLCLOSE dlclose
# define DLERROR dlerror
#endif

/** retrieve a function pointer from the loaded model, including error check */
#define INIT_FUNCTION(NAME) \
  m->NAME = (NAME ## TYPE*)DLSYM(m->handle, #NAME); \
  if (m->NAME == NULL) { \
    ssSetErrorStatus(S, "missing " #NAME); \
    return; \
  }

/**
 * @}
 */

/*#include <iftcl/If.h>*/
#include <tcl.h>
extern Tcl_Interp *If_Interp(); /**< Tcl interpreter used by the interface */

#include "fmiTypesPlatform.h"
#include "fmiFunctionTypes.h"

#include "simstruc.h"

/** prototype for initialization function exported by Hxi S-function */
typedef void (Hxi_SimStruct_init_t)(SimStruct *S);

/**
 * @name Additional FMI type definitions
 * @{
 */

static const char* BaseTypeNames[] = {"Real", "Boolean", "Integer", "String"};
#define NUM_BASETYPES (sizeof(BaseTypeNames)/sizeof(BaseTypeNames[0]))

/** Base types of FMI */
typedef enum {
  FMI_REAL,
  FMI_BOOLEAN,
  FMI_INTEGER,
  FMI_STRING
} BaseTypes;

/** Vector of FMU variables with different base types */
typedef struct {
  int nv;			/**< overall number of variables */
  BaseTypes *btv;		/**< base type per variable */
  int n[NUM_BASETYPES];		/**< number of variables per base type */
  fmiValueReference *vr[NUM_BASETYPES];	/**< value references per base type */
  fmiReal *r;			/**< Real values */
  fmiBoolean *b;		/**< Boolean values */
  fmiInteger *i;		/**< Integer values */
  fmiString *s;			/**< String values */
} Hxi_ModelVariables;

static const char* statusStrings[] = {
  "Message",
  "Warning",
  "Discard",
  "Error",
  "Fatal"
};

/**
 * @}
 */

/**
 * Instance data of a loaded FMU
 */
typedef struct {
  Tcl_Interp 		*interp; /*< used for model description and logging */

  char 			*fmuName;
  char 			*guid;
  char 			*generationTool;
  char 			*resourcesURI;
  char 			*messageStr;
  int 			messageStrLen;

  int 			np;
  int 			nxc;
  int 			nu;
  int 			ny;
  Hxi_ModelVariables 	*p;
  Hxi_ModelVariables 	*u;
  Hxi_ModelVariables 	*y;

  DLHANDLE 		handle;
  fmiComponent 		fmu;
  fmiCallbackFunctions  fmiCallbackFunctions;
  fmiBoolean 		initPending;
  /**< delay initialization until inputs are available */

  /**
   * @name FMI functions
   * @{
   */
  fmiInstantiateTYPE             	*fmiInstantiate;
  fmiEnterInitializationModeTYPE 	*fmiEnterInitializationMode;
  fmiExitInitializationModeTYPE  	*fmiExitInitializationMode;
  fmiNewDiscreteStatesTYPE       	*fmiNewDiscreteStates;
  fmiEnterContinuousTimeModeTYPE 	*fmiEnterContinuousTimeMode;
  fmiCompletedIntegratorStepTYPE 	*fmiCompletedIntegratorStep;
  fmiEnterEventModeTYPE          	*fmiEnterEventMode;
  fmiResetTYPE               		*fmiReset;
  fmiTerminateTYPE               	*fmiTerminate;
  fmiFreeInstanceTYPE            	*fmiFreeInstance;
  fmiSetTimeTYPE                 	*fmiSetTime;
  fmiSetContinuousStatesTYPE     	*fmiSetContinuousStates;
  fmiSetRealTYPE                 	*fmiSetReal;
  fmiSetBooleanTYPE              	*fmiSetBoolean;
  fmiSetIntegerTYPE              	*fmiSetInteger;
  fmiSetStringTYPE               	*fmiSetString;
  fmiGetContinuousStatesTYPE     	*fmiGetContinuousStates;
  fmiGetDerivativesTYPE          	*fmiGetDerivatives;
  fmiGetRealTYPE                 	*fmiGetReal;
  fmiGetBooleanTYPE              	*fmiGetBoolean;
  fmiGetIntegerTYPE              	*fmiGetInteger;
  fmiGetStringTYPE               	*fmiGetString;
  /**
   * @}
   */
} Hxi_ModelData;

/**
 * @name FMI callback functions
 * @{
 */

/** Log messages from FMU to Tcl stdout*/
static void callbackLogger(fmiComponentEnvironment ce, fmiString instanceName,
                           fmiStatus status, fmiString category,
                           fmiString message, ...)
{
  Hxi_ModelData *m = (Hxi_ModelData *)ce;
  va_list args;
  int n;

#ifndef WITH_LOGGING
  if (status == fmiOK)
    return;
#endif

  va_start(args, message);
  n = vsnprintf(m->messageStr, m->messageStrLen, message, args);
  /* enlarge m->messageStr up to max 16 KB */
  if (n >= m->messageStrLen && n < 16*1024) {
    free(m->messageStr);
    m->messageStr = malloc(n+1);
    m->messageStrLen = n+1;
    vsnprintf(m->messageStr, m->messageStrLen, message, args);
  }
  va_end(args);

  Tcl_VarEval(m->interp, "puts -nonewline {",  m->fmuName,
	      " ", statusStrings[status], "}", NULL);
  if (strlen(category) > 0 && strcmp(category, statusStrings[status]) != 0)
    Tcl_VarEval(m->interp, "puts -nonewline { ", category, ": }", NULL);
  else
    Tcl_Eval(m->interp, "puts -nonewline {: }");
  Tcl_VarEval(m->interp, "puts [::fmi::mapNames {", m->fmuName,
	      "} {", m->messageStr, "}]", NULL);
  Tcl_Eval(m->interp, "update idletasks");
}

/** Allocate memory for FMU */
static void* callbackAllocateMemory(size_t nobj, size_t size)
{
  return malloc(nobj*size);
}

/** Free allocated memory for FMU */
static void callbackFreeMemory(void* obj)
{
  free(obj);
}

/**
 * @}
 */

/**
 * @name Treatment of Hxi_ModelVariables
 * @{
 */

/** Free data structure used for Hxi_ModelVariables */
static void freeVariables(Hxi_ModelVariables *vars)
{
  int j;

  if (vars == NULL)
    return;

  if (vars->n[FMI_STRING] > 0)
    free((void *)(vars->s));
  if (vars->n[FMI_INTEGER] > 0)
    free(vars->i);
  if (vars->n[FMI_BOOLEAN] > 0)
    free(vars->b);
  if (vars->n[FMI_REAL] > 0)
    free(vars->r);
  for (j = NUM_BASETYPES - 1; j >= 0; j--) {
    if (vars->n[j] > 0)
      free(vars->vr[j]);
  }
  free(vars->btv);
  free(vars);
}

/** Retrieve Hxi_ModelVariables for one category from model description */
static Hxi_ModelVariables *getVariables(Hxi_ModelData *m, const char *category)
{
  int i, j, n, ref;
  char index[4*sizeof(fmiValueReference)]; /* upper bound for string length */
  Tcl_Obj *objPtr;
  const char *baseTypeName;
  int idx[NUM_BASETYPES];
  Hxi_ModelVariables *vars;

  /* obtain number of variables */
  if (Tcl_VarEval(m->interp, "llength ${::fmu::", m->fmuName,
		  "::", category, "References}", NULL) != TCL_OK
      || (objPtr = Tcl_GetObjResult(m->interp)) == NULL
      || Tcl_GetIntFromObj(m->interp, objPtr, &n) != TCL_OK) {
    return NULL;
  }

  /* allocate Hxi_ModelVariables data structure and obtain base types */
  vars = (Hxi_ModelVariables *)malloc(sizeof(Hxi_ModelVariables));
  vars->nv = n;
  vars->btv = (BaseTypes *)malloc(n*sizeof(BaseTypes));
  vars->r = NULL;
  vars->b = NULL;
  vars->i = NULL;
  vars->s = NULL;
  for (j = 0; j < NUM_BASETYPES; j++) {
    vars->vr[j] = NULL;
    vars->n[j] = 0;
  }
  for (i = 0; i < n; i++) {
    sprintf(index, "%d", i);
    if (Tcl_VarEval(m->interp, "lindex ${::fmu::", m->fmuName,
		    "::", category, "BaseTypes} ", index, NULL) != TCL_OK
	|| (baseTypeName = Tcl_GetStringResult(m->interp)) == NULL) {
      freeVariables(vars);
      return NULL;
    }
    for (j = 0; j < NUM_BASETYPES; j++) {
      if (baseTypeName[0] == BaseTypeNames[j][0]
	  && strcmp(baseTypeName, BaseTypeNames[j]) == 0) {
	  vars->btv[i] = j;
	  vars->n[j] ++;
	  break;
	}
    }
    if (j == NUM_BASETYPES) {
      /* unknown base type */
      freeVariables(vars);
      return NULL;
    }
  }

  /* allocate vectors per base type and fill in references */
  for (j = 0; j < NUM_BASETYPES; j++) {
    if (vars->n[j] > 0)
      vars->vr[j] =
	(fmiValueReference *)malloc((vars->n[j])*sizeof(fmiValueReference));
    idx[j] = 0;
  }
  if (vars->n[FMI_REAL] > 0)
    vars->r = (fmiReal *)malloc((vars->n[FMI_REAL])*sizeof(fmiReal));
  if (vars->n[FMI_BOOLEAN] > 0)
    vars->b = (fmiBoolean *)malloc((vars->n[FMI_BOOLEAN])*sizeof(fmiBoolean));
  if (vars->n[FMI_INTEGER] > 0)
    vars->i = (fmiInteger *)malloc((vars->n[FMI_INTEGER])*sizeof(fmiInteger));
  if (vars->n[FMI_STRING] > 0)
    vars->s = (fmiString *)malloc((vars->n[FMI_STRING])*sizeof(fmiString));
  for (i = 0; i < n; i++) {
    sprintf(index, "%d", i);
    if (Tcl_VarEval(m->interp, "lindex ${::fmu::", m->fmuName,
		    "::", category, "References} ", index, NULL) != TCL_OK
	|| (objPtr = Tcl_GetObjResult(m->interp)) == NULL
	|| Tcl_GetIntFromObj(m->interp, objPtr, &ref) != TCL_OK) {
      freeVariables(vars);
      return NULL;
    }
    j = vars->btv[i];
    if (idx[j] >= vars->n[j]) {
      /* inconsistent References vs. BaseTypes */
      freeVariables(vars);
      return NULL;
    }			
    vars->vr[j][idx[j]++] = ref;
  }

  return vars;
}

/** Write variable values to FMU */
static fmiStatus setValues(Hxi_ModelData *m, const Hxi_ModelVariables *vars)
{
  fmiStatus ret;
  if (vars->n[FMI_REAL] > 0) {
    ret = (*m->fmiSetReal)(m->fmu, vars->vr[FMI_REAL],
			   vars->n[FMI_REAL], vars->r);
    if (ret != fmiOK)
      return ret;
  }
  if (vars->n[FMI_BOOLEAN] > 0) {
    ret = (*m->fmiSetBoolean)(m->fmu, vars->vr[FMI_BOOLEAN],
			      vars->n[FMI_BOOLEAN], vars->b);
    if (ret != fmiOK)
      return ret;
  }
  if (vars->n[FMI_INTEGER] > 0) {
    ret = (*m->fmiSetInteger)(m->fmu, vars->vr[FMI_INTEGER],
			      vars->n[FMI_INTEGER], vars->i);
    if (ret != fmiOK)
      return ret;
  }
  if (vars->n[FMI_STRING] > 0) {
    ret = (*m->fmiSetString)(m->fmu, vars->vr[FMI_STRING],
			     vars->n[FMI_STRING], vars->s);
    if (ret != fmiOK)
      return ret;
  }
  return fmiOK;
}

/** Read variable values from FMU */
static fmiStatus getValues(Hxi_ModelData *m, Hxi_ModelVariables *vars)
{
  fmiStatus ret;
  if (vars->n[FMI_REAL] > 0) {
    ret = (*m->fmiGetReal)(m->fmu, vars->vr[FMI_REAL],
			   vars->n[FMI_REAL], vars->r);
    if (ret != fmiOK)
      return ret;
  }
  if (vars->n[FMI_BOOLEAN] > 0) {
    ret = (*m->fmiGetBoolean)(m->fmu, vars->vr[FMI_BOOLEAN],
			      vars->n[FMI_BOOLEAN], vars->b);
    if (ret != fmiOK)
      return ret;
  }
  if (vars->n[FMI_INTEGER] > 0) {
    ret = (*m->fmiGetInteger)(m->fmu, vars->vr[FMI_INTEGER],
			      vars->n[FMI_INTEGER], vars->i);
    if (ret != fmiOK)
      return ret;
  }
  return fmiOK;
}

/**
 * @}
 */

/**
 * @name S-function methods
 * @{
 */

#define MDL_CHECK_PARAMETERS
/**
 *  Validate parameters to verify they are okay.
 */
static void mdlCheckParameters(SimStruct *S)
{
  /* ToDo: should check parameter types */
}


/**
 *  Load FMU and initialize sizes of S-function
 *  (number of inputs, outputs, states, etc.).
 */
static void mdlInitializeSizes(SimStruct *S)
{
  Tcl_Interp *interp;
  Tcl_Obj *objPtr;
  char *fmuName;
  int doesDirectoryExist;
  int doesDescriptionExist;
  Hxi_ModelData *m;
  DLHANDLE handle;
  Hxi_SimStruct_init_t *Hxi_SimStruct_init_p;

  /*
   * Get Tcl interpreter
   */
  if ((interp = If_Interp()) == NULL
      || Tcl_InitStubs(interp, "8.0", 0) == NULL) {
    ssSetErrorStatus(S, "can't init Tcl");
    return;
  }

#ifdef WITH_LOGGING
  Tcl_Eval(interp, "puts mdlInitializeSizes");
#endif

  /*
   * Get name of FMU from path
   */
  if (Tcl_VarEval(interp, "::fmi::getName {", ssGetPath(S), "}",
                  NULL) != TCL_OK) {
    ssSetErrorStatus(S, "can't get FMU name from path");
    return;
  }
  fmuName = strdup(Tcl_GetStringResult(interp));

  /*
   * Test if model has been extracted before
   */
  if (Tcl_VarEval(interp, "::fmi::testExtracted {", ssGetPath(S), "}",
                  NULL) != TCL_OK
      || (objPtr = Tcl_GetObjResult(interp)) == NULL
      || Tcl_GetBooleanFromObj(interp, objPtr, &doesDirectoryExist) != TCL_OK) {
    ssSetErrorStatus(S, Tcl_GetStringResult(interp));
    return;
  }

  /*
   * Lookup existing model data
   */
  if (Tcl_VarEval(interp, "::set {::fmu::", fmuName,
                  "::modelData}", NULL) == TCL_OK) {
    if (sscanf(Tcl_GetStringResult(interp), "0x%p", &m) != 1) {
      ssSetErrorStatus(S, "can't get previously stored model instance");
      return;
    }
    doesDescriptionExist = 1;
  }
  else {
    m = NULL;
    doesDescriptionExist = 0;
  }

  /*
   * Release old model instance and extract new binary code
   */
  if (!doesDirectoryExist || !doesDescriptionExist) {
    if (m != NULL) {
      /* free previously allocated model instance */
      if (m->fmu != NULL)
        (*m->fmiFreeInstance)(m->fmu);
      free(m->guid);
      free(m->generationTool);
      free(m->resourcesURI);
      freeVariables(m->y);
      freeVariables(m->u);
      freeVariables(m->p);
      free(m->messageStr);
      free(m->fmuName);
      DLCLOSE(m->handle);
      free(m);
      m = NULL;
    }
    if (!doesDirectoryExist) {
      /* extract model from new FMU */
      if (Tcl_VarEval(interp,
                      "::fmi::extractModel {", ssGetPath(S), "}",
                      NULL) != TCL_OK
          || strcmp(Tcl_GetStringResult(interp), "") == 0) {
        ssSetErrorStatus(S, "can't extract model from FMU");
        return;
      }
    }
    /* read model description */
    if (Tcl_VarEval(interp, "::fmi::readModelDescription {", ssGetPath(S), "}",
                    NULL) != TCL_OK) {
      ssSetErrorStatus(S, "can't read model description of FMU");
      return;
    }
  }

  /*
   * Load binary model code
   */
  if (Tcl_VarEval(interp, "::fmi::getBinaryPath {", ssGetPath(S), "}",
		  NULL) != TCL_OK) {
    ssSetErrorStatus(S, "can't get path of FMU binary");
    return;
  }
  handle = DLOPEN(Tcl_GetStringResult(interp));
  if (!handle) {
    ssSetErrorStatus(S, DLERROR());
    return;
  }

  /*
   * Check for Hxi S-function, which is supported alternatively to FMI
   */
  Hxi_SimStruct_init_p =
    (Hxi_SimStruct_init_t*)DLSYM(handle, "Hxi_SimStruct_init");
  if (Hxi_SimStruct_init_p) {
    /* store null pointer to model data */
    if (Tcl_VarEval(interp, "::set {::fmu::", fmuName,
                    "::modelData} 0x00", NULL) != TCL_OK) {
      ssSetErrorStatus(S, "can't store model instance");
      return;
    }
    /* clear own references */
    free(fmuName);
    ssSetmdlCheckParameters(S, NULL);
    ssSetmdlInitializeSizes(S, NULL);
    ssSetmdlInitializeSampleTimes(S, NULL);
    ssSetmdlStart(S, NULL);
    ssSetmdlInitializeConditions(S, NULL);
    ssSetmdlUpdate(S, NULL);
    ssSetmdlDerivatives(S, NULL);
    ssSetmdlJacobian(S, NULL);
    ssSetmdlOutputs(S, NULL);
    ssSetmdlTerminate(S, NULL);
    /* redirect to Hxi S-function */
    (*Hxi_SimStruct_init_p)(S);
    Hxi_mdlInitializeSizes(S);
    return;
  }

  /*
   * create new model instance
   */
  if (m == NULL) {
    m = (Hxi_ModelData *)malloc(sizeof(Hxi_ModelData));
    memset(m, sizeof(Hxi_ModelData), 0);

    /*
     * Initialize FMU environment
     */
    m->interp = interp;
    m->fmuName = fmuName;
    m->messageStr = malloc(256);
    m->messageStrLen = 256;
    /* store pointer to model data for repeated calls */
    sprintf(m->messageStr, "0x%p", m);
    if (Tcl_VarEval(interp, "::set {::fmu::", fmuName,
                    "::modelData} ", m->messageStr, NULL) != TCL_OK) {
      ssSetErrorStatus(S, "can't store model instance");
      return;
    }
    /* get guid */
    if (Tcl_VarEval(interp, "::set {::fmu::", fmuName,
                    "::attributes(guid)}", NULL) != TCL_OK) {
      ssSetErrorStatus(S, "can't get guid of FMU");
      return;
    }
    m->guid = strdup(Tcl_GetStringResult(interp));
    /* get generation tool */
    if (Tcl_VarEval(interp, "::set {::fmu::", fmuName,
                    "::attributes(generationTool)}", NULL) != TCL_OK) {
      ssSetErrorStatus(S, "can't get generationTool of FMU");
      return;
    }
    m->generationTool = strdup(Tcl_GetStringResult(interp));
    /* get location of resources directory */
    if (Tcl_VarEval(interp, "::fmi::getResourcesURI {", ssGetPath(S), "}",
                    NULL) != TCL_OK) {
      ssSetErrorStatus(S, "can't get location of FMU resources");
      return;
    }
    m->resourcesURI = strdup(Tcl_GetStringResult(interp));

    /* 
     * Initialize parameters
     */
    if (Tcl_VarEval(interp, "llength ${::fmu::", fmuName,
                    "::parameterReferences}", NULL) != TCL_OK
        || (objPtr = Tcl_GetObjResult(interp)) == NULL
        || Tcl_GetIntFromObj(interp, objPtr, &m->np) != TCL_OK) {
      ssSetErrorStatus(S, Tcl_GetStringResult(interp));
      return;
    }
    if (!(m->p = getVariables(m, "parameter"))) {
      ssSetErrorStatus(S, "can't get parameters");
      return;
    }
    if (m->p->nv != m->np) {
      ssSetErrorStatus(S, "wrong parameter count");
      return;
    }

    /*
     * Initialize states, inputs and outputs
     */
    if (Tcl_VarEval(interp, "llength ${::fmu::", fmuName,
                    "::stateReferences}", NULL) != TCL_OK
        || (objPtr = Tcl_GetObjResult(interp)) == NULL
        || Tcl_GetIntFromObj(interp, objPtr, &m->nxc) != TCL_OK) {
      ssSetErrorStatus(S, "can't get number of states");
      return;
    }
    if (!(m->u = getVariables(m, "input"))) {
      ssSetErrorStatus(S, "can't get inputs");
      return;
    }
    m->nu = m->u->nv - m->u->n[FMI_STRING]; /* String inputs not supported */
    if (!(m->y = getVariables(m, "output"))) {
      ssSetErrorStatus(S, "can't get outputs");
      return;
    }
    m->ny = m->y->nv - m->y->n[FMI_STRING]; /* String outputs not supported */

    /* clean-up interp for potential upcoming error messages */
    Tcl_ResetResult(interp);

    /*
     * Initialize interface to binary
     */
    m->handle = handle;
    m->fmu = NULL;

    m->fmiCallbackFunctions.logger = &callbackLogger;
    m->fmiCallbackFunctions.allocateMemory = &callbackAllocateMemory;
    m->fmiCallbackFunctions.freeMemory = &callbackFreeMemory;
    m->fmiCallbackFunctions.stepFinished = NULL;
    m->fmiCallbackFunctions.componentEnvironment = m;

    INIT_FUNCTION(fmiInstantiate);
    INIT_FUNCTION(fmiEnterInitializationMode);
    INIT_FUNCTION(fmiExitInitializationMode);
    INIT_FUNCTION(fmiNewDiscreteStates);
    INIT_FUNCTION(fmiEnterContinuousTimeMode);
    INIT_FUNCTION(fmiCompletedIntegratorStep);
    INIT_FUNCTION(fmiEnterEventMode);
    INIT_FUNCTION(fmiReset);
    INIT_FUNCTION(fmiTerminate);
    INIT_FUNCTION(fmiFreeInstance);
    INIT_FUNCTION(fmiSetTime);
    INIT_FUNCTION(fmiSetContinuousStates);
    INIT_FUNCTION(fmiSetReal);
    INIT_FUNCTION(fmiSetBoolean);
    INIT_FUNCTION(fmiSetInteger);
    INIT_FUNCTION(fmiSetString);
    INIT_FUNCTION(fmiGetContinuousStates);
    INIT_FUNCTION(fmiGetDerivatives);
    INIT_FUNCTION(fmiGetReal);
    INIT_FUNCTION(fmiGetBoolean);
    INIT_FUNCTION(fmiGetInteger);
    INIT_FUNCTION(fmiGetString);
  }
  else
    free(fmuName);

  m->initPending = fmiFalse;

  /*
   * Initialize S-function
   */
  ssSetModelName(S, m->fmuName);
  ssSetNumSFcnParams(S, m->np);
  if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)) {
    if (ssGetSFcnParamsCount(S) != 0)
      return; /* Parameter mismatch will be reported by solver */
  }
  mdlCheckParameters(S);
  if (ssGetErrorStatus(S) != NULL) {
    return;
  }

  ssSetNumContStates(S, m->nxc);
  ssSetNumDiscStates(S, 0);

  if (!ssSetNumInputPorts(S, 1)) return;
  ssSetInputPortWidth(S, 0, m->nu);
  ssSetInputPortDirectFeedThrough(S, 0, 1);

  if (!ssSetNumOutputPorts(S, 1)) return;
  ssSetOutputPortWidth(S, 0, m->ny);

  ssSetNumSampleTimes(S, 1);
  ssSetNumRWork(S, 0);
  ssSetNumIWork(S, 0);
  ssSetNumPWork(S, 1);
  ssGetPWork(S)[0] = m; /* store pointer to model data */
  ssSetNumModes(S, 0);
  ssSetNumNonsampledZCs(S, 0);

  /* No long jumps */
  ssSetOptions(S, SS_OPTION_EXCEPTION_FREE_CODE);
}


/**
 *  Specifiy sample times.
 *  So far only one continuous sample time is used
 *  because an FMU treats discrete sample times internally.
 */
static void mdlInitializeSampleTimes(SimStruct *S)
{
#ifdef WITH_LOGGING
  Hxi_ModelData *m = (Hxi_ModelData *)ssGetPWork(S)[0];
  Tcl_Eval(m->interp, "puts mdlInitializeSampleTimes");
#endif
  ssSetSampleTime(S, 0, CONTINUOUS_SAMPLE_TIME);
  ssSetOffsetTime(S, 0, 0.0);
}


#define MDL_INITIALIZE_CONDITIONS
/**
 *  Initialize parameters and continuous states.
 */
static void mdlInitializeConditions(SimStruct *S)
{
  Hxi_ModelData *m = (Hxi_ModelData *)ssGetPWork(S)[0];

#ifdef WITH_LOGGING
  Tcl_Eval(m->interp, "puts mdlInitializeConditions");
#endif

  if (m->fmu != NULL) {
    if ((*m->fmiReset)(m->fmu) != fmiOK) {
      /* Free model instance if reset didn't work */
      (*m->fmiFreeInstance)(m->fmu);
      m->fmu = NULL;
    }
  }

  if (m->fmu == NULL) {
    m->fmu = (*m->fmiInstantiate)(m->fmuName, fmiModelExchange,
                                  m->guid, m->resourcesURI,
                                  &m->fmiCallbackFunctions,
                                  fmiFalse, fmiFalse);
    if (m->fmu == NULL) {
      ssSetErrorStatus(S, "can't instantiate FMU");
      return;
    }
  }

  m->initPending = fmiTrue; /* wait until inputs are available */

  /* get start values for continuous states */
  if((*m->fmiGetContinuousStates)(m->fmu, ssGetContStates(S), m->nxc)
     != fmiOK) {
    ssSetErrorStatus(S, "can't get continuous states of FMU");
    return;
  }
}


/**
 *  Set states and inputs and obtain outputs y = f(t,p,x,u).
 */
static void mdlOutputs(SimStruct *S, int_T tid)
{
  Hxi_ModelData *m = (Hxi_ModelData *)ssGetPWork(S)[0];
  fmiEventInfo eventInfo;
  fmiBoolean enterEventMode;
  fmiBoolean terminateSimulation;
  int i, j;
  int idx[NUM_BASETYPES];
  InputRealPtrsType uPtrs;
  mxArray *param;
  double *vals;

#ifdef WITH_LOGGING
  Tcl_Eval(m->interp, "puts mdlOutputs");
#endif

  if((*m->fmiSetTime)(m->fmu, ssGetT(S)) != fmiOK) {
    ssSetErrorStatus(S, "can't set time of FMU");
    return;
  }
  /* ToDo: set tunable parameters introduced with FMI 2.0 */
  if (m->initPending) {
    for (j = 0; j < NUM_BASETYPES; j++)
      idx[j] = 0;
    for (i = 0; i < m->p->nv; i++) {
      param = ssGetSFcnParam(S, i);
      j = m->p->btv[i];
      switch (j) {
      case FMI_REAL:
	m->p->r[idx[j]++] = *mxGetPr(param);
	break;
      case FMI_BOOLEAN:
	m->p->b[idx[j]++] = *mxGetPr(param) == 0.0? fmiFalse: fmiTrue;
	break;
      case FMI_INTEGER:
	m->p->i[idx[j]++] = (fmiInteger)(*mxGetPr(param));
	break;
      case FMI_STRING:
	m->p->s[idx[j]++] = mxArrayToString(param);
      default:
	break;
      }
    }
    if (setValues(m, m->p) != fmiOK)
      ssSetErrorStatus(S, "can't set parameters of FMU");
    idx[FMI_STRING] = 0;
    for (i = 0; i < m->p->nv; i++) {
      if (m->p->btv[i] == FMI_STRING)
	mxFree((void *)m->p->s[idx[FMI_STRING]++]);
    }
    if (ssGetErrorStatus(S))
      return;
  }

  /* set states */
  if(m->nxc > 0 &&
     (*m->fmiSetContinuousStates)(m->fmu, ssGetContStates(S), m->nxc)
     != fmiOK) {
    ssSetErrorStatus(S, "can't set continuous states of FMU");
    return;
  }

  /* set inputs */
  uPtrs = ssGetInputPortRealSignalPtrs(S, 0);
  for (j = 0; j < NUM_BASETYPES; j++)
    idx[j] = 0;
  for (i = 0; i < m->u->nv; i++) {
    j = m->u->btv[i];
    switch (j) {
    case FMI_REAL:
      m->u->r[idx[j]++] = *uPtrs[i];
      break;
    case FMI_BOOLEAN:
      m->u->b[idx[j]++] = *uPtrs[i] == 0.0? fmiFalse: fmiTrue;
      break;
    case FMI_INTEGER:
      m->u->i[idx[j]++] = (fmiInteger)(*uPtrs[i]);
      break;
    default:
      break;
    }
  }
  if (setValues(m, m->u) != fmiOK) {
    ssSetErrorStatus(S, "can't set inputs of FMU");
    return;
  }

  /* prepare event processing */
  if (m->initPending) {
    /* initialize FMU */
    m->initPending = fmiFalse;
    if ((*m->fmiEnterInitializationMode)(m->fmu) != fmiOK
        || (*m->fmiExitInitializationMode)(m->fmu) != fmiOK) {
      ssSetErrorStatus(S, "can't initialize FMU");
      return;
    }
    enterEventMode = fmiTrue;
  }
  else if (ssIsMinorTimeStep(S))
    /* no event processing in minor time steps */
    enterEventMode = fmiFalse;
  else {
    if (m->nxc > 0) {
      /* complete integrator step and optionally start event update */
      if ((*m->fmiCompletedIntegratorStep)(m->fmu, fmiTrue, &enterEventMode,
                                           &terminateSimulation) != fmiOK) {
        ssSetErrorStatus(S, "can't complete integrator step of FMU");
        return;
      }
      if (terminateSimulation) {
        ssSetErrorStatus(S, "terminate after integrator step of FMU");
        return;
      }
    }
    else
      /* always enter event mode if there are no continuous states */
      enterEventMode = fmiTrue;
    if (enterEventMode)
      (*m->fmiEnterEventMode)(m->fmu);
  }

  /* process events */
  if (enterEventMode) {
    eventInfo.newDiscreteStatesNeeded = fmiTrue;
    while (eventInfo.newDiscreteStatesNeeded) {
      (*m->fmiNewDiscreteStates)(m->fmu, &eventInfo);
      if (eventInfo.terminateSimulation) {
        ssSetErrorStatus(S, "terminate in event mode of FMU");
        return;
      }
    }
    if (m->nxc > 0)
      (*m->fmiEnterContinuousTimeMode)(m->fmu);
  }

  /* get outputs */
  vals = ssGetOutputPortRealSignal(S, 0);
  if (getValues(m, m->y) != fmiOK) {
    ssSetErrorStatus(S, "can't get outputs of FMU");
    return;
  }
  for (j = 0; j < NUM_BASETYPES; j++)
    idx[j] = 0;
  for (i = 0; i < m->y->nv; i++) {
    j = m->y->btv[i];
    switch (j) {
    case FMI_REAL:
      vals[i] = m->y->r[idx[j]++];
      break;
    case FMI_BOOLEAN:
      vals[i] = m->y->b[idx[j]++]? 1.0: 0.0;
      break;
    case FMI_INTEGER:
      vals[i] = m->y->i[idx[j]++];
      break;
    default:
      break;
    }
  }
}


#define MDL_UPDATE
/**
 *  Update discrete states x(k+1) = f(x(k),u(k)).
 *  Note: mdlOutputs must have been called before to set states and inputs 
 */
static void mdlUpdate(SimStruct *S, int_T tid)
{
  Hxi_ModelData *m = (Hxi_ModelData *)ssGetPWork(S)[0];

#ifdef WITH_LOGGING
  Tcl_Eval(m->interp, "puts mdlUpdate");
#endif

  if((*m->fmiGetContinuousStates)(m->fmu, ssGetContStates(S), m->nxc)
     != fmiOK) {
    ssSetErrorStatus(S, "can't get continuous states of FMU");
    return;
  }
}


#define MDL_DERIVATIVES
/**
 *  Obtain derivatives of continuous states dx = f(t,p,x,u).
 *  Note: mdlOutputs must have been called before to set states and inputs
 */
static void mdlDerivatives(SimStruct *S)
{
  Hxi_ModelData *m = (Hxi_ModelData *)ssGetPWork(S)[0];
#ifdef WITH_LOGGING
  Tcl_Eval(m->interp, "puts mdlDerivatives");
#endif

  if((*m->fmiGetDerivatives)(m->fmu, ssGetdX(S), m->nxc) != fmiOK) {
    ssSetErrorStatus(S, "can't get derivatives of FMU");
    return;
  }
}


/**
 *  Free model instance.
 */
static void mdlTerminate(SimStruct *S)
{
  Hxi_ModelData *m = (Hxi_ModelData *)ssGetPWork(S)[0];
#ifdef WITH_LOGGING
  Tcl_Eval(m->interp, "puts mdlTerminate");
#endif
  if((*m->fmiTerminate)(m->fmu) != fmiOK) {
    ssSetErrorStatus(S, "can't terminate FMU");
    return;
  }
}

/**
 * @}
 */

#ifdef  MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "simulink.c"      /* MEX-file interface mechanism */
#else
#include "cg_sfun.h"       /* Code generation registration function */
#endif


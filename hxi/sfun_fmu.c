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
#define INIT_FUNCTION(S, m, ID, NAME) \
  if (ID == 2) { \
    INIT_FUNCTION2(S, m, fmi2, NAME); \
  } \
  else { \
    INIT_FUNCTION2(S, m, fmi, NAME); \
  }

#define INIT_FUNCTION2(S, m, PREFIX, NAME) \
  m->fmi2 ## NAME = (fmi2 ## NAME ## TYPE*)DLSYM(m->handle, #PREFIX #NAME); \
  if (m->fmi2 ## NAME == NULL) { \
    ssSetErrorStatus(S, "missing " #PREFIX #NAME); \
    return; \
  }

#define GET_MODELDATA(S, m) \
  if (ssGetNumPWork(S) < 1) { \
    ssSetErrorStatus(S, "can't get model data"); \
    return; \
  } \
  m = (Hxi_ModelData *)ssGetPWork(S)[0]

/**
 * @}
 */

#include <iftcl/If.h>

#include "fmi2TypesPlatform.h"
#include "fmi2FunctionTypes.h"

#include "simstruc.h"

/** prototype for initialization function exported by Hxi S-function */
typedef void (Hxi_SimStruct_init_t)(SimStruct *S);

/**
 * @name Additional FMI type definitions
 * @{
 */

#define ARRAYSIZE(array) (sizeof(array)/sizeof(array[0]))
#define NUM_BASETYPES ARRAYSIZE(BaseTypeNames)
static const char* BaseTypeNames[] = {"Real", "Boolean", "Integer", "String"};

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
  fmi2ValueReference *vr[NUM_BASETYPES]; /**< value references per base type */
  fmi2Real *r;			/**< Real values */
  fmi2Boolean *b;		/**< Boolean values */
  fmi2Integer *i;		/**< Integer values */
  fmi2String *s;		/**< String values */
} Hxi_ModelVariables;

static const char* statusStrings[] = {
  "Message",
  "Warning",
  "Discard",
  "Error",
  "Fatal",
  "Pending"
};

static const char* logCategories[] = {
  "logAll",
  "logEvents", "logNonlinearSystems", "logDynamicStateSelection",
  "logStatusWarning", "logStatusPending", "logSingularLinearSystems",
  "logStatusDiscard", "logStatusError", "logStatusFatal"
};
// offsets into logCategories per log level 0 .. 4
static int logOffsets[] = {10, 7, 4, 1, 0};

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
  fmi2Component 	fmu;
  fmi2CallbackFunctions fmi2CallbackFunctions;
  fmi2Boolean 		initPending;
  /**< delay initialization until inputs are available */

  /**
   * @name FMI functions
   * @{
   */
  fmi2InstantiateTYPE             	*fmi2Instantiate;
  fmi2SetDebugLoggingTYPE             	*fmi2SetDebugLogging;
  fmi2EnterInitializationModeTYPE 	*fmi2EnterInitializationMode;
  fmi2ExitInitializationModeTYPE  	*fmi2ExitInitializationMode;
  fmi2NewDiscreteStatesTYPE       	*fmi2NewDiscreteStates;
  fmi2EnterContinuousTimeModeTYPE 	*fmi2EnterContinuousTimeMode;
  fmi2CompletedIntegratorStepTYPE 	*fmi2CompletedIntegratorStep;
  fmi2EnterEventModeTYPE          	*fmi2EnterEventMode;
  fmi2ResetTYPE               		*fmi2Reset;
  fmi2TerminateTYPE               	*fmi2Terminate;
  fmi2FreeInstanceTYPE            	*fmi2FreeInstance;
  fmi2SetTimeTYPE                 	*fmi2SetTime;
  fmi2SetContinuousStatesTYPE     	*fmi2SetContinuousStates;
  fmi2SetRealTYPE                 	*fmi2SetReal;
  fmi2SetBooleanTYPE              	*fmi2SetBoolean;
  fmi2SetIntegerTYPE              	*fmi2SetInteger;
  fmi2SetStringTYPE               	*fmi2SetString;
  fmi2GetContinuousStatesTYPE     	*fmi2GetContinuousStates;
  fmi2GetDerivativesTYPE          	*fmi2GetDerivatives;
  fmi2GetRealTYPE                 	*fmi2GetReal;
  fmi2GetBooleanTYPE              	*fmi2GetBoolean;
  fmi2GetIntegerTYPE              	*fmi2GetInteger;
  fmi2GetStringTYPE               	*fmi2GetString;
  /**
   * @}
   */
} Hxi_ModelData;

/**
 * @name FMI callback functions
 * @{
 */

/** Log messages from FMU to Tcl stdout*/
static void callbackLogger(fmi2ComponentEnvironment ce, fmi2String instanceName,
                           fmi2Status status, fmi2String category,
                           fmi2String message, ...)
{
  Hxi_ModelData *m = (Hxi_ModelData *)ce;
  va_list args;
  int n;

  /* format varargs into message string */
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

  /* replace variable references with names */
  Tcl_VarEval(m->interp, "::fmi::mapNames {", m->fmuName,
	      "} {", m->messageStr, "}", NULL);

  /* do the log */
  If_Log(statusStrings[status], "%s %s: %s", m->fmuName,
	 strcmp(category, statusStrings[status]) != 0? category: "",
	 Tcl_GetStringResult(m->interp));
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
  char index[4*sizeof(fmi2ValueReference)]; /* upper bound for string length */
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
	(fmi2ValueReference *)malloc((vars->n[j])*sizeof(fmi2ValueReference));
    idx[j] = 0;
  }
  if (vars->n[FMI_REAL] > 0)
    vars->r = (fmi2Real *)malloc((vars->n[FMI_REAL])*sizeof(fmi2Real));
  if (vars->n[FMI_BOOLEAN] > 0)
    vars->b = (fmi2Boolean *)malloc((vars->n[FMI_BOOLEAN])*sizeof(fmi2Boolean));
  if (vars->n[FMI_INTEGER] > 0)
    vars->i = (fmi2Integer *)malloc((vars->n[FMI_INTEGER])*sizeof(fmi2Integer));
  if (vars->n[FMI_STRING] > 0)
    vars->s = (fmi2String *)malloc((vars->n[FMI_STRING])*sizeof(fmi2String));
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
static fmi2Status setValues(Hxi_ModelData *m, const Hxi_ModelVariables *vars)
{
  fmi2Status ret;
  if (vars->n[FMI_REAL] > 0) {
    ret = (*m->fmi2SetReal)(m->fmu, vars->vr[FMI_REAL],
			    vars->n[FMI_REAL], vars->r);
    if (ret != fmi2OK)
      return ret;
  }
  if (vars->n[FMI_BOOLEAN] > 0) {
    ret = (*m->fmi2SetBoolean)(m->fmu, vars->vr[FMI_BOOLEAN],
			       vars->n[FMI_BOOLEAN], vars->b);
    if (ret != fmi2OK)
      return ret;
  }
  if (vars->n[FMI_INTEGER] > 0) {
    ret = (*m->fmi2SetInteger)(m->fmu, vars->vr[FMI_INTEGER],
			       vars->n[FMI_INTEGER], vars->i);
    if (ret != fmi2OK)
      return ret;
  }
  if (vars->n[FMI_STRING] > 0) {
    ret = (*m->fmi2SetString)(m->fmu, vars->vr[FMI_STRING],
			      vars->n[FMI_STRING], vars->s);
    if (ret != fmi2OK)
      return ret;
  }
  return fmi2OK;
}

/** Read variable values from FMU */
static fmi2Status getValues(Hxi_ModelData *m, Hxi_ModelVariables *vars)
{
  fmi2Status ret;
  if (vars->n[FMI_REAL] > 0) {
    ret = (*m->fmi2GetReal)(m->fmu, vars->vr[FMI_REAL],
			    vars->n[FMI_REAL], vars->r);
    if (ret != fmi2OK)
      return ret;
  }
  if (vars->n[FMI_BOOLEAN] > 0) {
    ret = (*m->fmi2GetBoolean)(m->fmu, vars->vr[FMI_BOOLEAN],
			       vars->n[FMI_BOOLEAN], vars->b);
    if (ret != fmi2OK)
      return ret;
  }
  if (vars->n[FMI_INTEGER] > 0) {
    ret = (*m->fmi2GetInteger)(m->fmu, vars->vr[FMI_INTEGER],
			       vars->n[FMI_INTEGER], vars->i);
    if (ret != fmi2OK)
      return ret;
  }
  return fmi2OK;
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
  int version_id;

  /*
   * Get Tcl interpreter
   */
  if ((interp = If_Interp()) == NULL
      || Tcl_InitStubs(interp, "8.0", 0) == NULL) {
    ssSetErrorStatus(S, "can't init Tcl");
    return;
  }

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
        (*m->fmi2FreeInstance)(m->fmu);
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
                      NULL) != TCL_OK) {
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
    if (strcmp(Tcl_GetStringResult(interp), "2") < 0) {
      Tcl_ResetResult(interp);
      ssSetErrorStatus(S, "require FMI version 2.0");
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
    m = (Hxi_ModelData *)calloc(1, sizeof(Hxi_ModelData));

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

    m->fmi2CallbackFunctions.logger = &callbackLogger;
    m->fmi2CallbackFunctions.allocateMemory = &calloc;
    m->fmi2CallbackFunctions.freeMemory = &free;
    m->fmi2CallbackFunctions.stepFinished = NULL;
    m->fmi2CallbackFunctions.componentEnvironment = m;

    /* support FMI 2.0 with prefix "fmi2" and 2.0 RC1 with prefix "fmi" */
    version_id = 2;
    if (DLSYM(m->handle, "fmi2Instantiate") == NULL
	&& DLSYM(m->handle, "fmiInstantiate") != NULL)
      version_id = 1;

    INIT_FUNCTION(S, m, version_id, Instantiate);
    INIT_FUNCTION(S, m, version_id, SetDebugLogging);
    INIT_FUNCTION(S, m, version_id, EnterInitializationMode);
    INIT_FUNCTION(S, m, version_id, ExitInitializationMode);
    INIT_FUNCTION(S, m, version_id, NewDiscreteStates);
    INIT_FUNCTION(S, m, version_id, EnterContinuousTimeMode);
    INIT_FUNCTION(S, m, version_id, CompletedIntegratorStep);
    INIT_FUNCTION(S, m, version_id, EnterEventMode);
    INIT_FUNCTION(S, m, version_id, Reset);
    INIT_FUNCTION(S, m, version_id, Terminate);
    INIT_FUNCTION(S, m, version_id, FreeInstance);
    INIT_FUNCTION(S, m, version_id, SetTime);
    INIT_FUNCTION(S, m, version_id, SetContinuousStates);
    INIT_FUNCTION(S, m, version_id, SetReal);
    INIT_FUNCTION(S, m, version_id, SetBoolean);
    INIT_FUNCTION(S, m, version_id, SetInteger);
    INIT_FUNCTION(S, m, version_id, SetString);
    INIT_FUNCTION(S, m, version_id, GetContinuousStates);
    INIT_FUNCTION(S, m, version_id, GetDerivatives);
    INIT_FUNCTION(S, m, version_id, GetReal);
    INIT_FUNCTION(S, m, version_id, GetBoolean);
    INIT_FUNCTION(S, m, version_id, GetInteger);
    INIT_FUNCTION(S, m, version_id, GetString);
  }
  else
    free(fmuName);

  m->initPending = fmi2False;

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
  ssSetSampleTime(S, 0, CONTINUOUS_SAMPLE_TIME);
  ssSetOffsetTime(S, 0, 0.0);
}


#define MDL_INITIALIZE_CONDITIONS
/**
 *  Initialize parameters and continuous states.
 */
static void mdlInitializeConditions(SimStruct *S)
{
  Hxi_ModelData *m;
  int logging;

  GET_MODELDATA(S, m);

  if (m->fmu != NULL) {
    if ((*m->fmi2Reset)(m->fmu) != fmi2OK) {
      /* Free model instance if reset didn't work */
      (*m->fmi2FreeInstance)(m->fmu);
      m->fmu = NULL;
    }
  }

  if (m->fmu == NULL) {
    m->fmu = (*m->fmi2Instantiate)(m->fmuName, fmi2ModelExchange,
				   m->guid, m->resourcesURI,
				   &m->fmi2CallbackFunctions,
				   fmi2False, fmi2False);
    if (m->fmu == NULL) {
      ssSetErrorStatus(S, "can't instantiate FMU");
      return;
    }
    /*
     * Initialize logging
     * ToDo: should check categories with modelDescription.xml
     */
    if (If_GetInt("mdl_logging", &logging) != IF_OK) {
      ssSetErrorStatus(S, If_ResultString());
      return;
    }
    if (logging > 0) {
      if (logging > ARRAYSIZE(logOffsets) - 1)
	logging = ARRAYSIZE(logOffsets) - 1;
      (*m->fmi2SetDebugLogging)(m->fmu, fmi2True,
				ARRAYSIZE(logCategories) - logOffsets[logging],
				&logCategories[logOffsets[logging]]);
    }
  }

  m->initPending = fmi2True; /* wait until inputs are available */
}


/**
 *  Set states and inputs and obtain outputs y = f(t,p,x,u).
 */
static void mdlOutputs(SimStruct *S, int_T tid)
{
  Hxi_ModelData *m;
  fmi2EventInfo eventInfo;
  fmi2Boolean enterEventMode;
  fmi2Boolean terminateSimulation;
  int i, j;
  int idx[NUM_BASETYPES];
  InputRealPtrsType uPtrs;
  mxArray *param;
  double *vals;

  GET_MODELDATA(S, m);

  if ((*m->fmi2SetTime)(m->fmu, ssGetT(S)) != fmi2OK) {
    ssSetErrorStatus(S, "can't set time of FMU");
    return;
  }
  /* Set all parameter values at initialization */
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
	m->p->b[idx[j]++] = *mxGetPr(param) == 0.0? fmi2False: fmi2True;
	break;
      case FMI_INTEGER:
	m->p->i[idx[j]++] = (fmi2Integer)(*mxGetPr(param));
	break;
      case FMI_STRING:
	m->p->s[idx[j]++] = mxArrayToString(param);
      default:
	break;
      }
    }
    if (setValues(m, m->p) != fmi2OK)
      ssSetErrorStatus(S, "can't set parameters of FMU");
    idx[FMI_STRING] = 0;
    for (i = 0; i < m->p->nv; i++) {
      if (m->p->btv[i] == FMI_STRING)
	mxFree((void *)m->p->s[idx[FMI_STRING]++]);
    }
    if (ssGetErrorStatus(S))
      return;
  }
  /* Set changed values else, to support tunable parameters of FMI 2.0. */
  else {
    for (j = 0; j < NUM_BASETYPES; j++)
      idx[j] = 0;
    for (i = 0; i < m->p->nv; i++) {
      param = ssGetSFcnParam(S, i);
      j = m->p->btv[i];
      switch (j) {
      case FMI_REAL:
        if (m->p->r[idx[j]] != *mxGetPr(param)) {
          m->p->r[idx[j]] = *mxGetPr(param);
          if ((*m->fmi2SetReal)(m->fmu, m->p->vr[FMI_REAL] + idx[j],
                                1, m->p->r + idx[j]) != fmi2OK) {
            ssSetErrorStatus(S, "can't set Real parameter of FMU");
            return;
          }
        }
        idx[j]++;
	break;
      case FMI_BOOLEAN:
        if (m->p->b[idx[j]] != (*mxGetPr(param) == 0.0? fmi2False: fmi2True)) {
          m->p->b[idx[j]] = *mxGetPr(param) == 0.0? fmi2False: fmi2True;
          if ((*m->fmi2SetBoolean)(m->fmu, m->p->vr[FMI_BOOLEAN] + idx[j],
                                   1, m->p->b + idx[j]) != fmi2OK) {
            ssSetErrorStatus(S, "can't set Boolean parameter of FMU");
            return;
          }
        }
        idx[j]++;
	break;
      case FMI_INTEGER:
        if (m->p->i[idx[j]] != (fmi2Integer)(*mxGetPr(param))) {
          m->p->i[idx[j]] = (fmi2Integer)(*mxGetPr(param));
          if ((*m->fmi2SetInteger)(m->fmu, m->p->vr[FMI_INTEGER] + idx[j],
                                   1, m->p->i + idx[j]) != fmi2OK) {
            ssSetErrorStatus(S, "can't set Integer parameter of FMU");
            return;
          }
        }
        idx[j]++;
	break;
      case FMI_STRING:
        /* ToDo: support tunable string parameters */
      default:
	break;
      }
    }
  }

  /* set states */
  if (m->nxc > 0 &&
      (*m->fmi2SetContinuousStates)(m->fmu, ssGetContStates(S), m->nxc)
      != fmi2OK) {
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
      m->u->b[idx[j]++] = *uPtrs[i] == 0.0? fmi2False: fmi2True;
      break;
    case FMI_INTEGER:
      m->u->i[idx[j]++] = (fmi2Integer)(*uPtrs[i]);
      break;
    default:
      break;
    }
  }
  if (setValues(m, m->u) != fmi2OK) {
    ssSetErrorStatus(S, "can't set inputs of FMU");
    return;
  }

  /* prepare event processing */
  if (m->initPending) {
    /* initialize FMU */
    m->initPending = fmi2False;
    if ((*m->fmi2EnterInitializationMode)(m->fmu) != fmi2OK
        || (*m->fmi2ExitInitializationMode)(m->fmu) != fmi2OK) {
      ssSetErrorStatus(S, "can't initialize FMU");
      return;
    }
    enterEventMode = fmi2True;
  }
  else if (ssIsMinorTimeStep(S))
    /* no event processing in minor time steps */
    enterEventMode = fmi2False;
  else {
    if (m->nxc > 0) {
      /* complete integrator step and optionally start event update */
      if ((*m->fmi2CompletedIntegratorStep)(m->fmu, fmi2True, &enterEventMode,
                                            &terminateSimulation) != fmi2OK) {
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
      enterEventMode = fmi2True;
    if (enterEventMode)
      (*m->fmi2EnterEventMode)(m->fmu);
  }

  /* process events */
  if (enterEventMode) {
    eventInfo.newDiscreteStatesNeeded = fmi2True;
    while (eventInfo.newDiscreteStatesNeeded) {
      (*m->fmi2NewDiscreteStates)(m->fmu, &eventInfo);
      if (eventInfo.terminateSimulation) {
        ssSetErrorStatus(S, "terminate in event mode of FMU");
        return;
      }
    }
    if (m->nxc > 0)
      (*m->fmi2EnterContinuousTimeMode)(m->fmu);
  }

  /* get outputs */
  vals = ssGetOutputPortRealSignal(S, 0);
  if (getValues(m, m->y) != fmi2OK) {
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
  Hxi_ModelData *m;

  GET_MODELDATA(S, m);

  if (m->nxc > 0 &&
      (*m->fmi2GetContinuousStates)(m->fmu, ssGetContStates(S), m->nxc)
      != fmi2OK) {
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
  Hxi_ModelData *m;

  GET_MODELDATA(S, m);

  if (m->nxc > 0 &&
      (*m->fmi2GetDerivatives)(m->fmu, ssGetdX(S), m->nxc) != fmi2OK) {
    ssSetErrorStatus(S, "can't get derivatives of FMU");
    return;
  }
}


/**
 *  Free model instance.
 */
static void mdlTerminate(SimStruct *S)
{
  Hxi_ModelData *m;

  GET_MODELDATA(S, m);

  if (m->fmi2Terminate != NULL && (*m->fmi2Terminate)(m->fmu) != fmi2OK) {
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

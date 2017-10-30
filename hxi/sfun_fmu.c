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
    Copyright (C) 1994--2017  Ruediger Franke

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

/** signature for initialization function exported by Hxi S-function */
typedef void (Hxi_SimStruct_init_t)(SimStruct *S);

/**
 * @name Additional FMI type definitions
 * @{
 */

typedef fmi2Status fmi2SetClockTYPE(fmi2Component, const int fmiInteger[],
                                    size_t, const fmi2Boolean tick[],
                                    const fmi2Boolean* subactive);
typedef fmi2Status fmi2SetIntervalTYPE(fmi2Component, const int fmiInteger[],
                                       size_t, const fmi2Real[]);

#define ARRAY_SIZE(array) (sizeof(array)/sizeof(array[0]))
#define NUM_BASETYPES ARRAY_SIZE(BaseTypeNames)
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
  long 			fmuTime;
  char 			*guid;
  char 			*generationTool;
  char 			*resourcesURI;
  int   		providesDirectionalDerivative;
  char 			*messageStr;
  int 			messageStrLen;
  int 			logging;

  int 			np;
  int 			nxd;
  int 			nxc;
  int 			nu;
  int 			ny;
  int 			nz;
  int			nc;
  Hxi_ModelVariables 	*p;
  Hxi_ModelVariables	*pre_x;
  double		*pre_x_vals;
  Hxi_ModelVariables	*dx;
  int 			*dx_perm;
  Hxi_ModelVariables	*x;
  Hxi_ModelVariables 	*u;
  Hxi_ModelVariables 	*y;
  fmi2Real		*z;
  fmi2Real		*pre_z;
  fmi2Integer		*cidx;
  fmi2Boolean		*ctck;
  fmi2Boolean		*csub;
  fmi2Real		*civl;

  DLHANDLE 		handle;
  fmi2Component 	fmu;
  fmi2CallbackFunctions fmi2CallbackFunctions;
  fmi2Boolean 		initPending;
  /**< delay initialization until inputs are available */
  fmi2Boolean 		nextEventTimeDefined;
  fmi2Real 		nextEventTime;

  /**
   * @name FMI functions
   * @{
   */
  fmi2InstantiateTYPE             	*fmi2Instantiate;
  fmi2SetDebugLoggingTYPE             	*fmi2SetDebugLogging;
  fmi2SetupExperimentTYPE             	*fmi2SetupExperiment;
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
  fmi2SetClockTYPE               	*fmi2SetClock;
  fmi2SetIntervalTYPE               	*fmi2SetInterval;
  fmi2GetContinuousStatesTYPE     	*fmi2GetContinuousStates;
  fmi2GetDerivativesTYPE          	*fmi2GetDerivatives;
  fmi2GetEventIndicatorsTYPE       	*fmi2GetEventIndicators;
  fmi2GetRealTYPE                 	*fmi2GetReal;
  fmi2GetBooleanTYPE              	*fmi2GetBoolean;
  fmi2GetIntegerTYPE              	*fmi2GetInteger;
  fmi2GetStringTYPE               	*fmi2GetString;
  fmi2GetDirectionalDerivativeTYPE 	*fmi2GetDirectionalDerivative;
  fmi2ValueReference               	*vrUnknown;
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
  int hasCategory = category != NULL && category[0] != '\0'
    && strcmp(category, statusStrings[status]) != 0;

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
  If_Log(statusStrings[status], "%s%s%s: %s", m->fmuName,
         hasCategory? " ": "", hasCategory? category: "",
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
  long fmuTime;
  int doesDirectoryExist;
  int doesDescriptionExist;
  int doesInstanceExist;
  Hxi_ModelData *m;
  DLHANDLE handle;
  Hxi_SimStruct_init_t *Hxi_SimStruct_init_p;
  int i;
  char tnStr[20];

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
   * Test if XML model description has been parsed before
   */
  if (!doesDirectoryExist) {
    doesDescriptionExist = 0;
  }
  else if (Tcl_VarEval(interp, "info exists ::fmu::", fmuName, "::attributes",
                       NULL) != TCL_OK
           || (objPtr = Tcl_GetObjResult(interp)) == NULL
           || Tcl_GetBooleanFromObj(interp, objPtr, &doesDescriptionExist) != TCL_OK) {
    ssSetErrorStatus(S, Tcl_GetStringResult(interp));
    return;
  }

  /*
   * Lookup existing model instance and its mtime
   */
  m = NULL;
  doesInstanceExist = 0;
  fmuTime = 0;
  sprintf(tnStr, "%d", ssGetOptions(S));
  if (Tcl_VarEval(interp, "::set {::fmu::", fmuName,
                  "::modelData", tnStr, "}", NULL) == TCL_OK) {
    if (sscanf(Tcl_GetStringResult(interp), "0x%p", &m) != 1) {
      ssSetErrorStatus(S, "can't get previously stored model instance");
      return;
    }
    if (m == NULL) {
      doesInstanceExist = 1; /* Hxi S-function */
    }
    else {
      /* get modification time */
      if (Tcl_VarEval(interp, "file mtime {", ssGetPath(S), "}" , NULL) != TCL_OK
          || (objPtr = Tcl_GetObjResult(m->interp)) == NULL
          || Tcl_GetLongFromObj(m->interp, objPtr, &fmuTime) != TCL_OK) {
        ssSetErrorStatus(S, "can't get mtime of FMU");
        return;
      }
      doesInstanceExist = (fmuTime == m->fmuTime);
    }
  }

  /*
   * Release old model instance and extract new binary code
   */
  if (!doesDirectoryExist || !doesDescriptionExist || !doesInstanceExist) {
    if (m != NULL) {
      /* free previously allocated model instance */
      if (m->fmu != NULL)
        (*m->fmi2FreeInstance)(m->fmu);
      free(m->guid);
      free(m->generationTool);
      free(m->resourcesURI);
      if (m->vrUnknown != NULL)
        free(m->vrUnknown);
      if (m->civl != NULL)
	free(m->civl);
      if (m->csub != NULL)
	free(m->csub);
      if (m->ctck != NULL)
	free(m->ctck);
      if (m->cidx != NULL)
	free(m->cidx);
      if (m->pre_z != NULL)
	free(m->pre_z);
      if (m->z != NULL)
	free(m->z);
      freeVariables(m->y);
      freeVariables(m->u);
      freeVariables(m->x);
      free(m->dx_perm);
      freeVariables(m->dx);
      freeVariables(m->pre_x);
      free(m->pre_x_vals);
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
    if (!doesDirectoryExist || !doesDescriptionExist) {
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
  }

  /*
   * Create new model instance if needed
   */
  if (m == NULL) {
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
                      "::modelData", tnStr, "} 0x00", NULL) != TCL_OK) {
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
     * Initialize FMU environment
     */
    m = (Hxi_ModelData *)calloc(1, sizeof(Hxi_ModelData));
    m->interp = interp;
    m->fmuName = fmuName;
    m->fmuTime = fmuTime;
    m->messageStr = malloc(256);
    m->messageStrLen = 256;
    /* store pointer to model data for repeated calls */
    sprintf(m->messageStr, "0x%p", m);
    if (Tcl_VarEval(interp, "::set {::fmu::", fmuName,
                    "::modelData", tnStr, "} ", m->messageStr, NULL) != TCL_OK) {
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
    /* get providesDirectionalDerivative attribute */
    if (Tcl_VarEval(interp, "::set {::fmu::", fmuName,
                    "::attributes(providesDirectionalDerivative)}", NULL) != TCL_OK
        || (objPtr = Tcl_GetObjResult(interp)) == NULL
        || Tcl_GetBooleanFromObj(interp, objPtr, &m->providesDirectionalDerivative) != TCL_OK) {
      m->providesDirectionalDerivative = 0;
    }

    /* 
     * Initialize parameters
     */
    if (!(m->p = getVariables(m, "parameter"))) {
      ssSetErrorStatus(S, "can't get parameters");
      return;
    }
    m->np = m->p->nv;

    /*
     * Initialize states, derivatives and previous states
     */
    if (!(m->pre_x = getVariables(m, "previous"))) {
      ssSetErrorStatus(S, "can't get previous states");
      return;
    }
    m->nxd = m->pre_x->nv;
    m->pre_x_vals = (fmi2Real *)calloc(m->nxd, sizeof(fmi2Real));
    if (!(m->dx = getVariables(m, "derivative"))) {
      ssSetErrorStatus(S, "can't get state derivatives");
      return;
    }
    if (m->dx->nv != m->dx->n[FMI_REAL]) {
      ssSetErrorStatus(S, "continuous states must be of type real");
      return;
    }
    m->nxc = m->dx->nv;
    if (!(m->x = getVariables(m, "state"))) {
      ssSetErrorStatus(S, "can't get states");
      return;
    }
    if (m->nxd + m->nxc != m->x->nv) {
      ssSetErrorStatus(S, "inconsistent number of states");
      return;
    }
    // hide continuous states because they have special access functions
    m->x->nv -= m->nxc;
    m->x->n[FMI_REAL] -= m->nxc;

    m->dx_perm = (int *)malloc(m->nxc * sizeof(int));
    for (i = 0; i < m->nxc; i++) {
      sprintf(m->messageStr, "%d", i);
      if (Tcl_VarEval(m->interp, "lindex ${::fmu::", m->fmuName,
                      "::derivativePermutation} ", m->messageStr, NULL) != TCL_OK
          || (objPtr = Tcl_GetObjResult(m->interp)) == NULL
          || Tcl_GetIntFromObj(m->interp, objPtr, &m->dx_perm[i]) != TCL_OK) {
        ssSetErrorStatus(S, "can't get state permutation of FMU");
        return;
      }
      if (m->dx_perm[i] < 0 || m->nxc <= m->dx_perm[i]) {
        ssSetErrorStatus(S, "got wrong state permutation of FMU");
        return;
      }
    }

    /*
     * Initialize inputs and outputs
     */
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

    /*
     * Initialize event indicators
     */
    if (Tcl_VarEval(interp, "::set {::fmu::", fmuName,
                    "::attributes(numberOfEventIndicators)}", NULL) != TCL_OK
        || (objPtr = Tcl_GetObjResult(interp)) == NULL
        || Tcl_GetIntFromObj(interp, objPtr, &m->nz) != TCL_OK) {
      ssSetErrorStatus(S, "can't get numberOfEventIndicators of FMU");
      return;
    }
    if (m->nz > 0) {
      m->z = (fmi2Real *)calloc(m->nz, sizeof(fmi2Real));
      m->pre_z = (fmi2Real *)calloc(m->nz, sizeof(fmi2Real));
    }

    /*
     * Initialize clocks
     */
    if (Tcl_VarEval(m->interp, "llength ${::fmu::", m->fmuName,
                    "::clockIntervals}", NULL) != TCL_OK
        || (objPtr = Tcl_GetObjResult(m->interp)) == NULL
        || Tcl_GetIntFromObj(m->interp, objPtr, &m->nc) != TCL_OK) {
      ssSetErrorStatus(S, "can't get number of clocks");
      return;
    }
    if (m->nc > 0) {
      m->cidx = (fmi2Integer *)calloc(m->nc, sizeof(fmi2Integer));
      m->ctck = (fmi2Boolean *)calloc(m->nc, sizeof(fmi2Boolean));
      m->csub = (fmi2Boolean *)calloc(m->nc, sizeof(fmi2Boolean));
      m->civl = (fmi2Real *)calloc(m->nc, sizeof(fmi2Real));
      for (i = 0; i < m->nc; i++) {
        m->cidx[i] = i + 1;
        m->ctck[i] = fmi2True;
        m->csub[i] = fmi2False;
      }
    }

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
    i = 2;
    if (DLSYM(m->handle, "fmi2Instantiate") == NULL
	&& DLSYM(m->handle, "fmiInstantiate") != NULL)
      i = 1;

    INIT_FUNCTION(S, m, i, Instantiate);
    INIT_FUNCTION(S, m, i, SetDebugLogging);
    INIT_FUNCTION(S, m, i, SetupExperiment);
    INIT_FUNCTION(S, m, i, EnterInitializationMode);
    INIT_FUNCTION(S, m, i, ExitInitializationMode);
    INIT_FUNCTION(S, m, i, NewDiscreteStates);
    INIT_FUNCTION(S, m, i, EnterContinuousTimeMode);
    INIT_FUNCTION(S, m, i, CompletedIntegratorStep);
    INIT_FUNCTION(S, m, i, EnterEventMode);
    INIT_FUNCTION(S, m, i, Reset);
    INIT_FUNCTION(S, m, i, Terminate);
    INIT_FUNCTION(S, m, i, FreeInstance);
    INIT_FUNCTION(S, m, i, SetTime);
    INIT_FUNCTION(S, m, i, SetContinuousStates);
    INIT_FUNCTION(S, m, i, SetReal);
    INIT_FUNCTION(S, m, i, SetBoolean);
    INIT_FUNCTION(S, m, i, SetInteger);
    INIT_FUNCTION(S, m, i, SetString);
    m->fmi2SetClock = (fmi2SetClockTYPE*)DLSYM(m->handle, "fmi2SetClock");
    m->fmi2SetInterval = (fmi2SetIntervalTYPE*)DLSYM(m->handle, "fmi2SetInterval");
    INIT_FUNCTION(S, m, i, GetContinuousStates);
    INIT_FUNCTION(S, m, i, GetDerivatives);
    INIT_FUNCTION(S, m, i, GetEventIndicators);
    INIT_FUNCTION(S, m, i, GetReal);
    INIT_FUNCTION(S, m, i, GetBoolean);
    INIT_FUNCTION(S, m, i, GetInteger);
    INIT_FUNCTION(S, m, i, GetString);
    if (m->providesDirectionalDerivative) {
      INIT_FUNCTION(S, m, i, GetDirectionalDerivative);
      m->vrUnknown = (fmi2ValueReference *)calloc(m->nxc + m->nxd + m->ny,
                                                  sizeof(fmi2ValueReference));
    }
  }
  else
    free(fmuName);

  m->initPending = fmi2True;

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

  ssSetNumDiscStates(S, m->nxd);
  ssSetNumContStates(S, m->nxc);

  if (!ssSetNumInputPorts(S, 1)) return;
  ssSetInputPortWidth(S, 0, m->nu);
  ssSetInputPortDirectFeedThrough(S, 0, 1);

  if (!ssSetNumOutputPorts(S, 1)) return;
  ssSetOutputPortWidth(S, 0, m->ny);

  ssSetNumSampleTimes(S, 1 + m->nc);
  ssSetNumRWork(S, 0);
  ssSetNumIWork(S, 0);
  ssSetNumPWork(S, 1);
  ssGetPWork(S)[0] = m; /* store pointer to model data */
  ssSetNumModes(S, 0);
  ssSetNumNonsampledZCs(S, 0);

  /* No long jumps */
  ssSetOptions(S, SS_OPTION_EXCEPTION_FREE_CODE);

  /* Setup Jacobian */
  if (Tcl_VarEval(m->interp, "::set ::fmu::", m->fmuName,
		  "::numberOfDependencies", NULL) != TCL_OK
      || (objPtr = Tcl_GetObjResult(m->interp)) == NULL
      || Tcl_GetIntFromObj(m->interp, objPtr, &i) != TCL_OK
      || i < 0) {
    ssSetErrorStatus(S, "can't get numberOfDependencies of FMU");
    return;
  }
  else
    ssSetJacobianNzMax(S, i);
}


/**
 *  Specifiy sample times.
 *  So far only one continuous sample time is used
 *  because an FMU treats discrete sample times internally.
 */
static void mdlInitializeSampleTimes(SimStruct *S)
{
  int i;
  Hxi_ModelData *m;

  GET_MODELDATA(S, m);

  ssSetSampleTime(S, 0, CONTINUOUS_SAMPLE_TIME);
  ssSetOffsetTime(S, 0, 0.0);
  for (i = 1; i <= m->nc; i++) {
    ssSetSampleTime(S, i, 1.0); // TODO: obtain actual interval
    ssSetOffsetTime(S, i, 0.0);
  }
}


#define MDL_INITIALIZE_CONDITIONS
/**
 *  Initialize the FMU.
 */
static void mdlInitializeConditions(SimStruct *S)
{
  Hxi_ModelData *m;
  int ret, logging;

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
    m->logging = 0;
  }

  /*
   * Initialize logging
   * ToDo: should check categories with modelDescription.xml
   */
  logging = 0;
  ret = IF_OK;
  #pragma omp master
  {
    ret = If_GetInt("mdl_logging", &logging);
  }
  if (ret != IF_OK) {
    ssSetErrorStatus(S, If_ResultString());
    return;
  }
  if (logging > ARRAY_SIZE(logOffsets) - 1)
    logging = ARRAY_SIZE(logOffsets) - 1;
  if (m->logging != logging) {
    if (logging > 0) {
      if ((*m->fmi2SetDebugLogging)(m->fmu, fmi2True,
           ARRAY_SIZE(logCategories) - logOffsets[logging],
           &logCategories[logOffsets[logging]]) > fmi2Warning) {
        ssSetErrorStatus(S, "can't set debug logging of FMU");
        return;
      }
    }
    else {
      if ((*m->fmi2SetDebugLogging)(m->fmu, fmi2False, 0, NULL) > fmi2Warning) {
        ssSetErrorStatus(S, "can't reset debug logging of FMU");
        return;
      }
    }
    m->logging = logging;
  }

  m->initPending = fmi2True; /* wait until inputs are available */
  m->nextEventTimeDefined = fmi2False;
}


/**
 *  Set parameters, states and inputs and obtain outputs y = f(t,p,x,u).
 */
static void mdlOutputs(SimStruct *S, int_T tid)
{
  Hxi_ModelData *m;
  fmi2EventInfo eventInfo;
  fmi2Boolean timeEvent;
  fmi2Boolean stateEvent;
  fmi2Boolean enterEventMode;
  fmi2Boolean terminateSimulation;
  int i, j;
  int idx[NUM_BASETYPES];
  InputRealPtrsType uPtrs;
  mxArray *param;
  double *vals, *xc;

  GET_MODELDATA(S, m);

  /* Set time after initialization */
  if (!m->initPending &&
      (*m->fmi2SetTime)(m->fmu, ssGetT(S)) != fmi2OK) {
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
  if (m->nxd > 0) {
    vals = ssGetRealDiscStates(S);
    for (j = 0; j < NUM_BASETYPES; j++)
      idx[j] = 0;
    for (i = 0; i < m->nxd; i++) {
      j = m->x->btv[i];
      switch (j) {
      case FMI_REAL:
        m->x->r[idx[j]++] = vals[i];
        break;
      case FMI_BOOLEAN:
        m->x->b[idx[j]++] = vals[i] == 0.0? fmi2False: fmi2True;
        break;
      case FMI_INTEGER:
        m->x->i[idx[j]++] = (fmi2Integer)(vals[i]);
        break;
      default:
        break;
      }
    }
    if (setValues(m, m->x) != fmi2OK) {
      ssSetErrorStatus(S, "can't set discrete states of FMU");
      return;
    }
  }
  if (m->nxc > 0) {
    vals = ssGetContStates(S);
    xc = m->dx->r;
    for (i = 0; i < m->nxc; i++)
      xc[m->dx_perm[i]] = *vals++;
    if ((*m->fmi2SetContinuousStates)(m->fmu, xc, m->nxc) != fmi2OK) {
      ssSetErrorStatus(S, "can't set continuous states of FMU");
      return;
    }
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
  timeEvent = fmi2False;
  stateEvent = fmi2False;
  if (m->initPending) {
    /* setup experiment and set start time */
    if ((*m->fmi2SetupExperiment)(m->fmu, fmi2False, 0.0,
                                  ssGetT(S), fmi2False, 0.0) != fmi2OK) {
      ssSetErrorStatus(S, "can't setup experiment for FMU");
      return;
    }
    /* initialize FMU */
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
      /* check for time and state event */
      if (m->nextEventTimeDefined && ssGetT(S) >= m->nextEventTime)
	timeEvent = fmi2True;
      if (m->nz > 0) {
	if ((*m->fmi2GetEventIndicators)(m->fmu, m->z, m->nz) != fmi2OK) {
	  ssSetErrorStatus(S, "can't get event indicators of FMU");
	  return;
	}
	for (i = 0; i < m->nz; i++) {
	  if ((m->z[i] >= 0.0 && m->pre_z[i] < 0.0) ||
	      (m->z[i] <= 0.0 && m->pre_z[i] > 0.0))
	    stateEvent = fmi2True;
	  m->pre_z[i] = m->z[i];
	}
      }
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
    if (m->nc > 0 && ssIsSampleHit(S, 1, tid) || m->nxc == 0)
      /* enter event mode in case of sample hit */
      /* always enter event mode if there are no continuous states */
      enterEventMode = fmi2True;
    if (enterEventMode || timeEvent || stateEvent)
      (*m->fmi2EnterEventMode)(m->fmu);
  }

  /* process events */
  if (enterEventMode || timeEvent || stateEvent) {
    if (m->nc > 0 && ssIsSampleHit(S, 1, tid) && m->fmi2SetClock != NULL) {
      /* set clock interval */
      if (m->fmi2SetInterval != NULL) {
        for (i = 0; i < m->nc; i++)
          m->civl[i] = ssGetSampleTime(S, 1 + i);
        if ((*m->fmi2SetInterval)(m->fmu, m->cidx, m->nc, m->civl) != fmi2OK) {
          ssSetErrorStatus(S, "can't set clock interval of FMU");
          return;
        }
      }
      /* don't update discrete states during continuous task */
      for (i = 0; i < m->nc; i++)
        m->csub[i] = ssIsContinuousTask(S, tid);
      if ((*m->fmi2SetClock)(m->fmu, m->cidx, m->nc, m->ctck, m->csub) != fmi2OK) {
        ssSetErrorStatus(S, "can't set clock of FMU");
        return;
      }
    }
    eventInfo.newDiscreteStatesNeeded = fmi2True;
    while (eventInfo.newDiscreteStatesNeeded) {
      (*m->fmi2NewDiscreteStates)(m->fmu, &eventInfo);
      if (eventInfo.terminateSimulation) {
        ssSetErrorStatus(S, "terminate in event mode of FMU");
        return;
      }
    }
    /* get next time event */
    m->nextEventTimeDefined = eventInfo.nextEventTimeDefined;
    if (m->nextEventTimeDefined)
      m->nextEventTime = eventInfo.nextEventTime;
    /* get initial event indicators */
    if (m->initPending) {
      if (m->nz > 0 &&
	  (*m->fmi2GetEventIndicators)(m->fmu, m->pre_z, m->nz) != fmi2OK) {
	ssSetErrorStatus(S, "can't get initial event indicators of FMU");
	return;
      }
      m->initPending = fmi2False;
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
 *  Update discrete states x(k) = f(x(k-1),u(k)).
 *  Note: mdlOutputs must have been called before to set states and inputs 
 */
static void mdlUpdate(SimStruct *S, int_T tid)
{
  Hxi_ModelData *m;
  double *vals, *xc;
  int idx[NUM_BASETYPES];
  int i, j;

  GET_MODELDATA(S, m);

  /* get states */
  if (m->nxd > 0) {
    vals = ssGetRealDiscStates(S);
    if (getValues(m, m->x) != fmi2OK) {
      ssSetErrorStatus(S, "can't get discrete states of FMU");
      return;
    }
    for (j = 0; j < NUM_BASETYPES; j++)
      idx[j] = 0;
    for (i = 0; i < m->nxd; i++) {
      j = m->x->btv[i];
      switch (j) {
      case FMI_REAL:
        vals[i] = m->x->r[idx[j]++];
        break;
      case FMI_BOOLEAN:
        vals[i] = m->x->b[idx[j]++]? 1.0: 0.0;
        break;
      case FMI_INTEGER:
        vals[i] = m->x->i[idx[j]++];
        break;
      default:
        break;
      }
    }
  }
  if (m->nxc > 0) {
    xc = m->dx->r;
    if ((*m->fmi2GetContinuousStates)(m->fmu, xc, m->nxc) != fmi2OK) {
      ssSetErrorStatus(S, "can't get continuous states of FMU");
      return;
    }
    vals = ssGetContStates(S);
    for (i = 0; i < m->nxc; i++)
      *vals++ = xc[m->dx_perm[i]];
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
  double *vals, *dx;
  int i;

  GET_MODELDATA(S, m);

  /* get derivatives */
  if (m->nxc > 0) {
    dx = m->dx->r;
    if ((*m->fmi2GetDerivatives)(m->fmu, dx, m->nxc) != fmi2OK) {
      ssSetErrorStatus(S, "can't get continuous states of FMU");
      return;
    }
    vals = ssGetdX(S);
    for (i = 0; i < m->nxc; i++)
      *vals++ = dx[m->dx_perm[i]];
  }
}

#define MDL_JACOBIAN
/**
 *  Obtain Jacobian J = d(dxc(t),xd(k),y)/d(xc(t),xd(k-1),u)
 */
static void mdlJacobian(SimStruct *S)
{
  Hxi_ModelData *m;
  double *pr; /* Jacobian elements */
  int *ir;    /* row indices */
  int *jc;    /* start index for each column */
  int i, j, jjac, jdx, offs; /* different indices */
  int idx[NUM_BASETYPES];
  InputRealPtrsType uPtrs;
  double *v, *w, dvj, vj_bak;
  fmi2ValueReference vrKnown[1];
  fmi2Real dvKnown[1];
  size_t nUnknown;

  GET_MODELDATA(S, m);

  /* get previous states */
  if (m->nxd > 0) {
    if (getValues(m, m->pre_x) != fmi2OK) {
      ssSetErrorStatus(S, "can't get previous states of FMU");
      return;
    }
    for (j = 0; j < NUM_BASETYPES; j++)
      idx[j] = 0;
    for (i = 0; i < m->nxd; i++) {
      j = m->x->btv[i];
      switch (j) {
      case FMI_REAL:
        m->pre_x_vals[i] = m->x->r[idx[j]++];
        break;
      case FMI_BOOLEAN:
        m->pre_x_vals[i] = m->x->b[idx[j]++]? 1.0: 0.0;
        break;
      case FMI_INTEGER:
        m->pre_x_vals[i] = m->x->i[idx[j]++];
        break;
      default:
        break;
      }
    }
  }

  /* get Jacobian variables */
  pr = ssGetJacobianPr(S);
  ir = ssGetJacobianIr(S);
  jc = ssGetJacobianJc(S);

  /* activate all clocks */
  for (i = 1; i <= m->nc; i++) {
    ssGetSampleHitPtr(S)[ssGetSampleTimeTaskID(S, i)] = 1;
  }

  /* obtain partial derivatives for state derivatives and outputs */
  if (m->nxc + m->ny > 0) {
    /* set clocks to subactive, i.e. no evaluation of states */
    ssGetSampleHitPtr(S)[ssGetSampleTimeTaskID(S, 0)] = 1;

    /* initialize current values */
    if (m->nxd > 0)
      memcpy(ssGetRealDiscStates(S), m->pre_x_vals, m->nxd*sizeof(fmi2Real));
    mdlOutputs(S, 0);
    if (m->nxc > 0)
      mdlDerivatives(S);
    if (m->nxd > 0)
      mdlUpdate(S, 0);

    if (m->providesDirectionalDerivative) {
      dvKnown[0] = 1.0;
      offs = m->nxc + m->nxd;
      for (j = 0; j < offs + m->nu; j++) {
        if (j < m->nxc)
          vrKnown[0] = m->x->vr[FMI_REAL][j];
        else if (j < offs)
          vrKnown[0] = m->pre_x->vr[FMI_REAL][j - m->nxc];
        else if (jc[j] == jc[j+1])
          continue; /* inactive input */
        else
          vrKnown[0] = m->u->vr[FMI_REAL][j - offs];
        nUnknown = 0;
        for (jdx = jc[j]; jdx < jc[j+1] && ir[jdx] < m->nxc; jdx++)
          m->vrUnknown[nUnknown++] = m->dx->vr[FMI_REAL][ir[jdx]];
        for (; jdx < jc[j+1] && ir[jdx] < offs; jdx++)
          /* list discrete states here because pr is filled directly */
          m->vrUnknown[nUnknown++] = m->x->vr[FMI_REAL][ir[jdx] - m->nxc];
        for (; jdx < jc[j+1] && ir[jdx] < m->nxc + m->ny; jdx++)
          m->vrUnknown[nUnknown++] = m->y->vr[FMI_REAL][ir[jdx] - offs];
        if (nUnknown > 0 && (*m->fmi2GetDirectionalDerivative)
            (m->fmu, m->vrUnknown, nUnknown, vrKnown, 1, dvKnown, &pr[jc[j]])
            != fmi2OK) {
          ssSetErrorStatus(S, "can't get directional derivative of FMU");
          return;
        }
      }
    }
    else {
      /* obtain finite differences */
      for (j = 0; j < m->nxc + m->nxd + m->nu; j++) {
        w = ssGetdX(S);
        for (jdx = jc[j]; jdx < jc[j+1] && ir[jdx] < m->nxc; jdx++)
          pr[jdx] = w[ir[jdx]];
        offs = m->nxc;
        w = ssGetRealDiscStates(S);
        for (; jdx < jc[j+1] && ir[jdx] < offs + m->nxd; jdx++)
          pr[jdx] = w[ir[jdx] - offs];
        offs += m->nxd;
        w = ssGetOutputPortRealSignal(S, 0);
        for (; jdx < jc[j+1]; jdx++)
          pr[jdx] = w[ir[jdx] - offs];
      }

      v = ssGetContStates(S);
      offs = m->nxc + m->nxd;
      for (j = 0; j < m->nxc; j++) {
        vj_bak = v[j];
        dvj = 1e-4 * fabs(vj_bak) + 1e-6;
        v[j] += dvj;

        if (m->nxd > 0)
          memcpy(ssGetRealDiscStates(S), m->pre_x_vals, m->nxd*sizeof(fmi2Real));
        mdlOutputs(S, 0);
        if (m->nxc > 0)
          mdlDerivatives(S);

        w = ssGetdX(S);
        for (jdx = jc[j]; jdx < jc[j+1] && ir[jdx] < m->nxc; jdx++)
          pr[jdx] = (w[ir[jdx]] - pr[jdx]) / dvj;
        for (; jdx < jc[j+1] && ir[jdx] < offs; jdx++)
          ;
        w = ssGetOutputPortRealSignal(S, 0);
        for (; jdx < jc[j+1] && ir[jdx] < offs + m->ny; jdx++)
          pr[jdx] = (w[ir[jdx] - offs] - pr[jdx]) / dvj;

        v[j] = vj_bak;
      }

      v = ssGetRealDiscStates(S);
      for (jjac = m->nxc, j = 0; j < m->nxd; jjac++, j++) {
        memcpy(v, m->pre_x_vals, m->nxd*sizeof(fmi2Real));
        vj_bak = v[j];
        dvj = 1e-4 * fabs(vj_bak) + 1e-6;
        v[j] += dvj;

        mdlOutputs(S, 0);
        if (m->nxc > 0)
          mdlDerivatives(S);

        w = ssGetdX(S);
        for (jdx = jc[jjac]; jdx < jc[jjac+1] && ir[jdx] < m->nxc; jdx++)
          pr[jdx] = (w[ir[jdx]] - pr[jdx]) / dvj;
        offs = m->nxc + m->nxd;
        for (; jdx < jc[jjac+1] && ir[jdx] < offs; jdx++)
          ;
        w = ssGetOutputPortRealSignal(S, 0);
        for (; jdx < jc[jjac+1] && ir[jdx] < offs + m->ny; jdx++)
          pr[jdx] = (w[ir[jdx] - offs] - pr[jdx]) / dvj;

        v[j] = vj_bak;
      }

      uPtrs = ssGetInputPortRealSignalPtrs(S, 0);
      for (jjac = m->nxc + m->nxd, j = 0; j < m->nu; jjac++, j++) {
        if (jc[jjac] == jc[jjac+1])
          continue; /* inactive input */
        vj_bak = *uPtrs[j];
        dvj = 1e-4 * fabs(vj_bak) + 1e-6;
        *uPtrs[j] += dvj;

        if (m->nxd > 0)
          memcpy(ssGetRealDiscStates(S), m->pre_x_vals, m->nxd*sizeof(fmi2Real));
        mdlOutputs(S, 0);
        if (m->nxc > 0)
          mdlDerivatives(S);

        w = ssGetdX(S);
        for (jdx = jc[jjac]; jdx < jc[jjac+1] && ir[jdx] < m->nxc; jdx++)
          pr[jdx] = (w[ir[jdx]] - pr[jdx]) / dvj;
        offs = m->nxc + m->nxd;
        for (; jdx < jc[jjac+1] && ir[jdx] < offs; jdx++)
          ;
        w = ssGetOutputPortRealSignal(S, 0);
        for (; jdx < jc[jjac+1] && ir[jdx] < offs + m->ny; jdx++)
          pr[jdx] = (w[ir[jdx] - offs] - pr[jdx]) / dvj;

        *uPtrs[j] = vj_bak;
      }
    }
  }

  /* obtain partial derivatives for discrete-time states */
  if (m->nxd > 0) {
    /* set clocks to active, i.e. evaluation of states */
    ssGetSampleHitPtr(S)[ssGetSampleTimeTaskID(S, 0)] = 0;

    /* initialize current values */
    memcpy(ssGetRealDiscStates(S), m->pre_x_vals, m->nxd*sizeof(fmi2Real));
    mdlOutputs(S, 0);
    mdlUpdate(S, 0);

    if (m->providesDirectionalDerivative) {
      dvKnown[0] = 1.0;
      offs = m->nxc + m->nxd;
      for (j = 0; j < m->nxc + m->nxd + m->nu; j++) {
        if (j < m->nxc)
          vrKnown[0] = m->x->vr[FMI_REAL][j];
        if (j < offs)
          vrKnown[0] = m->pre_x->vr[FMI_REAL][j - m->nxc];
        else if (jc[j] == jc[j+1])
          continue; /* inactive input */
        else
          vrKnown[0] = m->u->vr[FMI_REAL][j - offs];
        nUnknown = 0;
        for (jdx = jc[j]; jdx < jc[j+1] && ir[jdx] < m->nxc; jdx++)
          ;
        jjac = jdx;
        for (; jdx < jc[j+1] && ir[jdx] < offs; jdx++)
          m->vrUnknown[nUnknown++] = m->x->vr[FMI_REAL][ir[jdx] - m->nxc];
        if (nUnknown > 0 && (*m->fmi2GetDirectionalDerivative)
            (m->fmu, m->vrUnknown, nUnknown, vrKnown, 1, dvKnown, &pr[jjac])
            != fmi2OK) {
          ssSetErrorStatus(S, "can't get directional derivative of FMU states");
          return;
        }
      }
    }
    else {
      /* obtain finite differences */
      offs = m->nxc;
      for (j = 0; j < m->nxc + m->nxd + m->nu; j++) {
        for (jdx = jc[j]; jdx < jc[j+1] && ir[jdx] < offs; jdx++)
          ;
        w = ssGetRealDiscStates(S);
        for (; jdx < jc[j+1] && ir[jdx] < offs + m->nxd; jdx++)
          pr[jdx] = w[ir[jdx] - offs];
      }

      v = ssGetContStates(S);
      for (j = 0; j < m->nxc; j++) {
        vj_bak = v[j];
        dvj = 1e-4 * fabs(vj_bak) + 1e-6;
        v[j] += dvj;

        memcpy(ssGetRealDiscStates(S), m->pre_x_vals, m->nxd*sizeof(fmi2Real));
        mdlOutputs(S, 0);
        mdlUpdate(S, 0);

        for (jdx = jc[j]; jdx < jc[j+1] && ir[jdx] < offs; jdx++)
          ;
        w = ssGetRealDiscStates(S);
        for (; jdx < jc[j+1] && ir[jdx] < offs + m->nxd; jdx++)
          pr[jdx] = (w[ir[jdx] - offs] - pr[jdx]) / dvj;

        v[j] = vj_bak;
      }

      v = ssGetRealDiscStates(S);
      for (jjac = m->nxc, j = 0; j < m->nxd; jjac++, j++) {
        memcpy(v, m->pre_x_vals, m->nxd*sizeof(fmi2Real));
        vj_bak = v[j];
        dvj = 1e-4 * fabs(vj_bak) + 1e-6;
        v[j] += dvj;

        mdlOutputs(S, 0);
        mdlUpdate(S, 0);

        for (jdx = jc[jjac]; jdx < jc[jjac+1] && ir[jdx] < offs; jdx++)
          ;
        w = ssGetRealDiscStates(S);
        for (; jdx < jc[jjac+1] && ir[jdx] < offs + m->nxd; jdx++)
          pr[jdx] = (w[ir[jdx] - offs] - pr[jdx]) / dvj;

        v[j] = vj_bak;
      }

      uPtrs = ssGetInputPortRealSignalPtrs(S, 0);
      for (jjac = m->nxc + m->nxd, j = 0; j < m->nu; jjac++, j++) {
        if (jc[jjac] == jc[jjac+1])
          continue; /* inactive input */
        vj_bak = *uPtrs[j];
        dvj = 1e-4 * fabs(vj_bak) + 1e-6;
        *uPtrs[j] += dvj;

        memcpy(ssGetRealDiscStates(S), m->pre_x_vals, m->nxd*sizeof(fmi2Real));
        mdlOutputs(S, 0);
        mdlUpdate(S, 0);

        for (jdx = jc[jjac]; jdx < jc[jjac+1] && ir[jdx] < offs; jdx++)
          ;
        w = ssGetRealDiscStates(S);
        for (; jdx < jc[jjac+1] && ir[jdx] < offs + m->nxd; jdx++)
          pr[jdx] = (w[ir[jdx] - offs] - pr[jdx]) / dvj;

        *uPtrs[j] = vj_bak;
      }
    }
  }

  /* reset current values */
  /* set clocks to subactive, i.e. get outputs for previous states */
  ssGetSampleHitPtr(S)[ssGetSampleTimeTaskID(S, 0)] = 1;
  if (m->nxd > 0)
    memcpy(ssGetRealDiscStates(S), m->pre_x_vals, m->nxd*sizeof(fmi2Real));
  mdlOutputs(S, 0);
  if (m->nxc > 0)
    mdlDerivatives(S);
  if (m->nxd > 0)
    mdlUpdate(S, 0);

  /* deactivate all clocks */
  for (i = 1; i <= m->nc; i++) {
    ssGetSampleHitPtr(S)[ssGetSampleTimeTaskID(S, i)] = 0;
  }
}

/**
 *  Free model instance. Don't call fmi2Terminate because:
 *   - it is not allowed after a possible evaluation error
 *   - reuse the instance for subsequent evaluations
 */
static void mdlTerminate(SimStruct *S)
{
  /*
  Hxi_ModelData *m;

  GET_MODELDATA(S, m);

  if (!m->initPending && m->fmi2Terminate != NULL &&
      (*m->fmi2Terminate)(m->fmu) != fmi2OK) {
    ssSetErrorStatus(S, "can't terminate FMU");
    return;
  }
  */
}

/**
 * @}
 */

#ifdef  MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "simulink.c"      /* MEX-file interface mechanism */
#else
#include "cg_sfun.h"       /* Code generation registration function */
#endif

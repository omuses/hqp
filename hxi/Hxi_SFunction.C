/*
 * Hxi_SFunction.C:
 *   implementation of binary S-function interface for Hqp
 *
 * rf, 01/07/2005
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

#include "Hxi_SFunction.h"

#if defined(_MSC_VER) || defined(__MINGW32__)
#include <windows.h>
#else
#include <dlfcn.h>
#endif

// prototype for initialization function exported by binary object
typedef void (initFunction_t)(SimStruct *S);

//-------------------------------------------------------------------
SimStruct *Hxi_SimStruct_create()
{
  // allocate a SimStruct
  return new SimStruct();
}

//-------------------------------------------------------------------
void Hxi_SimStruct_destroy(SimStruct *S)
{
  if (S == NULL)
    return;

  // obtain handle and release S-function
#if defined(_MSC_VER) || defined(__MINGW32__)
  HMODULE handle = GetModuleHandle(ssGetPath(S));
  if (handle)
    FreeLibrary(handle);
#else
  void *handle = dlopen(ssGetPath(S), RTLD_LAZY);
  if (handle)
    dlclose(handle);
#endif

  // free memory
  delete S;
}

//-------------------------------------------------------------------
void mdlInitializeSizes(SimStruct *S)
{
  // get handle to binary S-function
#if defined(_MSC_VER) || defined(__MINGW32__)
  HMODULE handle = LoadLibrary(ssGetPath(S));
  if (!handle) {
    ssSetErrorStatus(S, "LoadLibrary failed");
    return;
  }
#else
  void *handle = dlopen(ssGetPath(S), RTLD_LAZY);
  if (!handle) {
    ssSetErrorStatus(S, dlerror());
    return;
  }
#endif

  // get pointer to entry point Hxi_SimStruct_init
  initFunction_t *initFunction_p;
#if defined(_MSC_VER) || defined(__MINGW32__)
  initFunction_p = (initFunction_t *)GetProcAddress(handle,
                                                    "Hxi_SimStruct_init");
  if (!initFunction_p) {
    ssSetErrorStatus(S, "GetProcAddress failed for Hxi_SimStruct_init");
    return;
  }
#else
  initFunction_p = (initFunction_t *)dlsym(handle, "Hxi_SimStruct_init");
  if (!initFunction_p) {
    ssSetErrorStatus(S, dlerror());
    return;
  }
#endif

  // initialize SimStruct
  (*initFunction_p)(S);

  // call S-function method mdlInitializeSizes
  (*(S->getmdlInitializeSizes()))(S);
}


//===================================================================

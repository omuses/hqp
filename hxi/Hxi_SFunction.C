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

#  if defined(HXI_WITH_MEX)
// prototype for MEX function
typedef void
(mexFunction_t)(int nlhs, mxArray *plhs[], int nrhs, mxArray *prhs[]);
#endif

//-------------------------------------------------------------------
extern "C" SimStruct *Hxi_SFunction_open(SimStruct *S)
{
  // get handle to binary S-function
#if defined(_MSC_VER) || defined(__MINGW32__)
  HMODULE handle = LoadLibrary(ssGetPath(S));
  if (!handle) {
    ssSetErrorStatus(S, "LoadLibrary failed");
    return S;
  }
#else
  void *handle = dlopen(ssGetPath(S), RTLD_LAZY);
  if (!handle) {
    ssSetErrorStatus(S, dlerror());
    return S;
  }
#endif

  // get pointer to entry point Hxi_SimStruct_init
  initFunction_t *initFunction_p;
#if defined(_MSC_VER) || defined(__MINGW32__)
  // try HXI S-function
  initFunction_p = (initFunction_t*)GetProcAddress(handle,
                                                   "Hxi_SimStruct_init");
  if (!initFunction_p) {
#  if !defined(HXI_WITH_MEX)
    ssSetErrorStatus(S, "GetProcAddress failed for Hxi_SimStruct_init");
    return S;
#  else
    // alternatively try MEX S-function
    mexFunction_t *mexFunction_p;
    mexFunction_p = (mexFunction_t*)GetProcAddress(handle, "mexFunction");
    if (!mexFunction_p) {
      ssSetErrorStatus(S, "GetProcAddress failed for Hxi_SimStruct_init"
                       " and mexFunction");
      return S;
    }
#  endif
  }
#else
  // try HXI S-function
  initFunction_p = (initFunction_t *)dlsym(handle, "Hxi_SimStruct_init");
  if (!initFunction_p) {
#  if !defined(HXI_WITH_MEX)
    ssSetErrorStatus(S, dlerror());
    return S;
#  else
    // alternatively try MEX S-function
    void *mexFunction_p = dlsym(handle, "mexFunction");
    if (!mexFunction_p) {
      ssSetErrorStatus(S, dlerror());
      return S;
    }
#  endif
  }
#endif

  if (initFunction_p)
    // initialize Hxi::SimStruct
    (*initFunction_p)(S);
#if defined(HXI_WITH_MEX)
  else
    S = NULL; // indicate that a MEX S-function is present
#endif

  return S;
}

//-------------------------------------------------------------------
extern "C" void Hxi_SFunction_close(SimStruct *S)
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
}


//===================================================================

/*
 * Hqp_DynLoad.C
 *   -- implementation
 *
 * rf, 3/13/97
 *  
 */

/*
    Copyright (C) 1994--2010  Ruediger Franke

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

#include "Hqp_DynLoad.h"

#if defined(_MSC_VER) || defined(__MINGW32__)
#include <windows.h>
#else
#include <dlfcn.h>
#endif

#include "Hqp.h";

//------------------------------------------------------------------------
Hqp_DynLoad::Hqp_DynLoad()
{
  _handle = NULL;
}

//------------------------------------------------------------------------
Hqp_DynLoad::~Hqp_DynLoad()
{
  if (_handle)
    close();
}

//------------------------------------------------------------------------
void Hqp_DynLoad::open(const char *pathname)
{
#if defined(_MSC_VER) || defined(__MINGW32__)
  _handle = (void *)LoadLibrary(pathname);
  if (!_handle) {
    m_error(E_INPUT, pathname);
  }
#else
  _handle = dlopen(pathname, RTLD_LAZY);
  if (!_handle) {
    m_error(E_INPUT, dlerror());
  }
#endif
}

//------------------------------------------------------------------------
void *Hqp_DynLoad::symbol(const char *name)
{
  void *sym = NULL;
#if defined(_MSC_VER) || defined(__MINGW32__)
  sym = (void *)GetProcAddress((HMODULE)_handle, name);
  if (!sym)
    m_error(E_NULL, name);
#else
  sym = dlsym(_handle, name);
  if (!sym)
    m_error(E_NULL, dlerror());
#endif
  return sym;
}

//------------------------------------------------------------------------
void Hqp_DynLoad::close()
{
#if defined(_MSC_VER) || defined(__MINGW32__)
  if (_handle)
    FreeLibrary((HMODULE)_handle);
#else
  if (_handle)
    dlclose(_handle);
#endif
}


//========================================================================

/*
 * Hqp_DynLoad.C
 *   -- implementation
 *
 * rf, 3/13/97
 *  
 */

/*
    Copyright (C) 1994--1998  Ruediger Franke

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
#include <dlfcn.h>

//------------------------------------------------------------------------
Hqp_DynLoad::Hqp_DynLoad()
{
  _handle = 0;
}

//------------------------------------------------------------------------
Hqp_DynLoad::~Hqp_DynLoad()
{
  if (_handle)
    close();
}

//------------------------------------------------------------------------
bool Hqp_DynLoad::open(const char *pathname)
{
  _handle = dlopen(pathname, RTLD_LAZY);
  return _handle? true: false;
}

//------------------------------------------------------------------------
void *Hqp_DynLoad::symbol(const char *name)
{
  return dlsym(_handle, name);
}

//------------------------------------------------------------------------
const char *Hqp_DynLoad::errmsg()
{
  return dlerror();
}

//------------------------------------------------------------------------
void Hqp_DynLoad::close()
{
  dlclose(_handle);
}


//========================================================================

/*
 *  If_Module.h
 *   - binds a pointer to a method to an interface command
 *   - class template (inline code to be compatible with most compilers)
 *
 *  rf, 7/19/94
 *
 */

/*
    Copyright (C) 1994--2000  Ruediger Franke

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

#ifndef If_Module_H
#define If_Module_H

#include "If_Command.h"
#include "If_Class.h"

// use this macro to create a new If_Module
//-----------------------------------------
#define IF_MODULE(id, ptr, base) \
  If_Module<base>(id, ptr, theIf_ClassList_##base)


template<class X>
class If_Module: public If_Command {

 protected:

  X	**_module;
  If_ClassList<X>*& _list;

  int 	invoke(int argc, char *argv[], char **result)
    {
      char *fallback;

      if (argc == 1) {
	if (*_module)
	  *result = (*_module)->name();
	else
	  *result = "No module chosen.";
	return IF_OK;
      }
      else if (argc == 2) {
	if (!_list) {
	  *result = "No modules available.";
	  return IF_ERROR;
	}

	// store the current module name
	if (*_module)
	  fallback = (*_module)->name();
	else
	  fallback = NULL;

	// delete the current module to free interface elements
	delete *_module;

	// create the requested module
	*_module = _list->createObject(argv[1]);
	if (!*_module) {
	  *result = (char *)_list->cand_names();
	  if (fallback)
	    *_module = _list->createObject(fallback);
	  return IF_ERROR;
	}
	return IF_OK;
      }

      *result = "If_Module: wrong # of args";
      return IF_ERROR;
    }	

 public:

  If_Module(const char *ifName, X **module, If_ClassList<X>*& list)
    :If_Command(ifName),
      _list(list)
      {
	_module = module;
      }
};


#endif

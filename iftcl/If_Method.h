/*
 *  If_Method.h
 *   - binds a pointer to a method to an interface command
 *   - class template (inline code to be compatible with most compilers)
 *
 *  rf, 7/19/94
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

#ifndef If_Method_H
#define If_Method_H

#include "If_Command.h"


template<class X>
class If_Method: public If_Command {

 protected:

  X	*_object;
  int	(X::*_method)(int, char *[], char **);

  int 	invoke(int argc, char *argv[], char **result)
    {
      return (_object->*_method)(argc, argv, result);
    }	

 public:

  If_Method(char *ifName,
	    int (X::*method)(int, char *[], char **), X *object)
    :If_Command(ifName)
      {
	_object = object;
	_method = method;
      }
};


#endif

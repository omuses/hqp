/*
 *  If_Proc.h
 *   - binds a function pointer to an interface command
 *
 *  rf, 7/20/94
 *
 */

/*
    Copyright (C) 1994--2002  Ruediger Franke

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

#ifndef If_Proc_H
#define If_Proc_H

#ifdef __GNUC__
// (standard C++ does not define #warning)
#warning "If_Proc is deprecated and will be removed; use If_Procedure or If_String instead!"
#endif

#include "If_Command.h"


typedef int If_Proc_t(int argc, char *argv[], char **result);

//--------------------------------
class IF_API If_Proc: public If_Command {

 protected:

  If_Proc_t	*_proc;

  int 	invoke(int argc, char *argv[], char **);

 public:

  If_Proc(const char *ifName, If_Proc_t *proc);
};


#endif

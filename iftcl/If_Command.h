/*
 *  If_Command.h
 *   - abstract base class for interface commands
 *
 *  rf, 7/20/94
 *
 */

/*
    Copyright (C) 1994--2001  Ruediger Franke

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

#ifndef If_Command_H
#define If_Command_H

#include "If_Element.h"

#define IF_CMD_ARGS	int, char *[], char **
#define IF_DEF_ARGS	int argc=0, char *argv[]=NULL, char **result=NULL


//-----------------------------------
class If_Command: public If_Element {

 protected:

  If_Command(const char *ifName);

  // interface to Tcl
  //-----------------
  static int 	tclCmd(ClientData, Tcl_Interp *, int argc, char *argv[]);

  // interface to derived classes
  //-----------------------------
  virtual int	invoke(int argc, char *argv[], char **result)=0;

 public:

  ~If_Command();
};


#endif

/*
 *  If_Element.C -- class definition
 *
 *  rf, 6/22/94
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

#include <stdlib.h>
#include <string.h>

#include "If_ListElement.h"
#include "If_Element.h"

//--------------------------------------------------------------------------
Tcl_Interp *theInterp = NULL;

//--------------------------------------------------------------------------
If_Element::If_Element(const char *ifName)
{
  _ifName = strdup(ifName);
}

//--------------------------------------------------------------------------
If_Element::~If_Element()
{
  free(_ifName);
}


//==========================================================================

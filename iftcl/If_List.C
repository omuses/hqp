/*
 * If_List.C -- class definition
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

#include <assert.h>

#include "If_ListElement.h"
#include "If_List.h"


//-----------------------------------------------------------------------------
If_List::If_List()
{
  _top= NULL;
}

//-----------------------------------------------------------------------------
If_List::~If_List()
{
  while (_top)
    delete _top;
}

//-----------------------------------------------------------------------------
void If_List::append(If_ListElement *element)
{
  assert(element->_list == NULL);  // may be appended to one list at once

  element->_prev= _top;
  element->_next= NULL;

  if (_top)
    _top->_next= element;
  _top= element;

  element->_list= this;
}

//-----------------------------------------------------------------------------
void If_List::release(If_ListElement *element)
{
//printf("list: %d, this: %d, element: %d\n", element->_list, this, element);
  if (!element->_list)
    return;

  assert(element->_list == this);  // is appended to this list
  
  if(element == _top)
    _top= element->_prev;
  
  if(element->_prev)
    element->_prev->_next= element->_next;

  if(element->_next)
    element->_next->_prev= element->_prev;

  element->_prev= element->_next= NULL;
  element->_list= NULL;
}


//=============================================================================

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
  _anchor = new If_ListElement();
}

//-----------------------------------------------------------------------------
If_List::~If_List()
{
  while (_anchor->_prev)
    delete _anchor->_prev;
  delete _anchor;
}

//-----------------------------------------------------------------------------
void If_List::append(If_ListElement *element)
{
  // element must not be appended to a list (can at most be appended to one)
  assert(element && !element->_prev && !element->_next);

  // NULL <n--anchor--p> <n--element--p> <n--el_old--p> ... 
  element->_next = _anchor;
  element->_prev = _anchor->_prev;
  _anchor->_prev = element;
  if (element->_prev)
    element->_prev->_next = element;
}

//-----------------------------------------------------------------------------
void If_List::release(If_ListElement *element)
{
  // search element in this list
  If_ListElement *el = top();
  while (el && el != element)
    el = prev(el);
  element = el;
  assert(element);  // element must be appended to this list

  if(element->_prev)
    element->_prev->_next = element->_next;

  if(element->_next)
    element->_next->_prev = element->_prev;

  element->_prev = element->_next = NULL;
}


//=============================================================================

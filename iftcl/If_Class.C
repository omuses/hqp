/*
 * If_Class.C --
 *   basic definitions for managing interface modules
 *
 * rf, 10/5/95
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

#include <malloc.h>
#include <stdio.h>

#include "If_Class.h"

//
// If_idList
//

If_idList::If_idList()
{
  _cand_names_size = 100;
  _cand_names = (char *)malloc(_cand_names_size);
  sprintf(_cand_names, "Candidates: ");
}

If_idList::~If_idList()
{
  free(_cand_names);
}

void If_idList::append(If_idListElement *element)
{
  size_t len, add;

  // store the element's name
  len = strlen(_cand_names) + 1; 	// one more for terminating null
  add = strlen(element->if_id()) + 2; 	// two more for seperator ", "
  if (len + add > _cand_names_size) {
    _cand_names_size += (add > 100)? add: 100;
    _cand_names = (char *)realloc(_cand_names, _cand_names_size);
  }
  if (top() != NULL)
    sprintf(_cand_names + len - 1, ", %s", (const char *)element->if_id());
  else
    sprintf(_cand_names + len - 1, "%s", (const char *)element->if_id());

  // append the element
  If_List::append(element);
}


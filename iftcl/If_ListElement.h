/*
 *  If_ListElement.h --
 *     - stores a pointer of a If_List
 *     - stores two pointer for If_List for concatting elements
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

#ifndef If_ListElement_H
#define If_ListElement_H

#include <string.h>

// fast string comparison
//-----------------------
inline int if_strcmp(const char *str1, const char *str2, int pos=0)
{
  return !(str1[pos] == str2[pos] && strcmp( str1, str2) == 0);
}

//--------------------
class If_ListElement {
   friend class If_List;

 protected:
   If_List		*_list;		// optional list of that element
   If_ListElement	*_prev;		// previous element
   If_ListElement	*_next;		// next element

 public:    
   If_ListElement();
   virtual ~If_ListElement();
};


#endif

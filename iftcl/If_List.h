/**
 *  @file If_List.h
 *     List of If_ListElement's
 *
 *   rf, 22/6/94
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

#ifndef If_List_H
#define If_List_H

#include "If.h"
#include "If_ListElement.h"

/**
 *  List of If_ListElement's. The list destructor destroys all elements.
 */
class IF_API If_List {

 protected:
  If_ListElement	*_anchor;	///< point to top of the list

 public:
  If_List(); 				///< constructor
  ~If_List(); 				///< destructor

  void append(If_ListElement *); 	///< append an element
  void release(If_ListElement *); 	///< release an element

  /** get top element */
  If_ListElement *top() {return _anchor->_prev;}

  /** get previous element */
  If_ListElement *prev(If_ListElement *elem) {return elem->_prev;}
};

#endif

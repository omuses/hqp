/*
 *
 *  If_List.h -- 
 *     - stores a pointer to a list of If_ListElement's
 *     - many list elements may be concatted
 *     - destructor deletes all elements
 *
 *   rf, 22/6/94, rev. 1.0
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

#ifndef If_List_H
#define If_List_H

#include "If_ListElement.h"

#ifndef NULL
#define NULL 0
#endif

class If_List {

  protected:
    If_ListElement	*_top;		// top of that list

  public:
    If_List();
    ~If_List();

    void		append(If_ListElement *);
    void		release(If_ListElement *);

    If_ListElement	*top() {return _top;}
    If_ListElement	*prev(If_ListElement *elem) {return elem->_prev;}
};

#endif

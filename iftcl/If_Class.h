/*
 * If_Class.h --
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

#ifndef If_Class_H
#define If_Class_H

#include "If_id.h"
#include "If_List.h"
#include "If_ListElement.h"

//
// macros for defining interface modules
//
// IF_BASE_DECLARE: usage in the header file of a base class
// IF_BASE_DEFINE:  usage in the implementation file of a base class
// IF_CLASS_DEFINE: usage in the implementation file of a derived class
// IF_CLASS_ALLOC:  manually allocate on machines without static constructors
//

#define IF_BASE_DECLARE(base) \
  If_ClassList<base> *&If_ClassList_##base(); \
  extern If_ClassList<base> *theIf_ClassList_##base;

#define IF_BASE_DEFINE(base) \
  If_ClassList<base> *theIf_ClassList_##base = NULL; \
  If_ClassList<base> *&If_ClassList_##base() { \
    return theIf_ClassList_##base; \
  }

#ifndef IF_CLASS_STATIC
#define IF_CLASS_DEFINE(id, type, base) \
  static If_Class<type,base> if_class_##type(id, If_ClassList_##base());
#else
#define IF_CLASS_DEFINE(id, type, base)
#define IF_CLASS_ALLOC(id, type, base) \
  new If_Class<type,base>(id, If_ClassList_##base());
#endif

//
// If_List
//   (see If_List.h)
// |
// If_idList
//   manage a string with If_id's of all elements of the list
// |
// If_ClassList<BASE>
//   one list should be defined for each base class of similar modules
//

class If_idListElement;

class IF_API If_idList: public If_List {
 protected:
  char		*_cand_names;
  size_t	_cand_names_size;

 public:
  If_idList();
  ~If_idList();

  void append(If_idListElement *);
  const char *cand_names() {return _cand_names;}
};

template <class BASE>
class If_BaseClass;

template <class BASE>
class If_ClassList: public If_idList {
 public:
  BASE *createObject(const If_id id)
    {
      If_BaseClass<BASE> *el;
      
      el = (If_BaseClass<BASE> *)top();
      while(el && el->if_id() != id)
        el = (If_BaseClass<BASE> *)prev(el);

      if (!el)
        return NULL;

      return el->createObject();
    }
};

//
// If_ListElement
//   (see If_ListElement.h)
// |
// If_idListElement
//   - store an If_id for the element
//
// |
// If_BaseClass<BASE>
//   - constructor creates the list and appends elements
//   - declare abstract method for object creation of type BASE
// |
// If_Class<TYPE, BASE>
//   - define abstract method for object creation, using type TYPE
//

class If_idListElement: public If_ListElement {
 protected:
  If_id _id;
 public:
  If_idListElement(If_id id)
    : _id(id) {}
  If_id if_id()
    {
      return _id;
    }
};

template <class BASE>
class If_BaseClass: public If_idListElement {
 public:
  If_BaseClass(If_id id, If_ClassList<BASE>*& list)
    : If_idListElement(id)
    {
      if (!list)
	list = new If_ClassList<BASE>;
      list->append(this);
    }
  virtual BASE *createObject() = 0;
};

template <class TYPE, class BASE>
class If_Class: public If_BaseClass<BASE> {
 public:
  If_Class(If_id id, If_ClassList<BASE>*& list)
    : If_BaseClass<BASE>(id, list) {};
  BASE *createObject()
    {
      return new TYPE;
    }
};


#endif

/*
 * fundamental declarations for a framework interface
 *
 * rf, 10/4/95
 */

/*
    Copyright (C) 1994--2000  Ruediger Franke

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

#ifndef If_id_H
#define If_id_H

#include <stdlib.h>
#include <string.h>

/*
 * If_id to address a module trough the interface
 */

class If_id {
 protected:
  char *_id;
 public:
  If_id(const char *id) {_id = strdup(id);}
  If_id(const If_id &id) {_id = strdup(id._id);}
  ~If_id() {free(_id);}
  operator const char* () {return _id;}
};

inline int operator == (If_id id1, If_id id2)
{
  return (strcmp(id1, id2) == 0);
}

inline int operator != (If_id id1, If_id id2)
{
  return (strcmp(id1, id2) != 0);
}

#endif

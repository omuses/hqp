/*
 * Hqp_DynLoad.h
 *   -- dynamically load symbols
 *
 * rf, 3/13/97
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

#ifndef Hqp_DynLoad_H
#define Hqp_DynLoad_H

class Hqp_DynLoad {

protected:
  void *_handle;

public:
  Hqp_DynLoad();
  ~Hqp_DynLoad();

  bool open(const char *pathname);
  void *symbol(const char *name);
  const char *errmsg();
  void close();
};


#endif

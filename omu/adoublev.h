/**
 * @file adoublev.h
 *    Basic operations for vectors of adouble variables (adoublev replacement).
 *
 * e_arnold, 2010-06-01
 *           2010-06-16 max_size, alloc
 *           2011-09-12 operator =
 */

/*
    Copyright (C) 2010--2014  Eckhard Arnold

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

#ifndef ADOUBLEV_H
#define ADOUBLEV_H

#ifdef OMU_WITH_ADOLC
#include <adolc/adouble.h>
#else
typedef double adouble;
#endif

class adoublev
{
 protected:
  adouble *v;
  int size;                         /* size of the vector */
  int max_size;                     /* max size of the vector */

 public:
  adoublev(void );
  ~adoublev(void );
  void alloc(int );                 /* allocate memory */
  int sz() const {return size;}     /* Get the size of the vector */
  adoublev& operator >>= (double* );
  adoublev& operator <<= (double* );
  adouble& operator[](int ) const;  /* Can access component like an array */
  adoublev& operator = (const adoublev& );
  adoublev& operator -= (const adoublev& );
  adoublev& operator += (const adoublev& );
};

#endif

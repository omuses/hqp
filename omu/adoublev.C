//
// adoublev.C --
//   -- class definition
//
// e_arnold, 2010-06-01
//           2010-06-16 max_size, alloc
//           2011-09-12 operator =
//

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

#include "adoublev.h"
#include "Meschach.h"

//--------------------------------------------------------------------------
adoublev::adoublev(void)
{
  size = max_size = 0;
  v = (adouble *) NULL;
}

//--------------------------------------------------------------------------
adoublev::~adoublev(void)
{ 
  if ( ( max_size > 0 ) && (v != (adouble *) NULL ) )
    delete[] v;
  size = max_size = 0;
  v = (adouble *) NULL;
}

//--------------------------------------------------------------------------
void adoublev::alloc(int sz)
{
  if ( sz <= max_size )
    size = sz;
  else {
    if ( ( max_size > 0 ) && (v != (adouble *) NULL ) )
      delete[] v;
    size = max_size = sz;
    v = new adouble[size]; 
  }
}

//--------------------------------------------------------------------------
adouble& adoublev::operator[](int i) const
{
  if ( i < 0 || i >= size )
    m_error(E_SIZES, "adoublev::operator[]");
  return v[i];
}

//--------------------------------------------------------------------------
adoublev& adoublev::operator <<= (double* y)
{
  int i;
#ifdef OMU_WITH_ADOLC
  for ( i = 0; i < size; i++ )
    v[i] <<= y[i];
#else
  m_error(E_NULL, "adoublev::operator <<=: was compiled without ADOL-C");
#endif
  return *this;
}

//--------------------------------------------------------------------------
adoublev& adoublev::operator >>= (double* y)
{
  int i;
#ifdef OMU_WITH_ADOLC
  for ( i = 0; i < size; i++ )
    v[i] >>= y[i];
#else
  m_error(E_NULL, "adoublev::operator >>=: was compiled without ADOL-C");
#endif
  return *this;
}

//--------------------------------------------------------------------------
adoublev& adoublev::operator = (const adoublev& y)
{
  int i;
  if(y.size != size)
    m_error(E_SIZES, "adoublev::operator =");

  for ( i = 0; i < size; ++i )
    v[i] = y.v[i];
  return *this;
}

//--------------------------------------------------------------------------
adoublev& adoublev::operator -= (const adoublev& y)
{
  int i;
  if(y.size != size)
    m_error(E_SIZES, "adoublev::operator -=");

  for ( i = 0; i < size; ++i )
    v[i] -= y.v[i];
  return *this;
}

//--------------------------------------------------------------------------
adoublev& adoublev::operator += (const adoublev& y)
{
  int i;
  if(y.size != size)
    m_error(E_SIZES, "adoublev::operator +=");

  for ( i = 0; i < size; ++i )
    v[i] += y.v[i];
  return *this;
}

//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

/*
 * t_mesch.h --
 *   - 'time-dependent' Meschach data structures:
 *
 *        TVECP   - vector of VECP
 *        TIVECP  - vector of IVECP
 *        PERMP   - PERM *
 *        TPERMP  - vector of PERMP
 *        TMATP   - vector of MATP
 *
 * E. Arnold  07/03/97
 *            2002-04-17 free() replaced by tfree()
 *            
 */

/*
    Copyright (C) 1997--2002  Eckhard Arnold

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


#ifndef T_MESCH_H
#define T_MESCH_H

#include "Meschach.h"

//----------   TVECP   ----------

class TVECP {
  
 private:
  int _kmax;
  VECP *_x;

 public:
  TVECP();
  ~TVECP();
  void tfree();
  void resize(int, int = 0);
  void resize(const IVECP);
  VECP & operator [] (int k)
    {
#ifdef DEBUG
      if ( ( k < 0 ) || ( k > _kmax ) )
	error(E_BOUNDS,"TVECP::operator[]");
#endif
      return _x[k];
    }
  void Print();
  void Print(char *s) { printf("TVECP %s\n", s); Print(); }
};

//----------   TIVECP   ----------

class TIVECP {
  
 private:
  int _kmax;
  IVECP *_x;

 public:
  TIVECP();
  ~TIVECP();
  void tfree();
  void resize(int, int = 0);
  IVECP & operator [] (int k)
    {
#ifdef DEBUG
      if ( ( k < 0 ) || ( k > _kmax ) )
	error(E_BOUNDS,"TIVECP::operator[]");
#endif
      return _x[k];
    }
  void Print();
  void Print(char *s) { printf("TIVECP %s\n", s); Print(); }
};

//----------   TPERMP   ----------

class TPERMP {
  
 private:
  int _kmax;
  PERMP *_x;

 public:
  TPERMP();
  ~TPERMP();
  void tfree();
  void resize(int, int = 0);
  PERMP & operator [] (int k)
    {
#ifdef DEBUG
      if ( ( k < 0 ) || ( k > _kmax ) )
	error(E_BOUNDS,"TPERMP::operator[]");
#endif
      return _x[k];
    }
  void Print();
  void Print(char *s) { printf("TPERMP %s\n", s); Print(); }
};

//----------   TMATP   ----------

class TMATP {
  
 private:
  int _kmax;
  MATP *_x;

 public:
  TMATP();
  ~TMATP();
  void tfree();
  void resize(int, int = 0, int = 0);
  MATP & operator [] (int k)
    {
#ifdef DEBUG
      if ( ( k < 0 ) || ( k > _kmax ) )
	error(E_BOUNDS,"TMATP::operator[]");
#endif
      return _x[k];
    }
  void Print();
  void Print(char *s) { printf("TMATP %s\n", s); Print(); }
};

#endif

/*
 * Hqp_HL_Gangster.h -- 
 *   - Hessian approximation according Broyden, Fletsher, Goldfarb, Shanno
 *   - modification for positive definiteness according Powell
 *   - uses Gangster operator to avoid fill-in
 *
 * rf, 7/19/94
 */

/*
    Copyright (C) 1994--2014  Ruediger Franke

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

#ifndef Hqp_HL_Gangster_H
#define Hqp_HL_Gangster_H

#include "Hqp_HL.h"


class Hqp_HL_Gangster: public Hqp_HL {

 protected:
  Real  _gamma;
 
  VEC	*_v;
  VEC	*_sQ;
  VEC	*_Qs;

 public:
  Hqp_HL_Gangster();
  ~Hqp_HL_Gangster();

  void update_Q(const VEC *s, const VEC *u, Hqp_SqpProgram *);

  const char *name() {return "Gangster";}
};  


#endif

/**
 * @file Hqp_impl.h
 *   internal declarations for HQP (Huge Quadratic Programming)
 *
 * rf, 5/28/94
 *
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

#ifndef Hqp_impl_H
#define Hqp_impl_H

#include "Hqp.h"

/** return HQP version string */
extern const char *hqp_version();

/** possible results of QP solver */
enum Hqp_Result {
  Hqp_Optimal = 0,
  Hqp_Feasible,
  Hqp_Infeasible,
  Hqp_Suboptimal,
  Hqp_Degenerate
};

/** array of result strings for Hqp_Result */
extern const char *hqp_result_strings[];


#endif

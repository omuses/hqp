/**
 * @file Hqp_omp.h
 *   interface to OpenMP
 *
 * rf, 6/29/17
 *
 */

/*
    Copyright (C) 1994--2017  Ruediger Franke

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

#ifndef Hqp_omp_H
#define Hqp_omp_H

#ifdef HQP_WITH_OMP
#include <omp.h>
// Note: may reduce CPUs with environment variable OMP_NUM_THREADS
#else
// Default OMP functions for single threading
inline int omp_get_max_threads()
{
  return 1;
}

inline int omp_get_thread_num()
{
  return 0;
}
#endif

#endif

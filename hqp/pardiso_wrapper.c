/**
 * @file pardiso_wrapper.C
 *   implement interface to MKL PARDISO
 *
 * rf, 2017/10/07
 *
 */

/*
    Copyright (C) 2006--2017   Hartmut Linke and Ruediger Franke

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

#include "pardiso_wrapper.h"

/** MKL PARDISO function */
pardiso_ft mkl_pds_pardiso;

void pardiso_wrapper(void *pt, long *maxfct, long *mnum, long *mtype,
                     long *phase, long *n, double *a, long *ia, long *ja,
                     long *perm, long *nrhs, long *iparm, long *msglvl,
                     double *b, double *x, long *error)
{
  mkl_pds_pardiso(pt, maxfct, mnum, mtype, phase, n, a, ia, ja, perm, nrhs,
                  iparm, msglvl, b, x, error);
}

/**
 * @file pardiso_wrapper.h
 *   wrap call to PARallel DIrect SOlver PARDISO to enable
 *   pre-compilation from mkl libs and dynamic load at runtime
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

#ifndef PARDISO_WRAPPER_H
#define PARDISO_WRAPPER_H

/** Signature of used PARDISO function (version 3) */
typedef void (pardiso_ft)
	(void *, long *, long *, long *, long *, long *,
	double *, long *, long *, long *, long *, long *,
	long *, double *, double *, long *);

/** define PARDISO_API when compiling a Dynamic Link Library (DLL) */
#ifndef PARDISO_API
#define PARDISO_API 
#endif

#ifdef __cplusplus
extern "C" {
#endif
  
/** Wrapper for PARDISO function */
PARDISO_API pardiso_ft pardiso_wrapper;

#ifdef __cplusplus
};
#endif

#endif

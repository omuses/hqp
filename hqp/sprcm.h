/*
 * Copyright (C) 1995 R"udiger Franke, all rights reserved.
 *
 * This code is based on the Meschach library (Copyright (C) 1993 
 * David E. Steward & Zbigniew Leyk, all rights reserved.).
 * The Meschach copyright conditions, as included bolow, extend to it.
 *
 **************************************************************************
 *
 *			     Meschach Library
 * 
 * This Meschach Library is provided "as is" without any express 
 * or implied warranty of any kind with respect to this software. 
 * In particular the authors shall not be liable for any direct, 
 * indirect, special, incidental or consequential damages arising 
 * in any way from use of the software.
 * 
 * Everyone is granted permission to copy, modify and redistribute this
 * Meschach Library, provided:
 *  1.  All copies contain this copyright notice.
 *  2.  All modified copies shall carry a notice stating who
 *      made the last modification and the date of such modification.
 *  3.  No charge is made for this software or works derived from it.  
 *      This clause shall not be construed as constraining other software
 *      distributed on the same medium as this software, nor is a
 *      distribution fee considered a charge.
 *
 **************************************************************************
 *
 * function prototypes for Reverse Cuthill McKee ordering of sparse matrices
 *
 * rf, 11/15/95
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

PERM *sp_symrcm(const SPMAT *Q, PERM *order);
PERM *sp_kktrcm(const SPMAT *Q, const SPMAT *A, const SPMAT *C,
		PERM *order);

SPMAT *pxinv_sprows(const PERM *px, const SPMAT *src, SPMAT *dst);
SPMAT *pxinv_spcols(const PERM *px, const SPMAT *src, SPMAT *dst);

VEC *pxinv_vec(const PERM *px, const VEC *src, VEC *dst);

IVEC *sp_rcm_scan(const SPMAT *Q, const SPMAT *A, const SPMAT *C,
		  IVEC *degree, IVEC *neigh_start, IVEC *neighs);
PERM *sp_rcm_order(const IVEC *degree, const IVEC *neigh_start,
		   const IVEC *neighs, PERM *order);
int sp_rcm_sbw(const IVEC *neigh_start, const IVEC *neighs, const PERM *order);

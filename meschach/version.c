
/**************************************************************************
**
** Copyright (C) 1993 David E. Steward & Zbigniew Leyk, all rights reserved.
**
**			     Meschach Library
** 
** This Meschach Library is provided "as is" without any express 
** or implied warranty of any kind with respect to this software. 
** In particular the authors shall not be liable for any direct, 
** indirect, special, incidental or consequential damages arising 
** in any way from use of the software.
** 
** Everyone is granted permission to copy, modify and redistribute this
** Meschach Library, provided:
**  1.  All copies contain this copyright notice.
**  2.  All modified copies shall carry a notice stating who
**      made the last modification and the date of such modification.
**  3.  No charge is made for this software or works derived from it.  
**      This clause shall not be construed as constraining other software
**      distributed on the same medium as this software, nor is a
**      distribution fee considered a charge.
**
***************************************************************************/


/*			Version routine			*/
/*	This routine must be modified whenever modifications are made to
	Meschach by persons other than the original authors
	(David E. Stewart & Zbigniew Leyk); 
	when new releases of Meschach are made the
	version number will also be updated
*/

#include	<stdio.h>

void	m_version()
{
	static char rcsid[] = "$Id: version.c,v 1.3 2002/04/23 14:17:35 rfranke Exp $";

	printf("Meshach matrix library version 1.2b\n");
	printf("RCS id: %s\n",rcsid);
	printf("Changes since 1.2a:\n");
	printf("\t Fixed bug in schur() for 2x2 blocks with real e-vals\n");
	printf("\t Fixed bug in schur() reading beyond end of array\n");
	printf("\t Fixed some installation bugs\n");
	printf("\t Fixed bugs & improved efficiency in spILUfactor()\n");
	printf("\t px_inv() doesn't crash inverting non-permutations\n");
	/**** List of modifications ****/
	/* Example below is for illustration only */
	/* printf("Modified by %s, routine(s) %s, file %s on date %s\n",
			"Joe Bloggs",
			"m_version",
			"version.c",
			"Fri Apr  5 16:00:38 EST 1994"); */
	/* printf("Purpose: %s\n",
			"To update the version number"); */
	printf("Changes made by Ruediger Franke for HQP\n");
	printf("\tmatrix.h, sparse.h: brought function prototypes to ANSI C\n");
	printf("\tmatrix.h: increased MAXDIM to M_MAX_INT to work around a bug in iv_free()\n");
	printf("\tvecop.c, v_resize(): fixed memory allocation bug for dim == 0\n");
	printf("\tpxop.c, pxinv_vec(): replaced with a new version\n");
	printf("\tmachine.h: introduced defining SPARSE_COL_ACCESS for conditional compilation of column accesses\n");
	printf("\tsparse.h, sparse.c, sparseio.c, sprow.c: excluded code for column accesses with conditional compilation\n");
	printf("\tmakefile: excluded spchfctr.o, splufctr.o, spbkp.o, and spswap.o from the library\n");
	/* rf, 2/17/98 */
	printf("\tcopy.c: _m_copy(), _v_copy() call resize() for any changing dimension, not only for increase\n");
	/* rf, 6/18/98 */
	printf("\tmatrixio.c: bm_finput(), bfin_vec(): call resize() for changing dimension, not only for NULL pointer\n");
	/* rf, 9/11/98 */
	printf("\tsparse.h, sparsdef.h: replaced pair with SPPAIR to avoid name conflict with STL (Par Winzell)\n");
	/* rf, 10/27/00 */
	printf("\thessen.c, Hfactor: check for numerical overflow after call to hvec() in hsehldr.c to avoid infinite loop in symmeig()\n");
	/* rf, 12/05/00 */
	printf("\tmove err.h to m_err.h (name conflict with /usr/include/err.h)\n");
	printf("\tm_err.h: rename macro \"catch\" to \"m_catch\" (conflict with C++ keyword)\n");
	/* rf, 04/02/01 */
	printf("\tmatrix.h, matdef.h: allocate at least 1 element in NEW_A (as some c/malloc's return NULL for 0 elements)\n");
	/* rf, 04/23/02 */
	printf("\tmakefile: add add splufctr to list of compiled objects\n");
	printf("\tsplufctr.c: disable spLUTsolve and spILUfactor via SPARSE_COL_ACCESS\n");
	printf("\n");
}

/* $Log: version.c,v $
/* Revision 1.3  2002/04/23 14:17:35  rfranke
/* include additional factorization routines
/*
/* Revision 1.2  2001/04/02 08:18:27  rfranke
/* 	- meschach/{matrix, matdef}.h: allocate at least 1 element in NEW_A,
/* 	  as some malloc implementations return NULL pointer for 0 elements
/* 	  (braught up again by Christan Schulz)
/*
/* Revision 1.1.1.1  2001/03/01 17:19:15  rfranke
/* Import of version 1.7b1 taken from odoaker.
/*
/* Revision 1.4  2000/12/05 09:42:57  rf
/* 12/05/00: rf: resolve meschach/err.h, catch conflict
/* 	- configure.in: change version to 1.7a5
/* 	- rename meschach/err.h to meschach/m_err.h
/* 	- m_err.h: rename macro catch to m_catch
/* 	- adapt meschach/{matrix.h, err.c, torture.c, ztorture.c, makefile,
/* 	  version.c}
/* 	- hqp/Meschach.h: don't need to redefine catch anymore
/* 	- hqp/Hqp_IpsMehrotra.C, hqp/HqpIpsFranke.C: use m_catch
/*
/* Revision 1.3  2000/11/01 21:39:02  rf
/* 11/01/00: rf:
/* 	- bug fix in meschach/hessen.c after call to hhvec in hsehldr.c
/* 	  to avoid infinite loop in symmeig()
/* 	- hqp/Hqp_HL_BFGS.C: change default _eigen_control to true
/* 	- omu/Omu_IntDASPK.C (realloc): more sensible allocation of memory
/* 	  when using direct solver
/*
/* Revision 1.2  2000/04/11 18:45:23  rf
/* removed nested comments at end of file
/*
 * Revision 1.1.1.1  2000/04/04 07:25:42  ea
 * CVS initial version 1.6a4
 *
 * Revision 1.9  1994/03/24  00:04:05  des
 * Added notes on changes to spILUfactor() and px_inv().
 *
 * Revision 1.8  1994/02/21  04:32:25  des
 * Set version to 1.2b with bug fixes in schur() and installation.
 *
 * Revision 1.7  1994/01/13  05:43:57  des
 * Version 1.2 update
 *

 * */

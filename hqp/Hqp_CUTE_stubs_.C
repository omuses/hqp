/*
 * Hqp_CUTE_stubs_.C -- 
 *   stub procedures for CUTE,
 *   needed as some Fortran compilers add underscore to symbol names
 *
 * rf, 3/19/97
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

#include <stdio.h>
#include <stdlib.h>

#ifndef fint
#define fint int
#endif

#ifndef freal
#define freal double
#endif

#ifndef fbool
#define fbool int
#endif

extern "C" void 
csize(fint *N, fint *M, fint *NNZJ, fint *NNZH);

extern "C" void 
csize_(fint *N, fint *M, fint *NNZJ, fint *NNZH)
{
  csize(N, M, NNZJ, NNZH);
}

extern "C" void
cinit(fint *N, fint *M, freal *X0, freal *BL, freal *BU, freal *INF,
      fbool *EQUATN, fbool *LINEAR, freal *V0, freal *CL, freal *CU,
      fbool *EFIRST, fbool *LFIRST, fbool *NVFRST);

extern "C" void
cinit_(fint *N, fint *M, freal *X0, freal *BL, freal *BU, freal *INF,
       fbool *EQUATN, fbool *LINEAR, freal *V0, freal *CL, freal *CU,
       fbool *EFIRST, fbool *LFIRST, fbool *NVFRST)
{
  cinit(N, M, X0, BL, BU, INF,
        EQUATN, LINEAR, V0, CL, CU,
        EFIRST, LFIRST, NVFRST);
}

extern "C" void
cfn(const fint *N, const fint *M, const freal *X,
    freal *F, const fint *LC, freal *C);

extern "C" void
cfn_(const fint *N, const fint *M, const freal *X,
     freal *F, const fint *LC, freal *C)
{
  cfn(N, M, X,
      F, LC, C);
}

extern "C" void
csgr(const fint *N, const fint *M, const fbool *GRLAGF,
     const fint *LV, const freal *V, const freal *X,
     fint *NNZSCJ, const fint *LSCJAC, freal *SCJAC,
     fint *INDVAR, fint *INDFUN);

extern "C" void
csgr_(const fint *N, const fint *M, const fbool *GRLAGF,
      const fint *LV, const freal *V, const freal *X,
      fint *NNZSCJ, const fint *LSCJAC, freal *SCJAC,
      fint *INDVAR, fint *INDFUN)
{
  csgr(N, M, GRLAGF,
       LV, V, X,
       NNZSCJ, LSCJAC, SCJAC,
       INDVAR, INDFUN);
}

extern "C" void
cscifg(const fint *N, const fint *I, const freal *X,
       freal *CI, fint *NNZSGC, const fint *LSGCI, freal *SGCI,
       fint *IVSGCI, fbool *GRAD);

extern "C" void
cscifg(const fint *N, const fint *I, const freal *X,
       freal *CI, fint *NNZSGC, const fint *LSGCI, freal *SGCI,
       fint *IVSGCI, fbool *GRAD);

extern "C" void
cscifg_(const fint *N, const fint *I, const freal *X,
	freal *CI, fint *NNZSGC, const fint *LSGCI, freal *SGCI,
	fint *IVSGCI, fbool *GRAD)
{
  cscifg(N, I, X,
	 CI, NNZSGC, LSGCI, SGCI,
	 IVSGCI, GRAD);
}

extern "C" void
csgrsh(const fint *N, const fint *M, const freal *X, const fbool *GRLAGF,
       const fint *LV, const freal *V,
       fint *NNZSCJ, const fint *LSCJAC, freal *SCJAC,
       fint *INDVAR, fint *INDFUN,
       fint *NNZSH, const fint *LSH, freal *SH,
       fint *IRNSH, fint *ICNSH);

extern "C" void
csgrsh_(const fint *N, const fint *M, const freal *X, const fbool *GRLAGF,
	const fint *LV, const freal *V,
	fint *NNZSCJ, const fint *LSCJAC, freal *SCJAC,
	fint *INDVAR, fint *INDFUN,
	fint *NNZSH, const fint *LSH, freal *SH,
	fint *IRNSH, fint *ICNSH)
{
  csgrsh(N, M, X, GRLAGF,
	 LV, V,
	 NNZSCJ, LSCJAC, SCJAC,
	 INDVAR, INDFUN,
	 NNZSH, LSH, SH,
	 IRNSH, ICNSH);
}

extern "C" void
csgreh(const fint *N, const fint *M, const freal *X, const fbool *GRLAGF,
       const fint *LV, const freal *V,
       fint *NNZSCJ, const fint *LSCJAC, freal *SCJAC,
       fint *INDVAR, fint *INDFUN, fint *NE, fint *IRNHI,
       const fint *LIRNHI, const fint *LE, fint *IPRNHI,
       freal *HI, const fint *LHI, fint *IPRHI, const fbool *BYROWS);

extern "C" void
csgreh_(const fint *N, const fint *M, const freal *X, const fbool *GRLAGF,
	const fint *LV, const freal *V,
	fint *NNZSCJ, const fint *LSCJAC, freal *SCJAC,
	fint *INDVAR, fint *INDFUN, fint *NE, fint *IRNHI,
	const fint *LIRNHI, const fint *LE, fint *IPRNHI,
	freal *HI, const fint *LHI, fint *IPRHI, const fbool *BYROWS)
{
  csgreh(N, M, X, GRLAGF,
	 LV, V,
	 NNZSCJ, LSCJAC, SCJAC,
	 INDVAR, INDFUN, NE, IRNHI,
	 LIRNHI, LE, IPRNHI,
	 HI, LHI, IPRHI, BYROWS);
}

extern "C" void
cwrtsn(const fint *N, const fint *M, const char *header,
       const freal *F, const freal *X, const freal *V);

extern "C" void
cwrtsn_(const fint *N, const fint *M, const char *header,
	const freal *F, const freal *X, const freal *V)
{
  cwrtsn(N, M, header,
	 F, X, V);
}


//=========================================================================

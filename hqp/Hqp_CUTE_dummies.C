/*
 * Hqp_CUTE_dummies.C -- 
 *   dummy procedures used if CUTE is not linked
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
csize(fint *N, fint *M, fint *NNZJ, fint *NNZH)
{
  fprintf(stderr, "Fatal: CUTE procedure CSIZE called but not linked!\n"); 
  exit(-1);
}

extern "C" void
cinit(fint *N, fint *M, freal *X0, freal *BL, freal *BU, freal *INF,
       fbool *EQUATN, fbool *LINEAR, freal *V0, freal *CL, freal *CU,
       fbool *EFIRST, fbool *LFIRST, fbool *NVFRST)
{
  fprintf(stderr, "Fatal: CUTE procedure CINIT called but not linked!\n"); 
  exit(-1);
}

extern "C" void
cfn(const fint *N, const fint *M, const freal *X,
    freal *F, const fint *LC, freal *C)
{
  fprintf(stderr, "Fatal: CUTE procedure CFN called but not linked!\n"); 
  exit(-1);
}

extern "C" void
csgr(const fint *N, const fint *M, const fbool *GRLAGF,
     const fint *LV, const freal *V, const freal *X,
     fint *NNZSCJ, const fint *LSCJAC, freal *SCJAC,
     fint *INDVAR, fint *INDFUN)
{
  fprintf(stderr, "Fatal: CUTE procedure CSGR called but not linked!\n"); 
  exit(-1);
}

extern "C" void
cscifg(const fint *N, const fint *I, const freal *X,
       freal *CI, fint *NNZSGC, const fint *LSGCI, freal *SGCI,
       fint *IVSGCI, fbool *GRAD)
{
  fprintf(stderr, "Fatal: CUTE procedure CSCIFG called but not linked!\n"); 
  exit(-1);
}

extern "C" void
csgrsh(const fint *N, const fint *M, const freal *X, const fbool *GRLAGF,
       const fint *LV, const freal *V,
       fint *NNZSCJ, const fint *LSCJAC, freal *SCJAC,
       fint *INDVAR, fint *INDFUN,
       fint *NNZSH, const fint *LSH, freal *SH,
       fint *IRNSH, fint *ICNSH)
{
  fprintf(stderr, "Fatal: CUTE procedure CSGRSH called but not linked!\n"); 
  exit(-1);
}

extern "C" void
csgreh(const fint *N, const fint *M, const freal *X, const fbool *GRLAGF,
       const fint *LV, const freal *V,
       fint *NNZSCJ, const fint *LSCJAC, freal *SCJAC,
       fint *INDVAR, fint *INDFUN, fint *NE, fint *IRNHI,
       const fint *LIRNHI, const fint *LE, fint *IPRNHI,
       freal *HI, const fint *LHI, fint *IPRHI, const fbool *BYROWS)
{
  fprintf(stderr, "Fatal: CUTE procedure CSGREH called but not linked!\n"); 
  exit(-1);
}

extern "C" void
cwrtsn(const fint *N, const fint *M, const char *header,
       const freal *F, const freal *X, const freal *V)
{
  fprintf(stderr, "Fatal: CUTE procedure CWRTSN called but not linked!\n"); 
  exit(-1);
}


//=========================================================================

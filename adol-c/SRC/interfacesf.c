#define _INTERFACESF_C_
#define _ADOLC_SRC_
/*
   ------------------------------------------------------------- 
   File interfacesf.c of ADOL-C version 1.8.0    as of Nov/30/98
   ------------------------------------------------------------- 
   Genuine Fortran callable C Interfaces to ADOL-C forward 
   & reverse calls.

   Last changes:
      981130 olvo:    newly created from driversc.c

   History of driversc.c:
      981126 olvo:    last check (p's & q's)

   --------------------------------------------------------------
*/


/****************************************************************************/
/*                                                                 INCLUDES */
#include "dvlparms.h"
#include "usrparms.h"
#include "interfaces.h"
#include "adalloc.h"
#include "fortutils.h"

#ifdef __cplusplus
extern "C" {
#endif


/*--------------------------------------------------------------------------*/
fint hos_forward_(fint* ftag,
                  fint* fm,
                  fint* fn,
                  fint* fd,
                  fint* fk,
                  fdouble* fbase,
                  fdouble* fx,
                  fdouble* fvalue,
                  fdouble* fy)
{
  int rc= -1;
  int tag=*ftag, m=*fm, n=*fn, d=*fd, k=*fk;
  double* base = myalloc1(n);
  double* value = myalloc1(m);
  double** X = myalloc2(n,d);
  double** Y = myalloc2(m,d);
  spread1(n,fbase,base);
  spread2(n,d,fx,X);
  rc= hos_forward(tag,m,n,d,k,base,X,value,Y);
  pack2(m,d,Y,fy);
  pack1(m,value,fvalue);
  free((char*)*X);free((char*)X);
  free((char*)*Y);free((char*)Y);
  free((char*)base); free((char*)value);
  return rc;
}

/*--------------------------------------------------------------------------*/
fint zos_forward_(fint* ftag,
                  fint* fm,
                  fint* fn,
                  fint* fk,
                  fdouble* fbase,
                  fdouble* fvalue)
{
  int rc=-1;
  int tag=*ftag, m=*fm, n=*fn, k=*fk;
  double* base=myalloc1(n);
  double* value = myalloc1(m);
  spread1(n,fbase,base);
  rc=zos_forward(tag,m,n,k,base,value);
  pack1(m,value,fvalue);
  free((char*)base); free((char*)value);
  return rc;
}

/*--------------------------------------------------------------------------*/
fint hov_forward_(fint* ftag,
                  fint* fm,
                  fint* fn,
                  fint* fd,
                  fint* fp,
                  fdouble* fbase,
                  fdouble* fx,
                  fdouble* fvalue,
                  fdouble* fy)
{
  int rc= -1;
  int tag=*ftag, m=*fm, n=*fn, d=*fd, p=*fp;
  double* base = myalloc1(n);
  double* value = myalloc1(m);
  double*** X = myalloc3(n,p,d);
  double*** Y = myalloc3(m,p,d);
  spread1(n,fbase,base);
  spread3(n,p,d,fx,X);
  rc= hov_forward(tag,m,n,d,p,base,X,value,Y);
  pack3(m,p,d,Y,fy);
  pack1(m,value,fvalue);
  free((char*)**X); free((char*)*X); free((char*)X);
  free((char*)**Y); free((char*)*Y); free((char*)Y);
  free((char*)base); free((char*)value);
  return rc;
}

/*--------------------------------------------------------------------------*/
fint fov_forward_(fint* ftag,
                  fint* fm,
                  fint* fn,
                  fint* fp,
                  fdouble* fbase,
                  fdouble* fx,
                  fdouble* fvalue,
                  fdouble* fy)
{
  int rc= -1;
  int tag=*ftag, m=*fm, n=*fn, p=*fp;
  double* base = myalloc1(n);
  double* value = myalloc1(m);
  double** X = myalloc2(n,p);
  double** Y = myalloc2(m,p);
  spread1(n,fbase,base);
  spread2(n,p,fx,X);
  rc= fov_forward(tag,m,n,p,base,X,value,Y);
  pack2(m,p,Y,fy);
  pack1(m,value,fvalue);
  free((char*)*X); free((char*)X);
  free((char*)*Y); free((char*)Y);
  free((char*)base); free((char*)value);
  return rc;
}

  
/*--------------------------------------------------------------------------*/
fint hos_reverse_(fint* ftag,
       		  fint* fm,
		  fint* fn,
		  fint* fd,
		  fdouble* fu,
		  fdouble* fz)
{
  int rc=-1;
  int tag=*ftag, m=*fm, n=*fn, d=*fd;
  double** Z = myalloc2(n,d+1);
  double* u = myalloc1(m);
  spread1(m,fu,u);
  rc=hos_reverse(tag,m,n,d,u,Z);
  pack2(n,d+1,Z,fz);
  free((char*)*Z); free((char*)Z);
  free((char*)u);
  return rc;
}

/*--------------------------------------------------------------------------*/
fint fos_reverse_(fint* ftag,
		  fint* fm,
		  fint* fn,
		  fdouble* fu,
		  fdouble* fz)
{
  int rc=-1;
  int tag=*ftag, m=*fm, n=*fn;
  double* u = myalloc1(m);
  double* Z = myalloc1(n);
  spread1(m,fu,u);
  rc=fos_reverse(tag,m,n,u,Z);
  pack1(n,Z,fz);
  free((char*)Z); free((char*)u);
  return rc;
}

/*--------------------------------------------------------------------------*/
fint hov_reverse_(fint* ftag,
		  fint* fm,
		  fint* fn,
		  fint* fd,
		  fint* fq,
		  fdouble* fu,
		  fdouble* fz)
{
  int rc=-1;
  int tag=*ftag, m=*fm, n=*fn, d=*fd, q=*fq;
  double** U = myalloc2(q,m);
  double*** Z = myalloc3(q,n,d+1);
  short ** nop = 0;
  spread2(q,m,fu,U);
  rc=hov_reverse(tag,m,n,d,q,U,Z,nop);
  pack3(q,n,d+1,Z,fz);
  free((char*)**Z); free((char*)*Z); free((char*)Z);
  free((char*)*U); free((char*)U);
  return rc;
}

/*--------------------------------------------------------------------------*/
fint fov_reverse_(fint* ftag,
		  fint* fm,
		  fint* fn,
		  fint* fq,
		  fdouble* fu,
		  fdouble* fz)
{
  int rc=-1;
  int tag=*ftag, m=*fm, n=*fn, q=*fq;
  double** U = myalloc2(q,m);
  double** Z = myalloc2(q,n);
  spread2(q,m,fu,U);
  rc=fov_reverse(tag,m,n,q,U,Z);
  pack2(q,n,Z,fz);
  free((char*)*Z); free((char*)Z);
  free((char*)*U); free((char*)U);
  return rc;
}


/****************************************************************************/
/*                                                               THAT'S ALL */
#ifdef __cplusplus
}
#endif

#undef _ADOLC_SRC_
#undef _INTERFACESF_C_

#define _DRIVERSC_C_
#define _ADOLC_SRC_
/*
  ------------------------------------------------------------------------
  File driversc.c of ADOL-C version 1.8.7                as of Feb/28/2000
  ------------------------------------------------------------------------
   Easy to use drivers for optimization and nonlinear equations
   (Implementation of the C/C++ callable interfaces).

  Last changed : 
  20000228 olvo:    corrected comment at lagra_hess_vec
    990622 olvo:    jacobian(..) makes decision to use forward only or
                    a reverse sweep
    981130 olvo:    newly created from old wersion (which was splitted)
    981126 olvo:    last check (p's & q's)
    981020 olvo:    deleted debug messages in tensor

  ------------------------------------------------------------------------
*/


/****************************************************************************/
/*                                                                 INCLUDES */
#include "dvlparms.h" /* Developers Parameters */
#include "usrparms.h" /* Users Parameters      */
#include "drivers.h"
#include "interfaces.h"
#include "adalloc.h"

#include <stdio.h>
#include <math.h>

#ifdef __cplusplus
#include <malloc.h>
extern "C" {
#endif


/****************************************************************************/
/*                                                                   MACROS */

#define maxinc(a,b) if ((a) < (b)) (a) = (b)
#define mindec(a,b) if ((a) > (b)) (a) = (b)


/****************************************************************************/
/*                         DRIVERS FOR OPTIMIZATION AND NONLINEAR EQUATIONS */

/*--------------------------------------------------------------------------*/
/*                                                                 function */
/* function(tag, m, n, x[n], y[m])                                          */
                                                                          
int function(short tag,
             int m,
	     int n,
	     double* argument,
	     double* result)
{
  int rc= -1;
  
  rc= zos_forward(tag,m,n,0,argument,result);

  return rc;
}

/*--------------------------------------------------------------------------*/
/*                                                                 gradient */
/* gradient(tag, n, x[n], g[n])                                             */

static double one = 1.0;

int gradient(short tag,
	     int n,
	     double* argument,
	     double* result)
{
  int rc= -1;
  
  rc = zos_forward(tag,1,n,1,argument,result);
  if(rc < 0)
    return rc;
  maxinc(rc, fos_reverse(tag,1,n,&one,result));
  return rc;
}


/*--------------------------------------------------------------------------*/
/*                                                                          */
/* vec_jac(tag, m, n, repeat, x[n], u[m], v[n])                             */

int vec_jac(short tag,
	    int m,
	    int n,
	    int repeat,
            double* argument,
            double* lagrange,
            double* row)
{
  int rc= -1;
  static double *y;
  static int maxm;
  if(m > maxm)
    {
      if(maxm)
        free((char*)y);
      y = myalloc1(m);
      maxm = m;
    } 
  if(!repeat) {
    rc = zos_forward(tag,m,n,1, argument, y);
    if(rc < 0)
     return rc;
  } 
  maxinc(rc, fos_reverse(tag,m,n,lagrange,row)); 
  return rc;
}


/*--------------------------------------------------------------------------*/
/*                                                                 jacobian */
/* jacobian(tag, m, n, x[n], J[m][n])                                       */

int jacobian(short tag,
	     int depen,
	     int indep,
	     double *argument,
	     double **jacobian)
{ int rc;
  static int nmmax; // dim for I
  static int mmax;  // dim for result
  static double *result, **I; 

  if (depen > mmax)
  { if (mmax)
      myfree1(result);
    result = myalloc1(mmax = depen); 
  }

  if (indep/2 < depen)
  { if (indep > nmmax)   
    { if (nmmax) 
        myfreeI2(nmmax,I);
      I = myallocI2(nmmax = indep);
    }
    rc = fov_forward(tag,depen,indep,indep,argument,I,result,jacobian); 
  }
  else
  { if (depen > nmmax)   
    { if (nmmax) 
        myfreeI2(nmmax,I);
      I = myallocI2(nmmax = depen);
    }
    rc = zos_forward(tag,depen,indep,1,argument,result);
    if (rc < 0) 
      return rc;
    maxinc(rc,fov_reverse(tag,depen,indep,depen,I,jacobian));
  }
  return rc;
}


/*--------------------------------------------------------------------------*/
/*                                                                          */
/* jac_vec(tag, m, n, x[n], v[n], u[m]);                                    */

int jac_vec(short tag,
            int m, 
            int n,
            double* argument,
            double* tangent,
            double* column)
{
  int rc= -1;
  static double *y;
  static int maxm;
  if(m > maxm )
    {
      if(maxm)
        myfree1(y);
      y = myalloc1(m);
      maxm = m;     
    } 
  
  rc = fos_forward(tag, m, n, 0, argument, tangent, y, column);
  
  return rc;
}


/*--------------------------------------------------------------------------*/
/*                                                                          */
/* hess_vec(tag, n, x[n], v[n], w[n])                                       */

int hess_vec(short tag,
	     int n,
	     double *argument,
	     double *tangent,
	     double *result)
{
  return lagra_hess_vec(tag,1,n,argument,tangent,&one,result);
}


/*--------------------------------------------------------------------------*/
/*                                                                  hessian */
/* hessian(tag, n, x[n], lower triangle of H[n][n])                         */

int hessian(short tag,
	    int n,
	    double* argument,
	    double** hess)
{
  int rc= 3;
  int i,j;
  double *v = myalloc1(n);
  double *w = myalloc1(n);
  for(i=0;i<n;i++)
    v[i] = 0;
  for(i=0;i<n;i++)
    {
      v[i] = 1;   
      mindec(rc,hess_vec(tag,n,argument,v,w));
      if( rc < 0){
        free((char *)v); 
        free((char *) w); 
        return rc;
      }
      for(j=0;j<=i;j++)
	hess[i][j] = w[j];
      v[i] = 0;
    }

  free((char *)v); 
  free((char *) w); 
  return rc;
  /* Note that only the lower triangle of hess is filled */
}


/*--------------------------------------------------------------------------*/
/*                                                                          */
/* lagra_hess_vec(tag, m, n, x[n], v[n], u[m], w[n])                        */

int lagra_hess_vec(short tag,
	           int m,
	           int n,
	           double *argument,
	           double *tangent,
	           double *lagrange,
	           double *result)
{ 
  int rc=-1;
  int i;
  int degree = 1;
  int keep = degree+1;
  static double **X, *y, *y_tangent;
  static int maxn, maxm;
  
  if (n > maxn || m > maxm)
    {
      if (X){  free((char*)*X); free((char*)X);} 
      X = myalloc2(n,2);
      maxn = n;
      if (y) free((char*)y);
      if (y_tangent) free((char*)y_tangent);
      y         = myalloc1(m);
      y_tangent = myalloc1(m);
    }
 
  rc = fos_forward(tag,m,n,keep, argument, tangent, y, y_tangent);

  if(rc < 0) 
     return rc;
  
  for(i=0;i<n;i++)
    {
      X[i][0] = argument[i];
      X[i][1] = tangent[i];
    }
   
  maxinc(rc,hos_reverse(tag,m,n,degree,lagrange,X));
  
  for(i=0;i<n;i++)
    result[i] = X[i][1];
  
  return rc;
}


/****************************************************************************/
/*                                                       OLD TENSOR DIRVERS */

/*--------------------------------------------------------------------------*/
void swap( int* i, int* j)
{ int si; if(*i < *j) {si = *i; *i= *j; *j = si;}; }

int hessloc(int i, int j)
{
/* swap(i,j); */
return i*(i+1)/2 + j;
}

int tensloc(int i, int j, int k)
{
/* swap(i,j); swap(j,k); swap(i,j); */
return (i*(i+1)*(i+2)/6 + j*(j+1)/2+k);
}

int tensor( int tag, 
	    int m,
	    int n,
	    double* argument,
	    double* functions,
	    double** gradients,
	    double** hessians,     /* in symmetric storage mode */
            double** tensors       /* in symmetric storage mode */
	  )
{ 
int rc=-1;
static int nt,i,ii,j,jj,k,l,li,lj,lij,lji,lli,lk,llii,
	llij,mm,nn,lijk,ljk,lik,lkj,lki;
double ***X, ***Y, temp;
nt = tensloc(n-1,n-1,n-1)+1;
X = myalloc3(n,nt,3);
Y = myalloc3(m,nt,3);
l = 0;
/* Seed the input Taylor Series with the directions (e_i+e_j+e_k)/3 */ 
for(i=0;i<n;i++)
  for(j=0;j<=i;j++)
    for(k=0;k<=j;k++)
      {
      for(nn=0;nn<n;nn++) 
	 for(jj=0;jj<3;jj++)
	    X[nn][l][jj] = 0;
      *X[i][l] += 1.0/3;
      *X[j][l] += 1.0/3;
      *X[k][l] += 1.0/3;
      l++;
      }
/* Propagate the Taylor series in the forward mode */
rc=hov_forward(tag,m,n,3,l,argument,X,functions,Y);
for(i=0;i<n;i++)
   {
   li = tensloc(i,i,i);
   lli = hessloc(i,i);
/* Copy over the pure derivatives along the axes */
   for(mm=0;mm<m;mm++)
     {
     gradients[mm][i] = Y[mm][li][0];
     hessians[mm][lli] = 2*Y[mm][li][1];
     tensors[mm][li] = 6*Y[mm][li][2];
     }
   lij = tensloc(i,i,0);
   llij = hessloc(i,0);
   for(j=0;j<i;j++) 
      {
      lj = tensloc(j,j,j);
      lji = tensloc(i,j,j);
/* Look at each coordinate plane and interpolate mixed derivatives */
      for(mm=0;mm<m;mm++)
      {
        hessians[mm][llij] = (9*(Y[mm][lij][1]+Y[mm][lji][1])
                              -5*(Y[mm][li][1]+Y[mm][lj][1]))/4;
        tensors[mm][lij] = 9*(2*Y[mm][lij][2]-Y[mm][lji][2])
			     +2*Y[mm][lj][2]-5*Y[mm][li][2];
        tensors[mm][lji] = 9*(2*Y[mm][lji][2]-Y[mm][lij][2])
			     +2*Y[mm][li][2]-5*Y[mm][lj][2];
	*Y[mm][lij] = Y[mm][li][2]+Y[mm][lj][2]
		      - 4.5*(Y[mm][lij][2]+Y[mm][lji][2]);
       }
	lijk = tensloc(i,j,0);
	lik = tensloc(i,i,0);
	ljk = tensloc(j,j,0);
        for(k=0;k<j;k++)
	 {
/* Look at each coordinate triple and interpolate mixed derivative */
          for(mm=0;mm<m;mm++)
	   tensors[mm][lijk] = 27*Y[mm][lijk][2] + *Y[mm][lij]
		             + *Y[mm][lik] + *Y[mm][ljk];
          lijk++; lik++; ljk ++;
       }
       lij++; llij++;
    }
  }
  free((char*)**Y); free((char*)*Y); free((char*)Y);
  free((char*)**X); free((char*)*X); free((char*)X);
  return rc;
}


/****************************************************************************/
/*                                                               THAT'S ALL */

#ifdef __cplusplus
}
#endif

#undef _ADOLC_SRC_
#undef _DRIVERSC_C_

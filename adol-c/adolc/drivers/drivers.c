/*---------------------------------------------------------------------------- 
 ADOL-C -- Automatic Differentiation by Overloading in C++
 File:     drivers/drivers.c
 Revision: $Id: drivers.c,v 1.1 2004/10/13 14:18:12 e_arnold Exp $
 Contents: Easy to use drivers for optimization and nonlinear equations
           (Implementation of the C/C++ callable interfaces).

 Copyright (c) 2003
               Technical University Dresden
               Department of Mathematics
               Institute of Scientific Computing
  
 This file is part of ADOL-C. This software is provided under the terms of
 the Common Public License. Any use, reproduction, or distribution of the
 software constitutes recipient's acceptance of the terms of this license.
 See the accompanying copy of the Common Public License for more details.

 History:
          20030306 andrea:  change maxinc to mindec for rc computation
          20030303 andrea:  new hess_mat(..), new hessian2(..) 
          19990622 olvo:    jacobian(..) makes decision whether to use the
                            forward or reverse mode
----------------------------------------------------------------------------*/
#include "../drivers/drivers.h"
#include "../interfaces.h"
#include "../adalloc.h"

#include <math.h>
#include <malloc.h>

BEGIN_C_DECLS

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
  mindec(rc, fos_reverse(tag,1,n,&one,result));
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
  mindec(rc, fos_reverse(tag,m,n,lagrange,row)); 
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
  static int nmmax; /* dim for I */
  static int mmax;  /* dim for result */
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
    mindec(rc,fov_reverse(tag,depen,indep,depen,I,jacobian));
  }
  return rc;
}

/*--------------------------------------------------------------------------*/
/*                                                           jacobian_partx */
/* jacobian_partx(tag, m, n, x[n][], J[m][n][])                             */

int jacobian_partx(short tag,
           int depen,
           int n,
             int *ndim,
           double **x,
           double ***J)
{ int rc;
  double *argument, **jac; 
  int i,j,k,ind,indep;

  indep = 0;
  for(i=0;i<n;i++)
    indep += ndim[i];
  
  argument = myalloc1(indep);
  jac = myalloc2(depen,indep);
  ind = 0;
  for(i=0;i<n;i++)
    for(j=0;j<ndim[i];j++)
     {
      argument[ind] = x[i][j];
      ind++;
     }

  rc = jacobian(tag,depen,indep,argument,jac);

  for(i=0;i<depen;i++)
   {
    ind = 0;
    for(j=0;j<n;j++)
      for(k=0;k<ndim[j];k++)
       {
        J[i][j][k] = jac[i][ind];
        ind++;
       }
   }

  return rc;
}

/*--------------------------------------------------------------------------*/
/*                                                                  jac_vec */
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
/*                                                                 hess_vec */
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
/*                                                                 hess_mat */
/* hess_mat(tag, n, q, x[n], V[n][q], W[q][n])                              */
int hess_mat(short tag,
	     int n,
	     int q,
	     double *argument,
	     double **tangent,
	     double **result)
{
  int rc;
  int i,j;

  double*** Xppp = myalloc3(n,q,1);   /* matrix on right-hand side  */
  double*   y    = myalloc1(1);       /* results of function evaluation */  
  double*** Yppp = myalloc3(1,q,1);   /* results of hos_wk_forward  */  
  double*** Zppp = myalloc3(q,n,2);   /* result of Up x H x XPPP */            
  double**  Upp  = myalloc2(1,2);     /* vector on left-hand side */
 
  for (i=0; i<n; i++)
    for (j=0;j<q;j++)
       Xppp[i][j][0] = tangent[i][j]; 

  Upp[0][0] = 1; Upp[0][1] = 0;

  rc = hov_wk_forward(tag,1,n,1,2,q,argument,Xppp,y,Yppp);   
  mindec(rc,hos_ov_reverse(tag,1,n,1,q,Upp,Zppp)); 
  
  for (i=0; i<q; i++)
    for (j=0;j<n;j++)
      result[i][j] = Zppp[i][j][1]; 

  myfree2(Upp);
  myfree3(Zppp);
  myfree3(Yppp);
  myfree1(y);
  myfree3(Xppp);
  return rc;
}

/*--------------------------------------------------------------------------*/
/*                                                                  hessian */
/* hessian(tag, n, x[n], lower triangle of H[n][n])                         */
/* uses Hessian-vector product                                              */
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
/*                                                                 hessian2 */
/* hessian2(tag, n, x[n], lower triangle of H[n][n])                        */
/* uses Hessian-matrix product                                              */
int hessian2(short tag,
	    int n,
	    double* argument,
	    double** hess)
{
  int rc;
  int i,j;

  double*** Xppp = myalloc3(n,n,1);   /* matrix on right-hand side  */
  double*   y    = myalloc1(1);       /* results of function evaluation */  
  double*** Yppp = myalloc3(1,n,1);   /* results of hos_wk_forward  */  
  double*** Zppp = myalloc3(n,n,2);   /* result of Up x H x XPPP */            
  double**  Upp  = myalloc2(1,2);     /* vector on left-hand side */
 
  for (i=0; i<n; i++)
   {
    for (j=0;j<n;j++)
       Xppp[i][j][0] = 0;
    Xppp[i][i][0] = 1;
   }

  Upp[0][0] = 1; Upp[0][1] = 0;

  rc = hov_wk_forward(tag,1,n,1,2,n,argument,Xppp,y,Yppp);   
  mindec(rc,hos_ov_reverse(tag,1,n,1,n,Upp,Zppp)); 
  
  for (i=0; i<n; i++)
    for (j=0;j<=i;j++)
      hess[i][j] = Zppp[i][j][1]; 

  myfree2(Upp);
  myfree3(Zppp);
  myfree3(Yppp);
  myfree1(y);
  myfree3(Xppp);
  return rc; 
  /* Note that only the lower triangle of hess is filled */
}

/*--------------------------------------------------------------------------*/
/*                                                           lagra_hess_vec */
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
      maxn = n; maxm=m; /* ov20020116 set maxm to m */
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
   
  mindec(rc,hos_reverse(tag,m,n,degree,lagrange,X));
  
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
static int nt,i,j,jj,k,l,li,lj,lij,lji,lli,llij,mm,nn,lijk,ljk,lik;
double ***X, ***Y;
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

END_C_DECLS

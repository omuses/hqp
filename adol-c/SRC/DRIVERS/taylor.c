#define _TAYLORC_C_
#define _ADOLC_SRC_
/*
   ------------------------------------------------------------- 
   File taylor.c of ADOL-C version 1.8.1         as of Feb/16/99
   ------------------------------------------------------------- 

   Last changes:
      990216  olvo/andrea:   mindec in tensor_eval geklammert
      981214  olvo/walther:  ec Xhelp
      981130  olvo:     last check (includes ...)
      981120  olvo/walther:  return values
      980914  olvo:          Jac_solv --> jac_solv
      980814  olvo:          integral constant expressions in
                             arrays (ANSI-C) 
      980806  walther:       (1) access to tensors 
                             (2) seems to be finished 
      980804  olvo/walther:  debugging --> finished? 

   ------------------------------------------------------------- 
*/

/****************************************************************************/
/*                                                                 INCLUDES */
#include "dvlparms.h"
#include "usrparms.h"
#include "taylor.h"
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
#define compsize >
#define mindec(a,b) if ((a) > (b)) (a) = (b)


/****************************************************************************/
/*                                                              STRUCT ITEM */
struct item {
  int a;                 /* address in array of derivatives */ 
  int b;                 /* absolute value of the correspondig multiindex i */
  double c;              /* value of the coefficient c_{i,j} */
  struct item *next;     /* next item */
 };


/****************************************************************************/
/*                                                     DEALLOCATE COEFFLIST */
void freecoefflist( int dim, struct item *coeff_list )
{ int i;
  struct item *ptr1;
  struct item *ptr2;

  for (i=0; i<dim; i++)  /* sum over all multiindices jm with |jm| = d */
  { ptr1 = &coeff_list[i];
    ptr1 = ptr1->next;
    while (ptr1 != NULL) 
    { ptr2 = ptr1->next;
      free((char *) ptr1);    
      ptr1 = ptr2;
    }
  }
}


/****************************************************************************/
/*                                            ALLOCATE/DEALLOCATE COEFFLIST */
double* tensoriglob;

/*--------------------------------------------------------------------------*/
/* Allcoate space for symmetric derivative tensors
   of up to order d in n variables, derivatives are  */
void* tensorpoint( int n, int d )
{ int i;
  void* t;

  if (d == 1)
  { t = (void*) tensoriglob;
    tensoriglob += n+1;
  }
  else
  { t = (void*) malloc((n+1)*sizeof(void*));
    for (i=0; i<=n; i++)
      ((void**)t)[i] = (void*) tensorpoint(i,d-1);
  }
  return t;
}

/*--------------------------------------------------------------------------*/
void** tensorsetup( int m, int n, int d, double** tensorig )
{ int i;
  void** t = (void**) malloc(m*sizeof(void*));

  for (i=0; i<m; i++)
  { tensoriglob = tensorig[i];
    t[i] = (void*) tensorpoint(n,d);
  }
  return t;
}

/*--------------------------------------------------------------------------*/
/* Deallocate space for symmetric derivative tensors
   of up to order d in n variables  */
void freetensorpoint( int n, int d, double** tensor )
{ int i;
  double* t;

  if (d > 2)
    for(i=0;i<=n;i++)
    { t = tensor[i];
      freetensorpoint(i,d-1,(double **) t);
      free((char *) t);
    }
}

/*--------------------------------------------------------------------------*/
void freetensor( int m, int n, int d, double** tensor )
{ int i;
  double* t;

  for (i=0; i<m; i++)
  { t = tensor[i];
    freetensorpoint(n,d,(double **)t);
    free((char *) t);
  }
}


/****************************************************************************/
/*                                                           SOME UTILITIES */

/*--------------------------------------------------------------------------*/
int binomi( int a, int b )
{ int i;
  int result = 1;

  for (i=1; i<=b; i++)
    result = (result*(a-i+1))/i; /* olvo 980804 () */ 
  return result;
 }

/*--------------------------------------------------------------------------*/
double dbinomi( double a, int b )
{ int i;
  double result = 1.0;

  for (i=1; i<=b; i++)
    result = result*(a-i+1)/i;
  return result;
}

/*--------------------------------------------------------------------------*/
void convert( int p, int d, int *im, int *multi )
{ int i;
 
  for (i=0; i<p; i++)
    multi[i] = 0;
  for (i=0; i<d; i++)
    if (im[i]) /* olvo/walther 980804 new tl */
      multi[im[i]-1] += 1;
}

/*--------------------------------------------------------------------------*/
int address( int d, int* im )
{ int i, 
      add = 0;

  for (i=0; i<d; i++)
    add += binomi(im[i]+i,i+1);
  return add;
}

/*--------------------------------------------------------------------------*/
void multiind( int p, int m, int d, int **multi )   
{ int i,j,num_jm,rem,add;
  int *jm_it = (int*) malloc(m*sizeof(int));
  int *im    = (int*) malloc(d*sizeof(int));
  
  for (i=0; i<m-1; i++)
    jm_it[i] = 1;
  jm_it[m-1] = 0;
  
  num_jm = binomi(p+m-1,m);
  for (i=0; i<num_jm; i++)
  { for (j=0; j<d-m; j++)
      im[j]=0;
    jm_it[m-1] = jm_it[m-1] +1;
    for (j=m-2; j>=0; j--)
      jm_it[j] = jm_it[j] + jm_it[j+1]/(p+1);
    im[d-m] = jm_it[0];
    for (j=1; j<m; j++)
    { if (jm_it[j] > p) 
        jm_it[j] = jm_it[j-1];
      im[d-m+j] = jm_it[j];
    }
    add = address(d,im);
    convert(p,d,im,multi[i]);
    multi[i][p]=add;
  }
  free((char*) im);
  free((char*) jm_it);
}


/****************************************************************************/
/*                                                    EVALUATE COEFFICIENTS */
void coeff( int p, int d, struct item *coeff_list ) 
{ int m,i,ii,j,k,kk,l,dim,dimi,diml,q;
  int    ***multi;   /* multiindex */  
  double d_real;
  double **bin_coeff;
  double *factor; 
  double addend, rem, coeffpart1, coeffpart2;
  struct item **ptr;

  dim  = binomi(p+d-1,d);
  dimi = binomi(p+d,d);
  factor    = (double *)       malloc(sizeof(double)*dim);
  bin_coeff = (double **)      malloc(sizeof(double*)*dim); 
  ptr       = (struct item **) malloc(sizeof(struct item*)*dim);

  /* generate all multiindices i with 0 <= |i| <= d and 
     store them in multi */

  multi = (int ***) malloc((d+2)*sizeof(int**));
  d_real = d;
  for (m=1; m<=d; m++)
  { dim = binomi(p+m-1,m);
    multi[m] = (int **) malloc(dim*sizeof(int*)); 
    for(j=0; j<dim; j++)  /* address in multi[m][j][p] */
      multi[m][j] = (int *) malloc(sizeof(int)*(p+1));   
    multiind(p,m,d,multi[m]); 
  }  

  /* Calculate all coefficients c_{i,j} 
     and store them in coeff_list */

  for (m=0; m<dim; m++)
  { ptr[m]       = NULL;
    factor[m]    = 0;
    bin_coeff[m] = (double *) malloc(dimi*sizeof(double));
  }
  for (i=1; i<=d; i++) 
  { /* for-loop over the absolute value of the multiindex i */
    dimi = binomi(p+i-1,i);
    coeffpart1 = 1;     /* calculation of (|i|/d)^{|i|} */
    for (kk=0; kk<i; kk++)
      coeffpart1 = (coeffpart1 * i) / d_real; 
    for (ii=0; ii<dimi; ii++) 
    { /* loop over all multiindices with absolute value i */
      for (l=1; l<i; l++)
      { /* for-loop over the absolute value = l < i */
        diml = binomi(p+l-1,l);
        coeffpart2 = 1;      /* calculation of (|k|/d)^{|i|} */
        for (kk=0; kk<i; kk++) 
          coeffpart2 = (coeffpart2 * l) / d_real;
        for (kk=0; kk<i-l; kk++) /* calculation of (|k|/d)^{|i|}*(-1)^{|i-k|}*/
          coeffpart2 *= -1;
        for (k=0; k<diml; k++)
        { /* loop over all multiindices k with |k| = l */ 
          if (multi[i][ii][p] > multi[l][k][p] ) 
          { /* else: k > i -->  hence addend = 0 */
            addend = 1;
            for (kk=0; kk<p; kk++) /* calculation of i over k */
              addend *= binomi(multi[i][ii][kk],multi[l][k][kk]);
            if (addend != 0) 
            { /* i over k not equal 0, i.e. k < i */
              addend *= coeffpart2;
              for (m=0; m<dim; m++) 
                factor[m] =factor[m] + addend*bin_coeff[m][multi[l][k][p]-1];
            } 
          } 
        }  /* loop over all multiindices k with |k| = l */
      }  /* for-loop over l */ 
         /* nun k = i */
      for (m=0; m<dim; m++)
      { addend = 1;
        for (kk=0; kk<p; kk++)
          addend *= dbinomi((d_real*multi[i][ii][kk])/i,multi[d][m][kk]);
        bin_coeff[m][multi[i][ii][p]-1] = addend;
        factor[m] = factor[m] + addend * coeffpart1; 
        if (factor[m] != 0)  /* olvo 980804 [m] */
        { if (ptr[m] == NULL)
            ptr[m] = &coeff_list[m];
          else
          { ptr[m]->next = (struct item *) malloc(sizeof(struct item));
            ptr[m] = ptr[m]->next;
          }
          ptr[m]->a = multi[i][ii][p];
          ptr[m]->b = i;
          ptr[m]->c = factor[m];
          factor[m] = 0;
        }  
      }
    } 
  } /* end for-loop over all multiindices i */
  for (m=0; m<dim; m++)
     ptr[m]->next = NULL; 
  for (m=1; m<=d; m++)
  { dim = binomi(p+m-1,m);
    for(j=0; j<dim; j++)  /* address in multi[m][j][p] */
      free((char *) multi[m][j]);
    free((char *) multi[m]);
  } 
  free((char *) multi);
  dim  = binomi(p+d-1,d);
  for(m=0;m<dim;m++)
    free((char *) bin_coeff[m]);
  free((char *) bin_coeff);
  free((char *) factor);
  free((char *) ptr);
}


/****************************************************************************/
/*                                                           MORE UTILITIES */

/*--------------------------------------------------------------------------*/
void multma3vec2( int n, int p, int d, int bd, 
                  double ***X, double **S, int **jm )
{ int i,j,k;
  double sum;

  for (i=0; i<n; i++)
    for (k=0; k<bd; k++)
    { sum = 0;
      for (j=0; j<p; j++) 
        sum += S[i][j]*jm[k][j];
      X[i][k][0] = sum;
      for (j=1; j<d; j++)
        X[i][k][j] = 0;
    } 
}

/*--------------------------------------------------------------------------*/
void multma2vec2( int n, int p, int bd, double **X, double **S, int **jm )
{ int i,j,k;
  double sum;

  for (i=0; i<n; i++)
    for (k=0; k<bd; k++)
    { sum = 0;
      for (j=0; j<p; j++)
        sum += S[i][j]*jm[k][j];
      X[i][k] = sum;
    }
}

/*--------------------------------------------------------------------------*/
void multma2vec1( int n, int p, int d, double **X, double **S, int *jm )
{ int i,j,k;
  double sum;

  for (i=0; i<n; i++) 
  { sum = 0;
    for (j=0; j<p; j++)
      sum += S[i][j]*jm[j];
    X[i][1] = sum;
    for (j=2; j<d; j++)
      X[i][j] = 0;
  }
}


/****************************************************************************/
/*----------------------------------------------------------------------------
  From:
  Olaf Vogel 
  Institute of Scientific Computing
  TU Dresden
  diploma thesis "Zur Berechnung von Rand und Randfaltungen"
  File:    dvGausz.C
  LU-Factorisation and triangular backward/forward substitution
   
----------------------------------------------------------------------------*/

/* test if zero */
#define ZERO 1.0E-15

/*--------------------------------------------------------------------------*/
int LUFactorization( double** J, int n, int* RI, int* CI )
{ int i, j, k, cIdx, rIdx, h;
  double v;
    
  for (i=0; i<n; i++)
    RI[i]=i;
  for (j=0; j<n; j++)
    CI[j]=j;
  /* n Gausz-steps with full Pivoting */            
  for (k=0; k<n; k++)
  { v=0.0; cIdx=rIdx=0;
    /* Pivotsearch */
    for (i=k; i<n; i++)
      for (j=k; j<n; j++)
        if (fabs(J[RI[i]][CI[j]])>v) 
          { v=fabs(J[RI[i]][CI[j]]);
            rIdx=i; cIdx=j;
          }
    if (ZERO > v) 
    { fprintf(DIAG_OUT,
              "Error:LUFactorisation(..): no Pivot in step %d (%le)\n",k+1,v);
      return -(k+1);
    }
    /* row and column change resp. */
    if (rIdx > k)
    { h=RI[k]; 
      RI[k]=RI[rIdx]; 
      RI[rIdx]=h;
    }
    if (cIdx > k)
    { h=CI[k]; 
      CI[k]=CI[cIdx]; 
      CI[cIdx]=h;
    }
    /* Factorisation step */
    for (i=k+1; i<n; i++) 
    { /* L-part */
      J[RI[i]][CI[k]]/=J[RI[k]][CI[k]];
      /* R-part */
      for (j=k+1; j<n; j++)
        J[RI[i]][CI[j]]-=J[RI[i]][CI[k]]*J[RI[k]][CI[j]];
    }
  }
  return k; 
}

/*--------------------------------------------------------------------------*/
void GauszSolve( double** J, int n, int* RI, int* CI, double* b )
{ double* tmpZ;
  int i,j;

  tmpZ = myalloc1(n);
  for (i=0; i<n; i++)
  { tmpZ[i]=b[RI[i]];
    for (j=0; j<i; j++)
      tmpZ[i]-=J[RI[i]][CI[j]]*tmpZ[j];
  }
  for (i=n-1; i>=0; i--)
  { b[CI[i]]=tmpZ[i];
    for (j=i+1; j<n; j++)
      b[CI[i]]-=J[RI[i]][CI[j]]*b[CI[j]];
    b[CI[i]]/=J[RI[i]][CI[i]];
  }
  free(tmpZ);
}


/****************************************************************************/
int jac_solv( unsigned short tag, int n, double* x, double* b, 
              unsigned short sparse, unsigned short mode )
{ static double **J;
  static double **I;
  static double *y;
  static double *xold;
  static int* ri;
  static int* ci;
  static int nax,tagold,modeold,cgd;
  int i,j;
  int rc = 3;

  if ((n != nax) || (tag != tagold)) {
    if (nax)
    { free(*J); free(J);
      free(*I); free(I);
      free(xold); free(ri); free(ci); free(y);
    }
    J = myalloc2(n,n);
    I = myalloc2(n,n);
    y = myalloc1(n);
    
    xold = myalloc1(n);
    ri = (int*)malloc(n*sizeof(int));
    ci = (int*)malloc(n*sizeof(int));
    for (i=0; i<n; i++)
    { xold[i] = 0;
      for (j=0;j<n;j++)
        I[i][j]=(i==j)?1.0:0.0;
    }
    cgd = 1;
    modeold = 0;
    nax = n;
    tagold = tag;
  }
  if (cgd == 0)
    for (i=0; i<n; i++)
      if (x[i] != xold[i])
        cgd = 1;
  if (cgd == 1)
    for (i=0; i<n; i++)
      xold[i] = x[i];
  switch(mode) {
    case 0:
      mindec(rc,zos_forward(tag,n,n,1,x,y));
      mindec(rc,fov_reverse(tag,n,n,n,I,J));
      break;
    case 1:
      if ((modeold == 0) || (cgd == 1))
      { mindec(rc,zos_forward(tag,n,n,1,x,y));
        mindec(rc,fov_reverse(tag,n,n,n,I,J));
      }
      if (LUFactorization(J,n,ri,ci) < 0)
        return -3;
      modeold = 1;
      break;
    case 2:
      if ((modeold < 1) || (cgd == 1))
      { mindec(rc,zos_forward(tag,n,n,1,x,y));
        mindec(rc,fov_reverse(tag,n,n,n,I,J));
        if (LUFactorization(J,n,ri,ci) < 0)
          return -3;
      }
      GauszSolve(J,n,ri,ci,b);
      modeold = 2;
      break;
  }
  cgd = 0;
  return rc;
}


/****************************************************************************/
int inverse_Taylor_prop( unsigned short tag, int n, int d, 
                          double** Y, double** X )
{ int i,j,l,q;
  static double **I;
  register double bi;
  static double** Xhelp;
  static double** W;
  static double* xold;
  static double ***A;
  static double *w;
  static int *dd;
  static double *b;
  static int nax,dax,bd,cgd;
  static short **nonzero;
  short* nz;
  double* Aij;
  double* Xj;
  int ii, di, da, Di;
  int rc = 3;

  /* Re/Allocation Stuff */
  if ((n != nax) || (d != dax))
  { if (nax)
    { free(**A); free(*A); free(A);
      free(*I); free(I);
      free(*W); free(W);
      free(*Xhelp); free(Xhelp);
      free(w); 
      free(xold);
      free(*nonzero); free(nonzero);
      free(dd); free(b);
    }
    A = myalloc3(n,n,d+1);
    I = myalloc2(n,n);
    W = myalloc2(n,d);
    Xhelp = myalloc2(n,d);
    w = myalloc1(n);
    dd = (int*)malloc((d+1)*sizeof(int));
    b  = (double*)malloc(n*sizeof(double));
    xold = (double*)malloc(n*sizeof(double));
    nonzero = (short**)malloc(n*sizeof(short*)); 
    nz = (short*)malloc(n*n*sizeof(short)); 
    for (i=0; i<n; i++)
    { nonzero[i] = nz;
      nz = nz + n;
      xold[i] = 0;
      for (j=0; j<n; j++)
        I[i][j]=(i==j)?1.0:0.0;
    }
    cgd = 1;
    nax=n;
    dax=d;
    dd[0] = d+1;
    i = -1;
    while(dd[++i] > 1)
      dd[i+1] = (int)ceil(dd[i]*0.5);
    bd = i+1;
  }
  if (cgd == 0)
    for (i=0; i<n; i++)
      if (X[i][0] != xold[i])
        cgd = 1;
  if (cgd == 1)
  { cgd = 0;
    for (i=0; i<n; i++)
      xold[i] = X[i][0];
    mindec(rc,jac_solv(tag,n,xold,b,0,1));
    if (rc == -3)
      return -3;
  }
  ii = bd;
  for (i=0; i<n; i++)
    for (j=0; j<d; j++)  
      Xhelp[i][j] = X[i][j+1];

  while (--ii > 0)
  { di = dd[ii-1]-1;
    Di = dd[ii-1]-dd[ii]-1;
    mindec(rc,hos_forward(tag,n,n,di,Di+1,xold,Xhelp,w,W));
    mindec(rc,hov_reverse(tag,n,n,Di,n,I,A,nonzero));
    da = dd[ii];
    for (l=da; l<dd[ii-1]; l++)
    { for (i=0; i<n; i++)  
      { if (l == 0)
          bi = w[i]-Y[i][0];
        else
          bi = W[i][l-1]-Y[i][l];
        for (j=0; j<n; j++)
          if (nonzero[i][j]>1)
          { Aij = A[i][j];
            Xj = X[j]+l;
            for (q=da; q<l; q++)
              bi += (*(++Aij))*(*(--Xj));
	  }
	b[i] = -bi;
      }
      mindec(rc,jac_solv(tag,n,xold,b,0,2));
      if (rc == -3)
        return -3;
      for (i=0; i<n; i++)
      { X[i][l] += b[i];
        /* 981214 new nl */
        Xhelp[i][l-1] += b[i];
      }
    }
  } 
  return rc;
}
 

/****************************************************************************/
int inverse_tensor_eval( int tag, int n, int d, int p, 
                          double *x, double **tensor, double** S )
{ static int dim;
  static int dold,pold;
  static struct item *coeff_list;
  int i,j,k,dimten;
  int *it = (int*) malloc(d*sizeof(int));
  double** X;
  double** Y;
  int *jm;
  double *y = (double*) malloc(n*sizeof(double));
  struct item *ptr;
  int rc = 3;

  dimten=binomi(p+d,d);
  for(i=0;i<n;i++)
    for(j=0;j<dimten;j++)
      tensor[i][j] = 0;
  mindec(rc,zos_forward(1,n,n,0,x,y));
  if (d > 0)
  { if ((d != dold) || (p != pold))
    { if (pold)
      { /* olvo 980728 */
        dim = binomi(pold+dold-1,dold);
        freecoefflist(dim,coeff_list);
        free((char*) coeff_list);
      }
      dim = binomi(p+d-1,d);
      coeff_list = (struct item *) malloc(sizeof(struct item)*dim);
      coeff(p,d, coeff_list);
      dold = d;
      pold = p;
    }
    jm = (int *)malloc(sizeof(int)*p);
    X = myalloc2(n,d+1);
    Y = myalloc2(n,d+1);
    for (i=0; i<n; i++)
    { X[i][0] = x[i];
      for (j=1; j<d; j++)
        X[i][j] = 0;
      Y[i][0] = y[i];
    }
    if (d == 1)
    { it[0] = 0;
      for (i=0; i<dim; i++)  /* sum over all multiindices jm with |jm| = d */
      { it[0] = it[0]+1;
        convert(p,d,it,jm);
        ptr = &coeff_list[i];
        multma2vec1(n,p,d,Y,S,jm);
        mindec(rc,inverse_Taylor_prop(tag,n,d,Y,X));
        if (rc == -3)
          return -3;              
        do 
        { for(j=0;j<n;j++)
            tensor[j][ptr->a] += X[j][ptr->b]*ptr->c;
          ptr = ptr->next;
        } while (ptr != NULL);
      }
    }
    else
    { for (i=0; i<d-1; i++)
        it[i] = 1;
      it[d-1] = 0;
      for (i=0; i<dim; i++)  /* sum over all multiindices jm with |jm| = d */ 
      { it[d-1] = it[d-1]+1;
        for (j=d-2; j>=0; j--)
          it[j] = it[j] + it[j+1]/(p+1);
        for (j=1; j<d; j++) 
          if (it[j] > p) it[j] = it[j-1];
        convert(p,d,it,jm);
        multma2vec1(n,p,d,Y,S,jm); /* Store S*jm in Y */
        mindec(rc,inverse_Taylor_prop(tag,n,d,Y,X));
        if (rc == -3)
          return -3;              
        ptr = &coeff_list[i];
        do 
        { for(j=0;j<n;j++)
            tensor[j][ptr->a] += X[j][ptr->b]*ptr->c;
          ptr = ptr->next;
        } while (ptr != NULL);
      }
    }
    free((char*) jm); 
    free((char*) *X); free((char*) X);
    free((char*) *Y); free((char*) Y);
   }
  for(i=0;i<n;i++)
    tensor[i][0] = x[i];
  free((char*) y); 
  free((char*) it); 
  return rc;
}


/****************************************************************************/
int tensor_eval( int tag, int m, int n, int d, int p, 
                  double* x, double **tensor, double **S )
{ static int bd,dim;
  static int dold,pold;
  static struct item *coeff_list;
  int i,j,k,dimten,ctr;
  int **jm, jmbd;
  int *it = (int*) malloc(d*sizeof(int));
  double *y = (double*) malloc(m*sizeof(double));
  double*** X;
  double*** Y;
  struct item *ptr[10];
  int rc = 3;

  dimten=binomi(p+d,d);
  for (i=0; i<m; i++)
    for (j=0; j<dimten; j++)
      tensor[i][j] = 0;
  if (d == 0) 
  { mindec(rc,zos_forward(1,m,n,0,x,y));
  }
  else
  { if ((d != dold) || (p != pold)) 
    { if (pold)
      { dim = binomi(pold+dold-1,dold);
        freecoefflist(dim,coeff_list);
        free((char*) coeff_list);
      }
      dim = binomi(p+d-1,d);
      if (dim < 10)
        bd = dim;
      else
        bd = 10;
      coeff_list = (struct item *) malloc(sizeof(struct item)*dim);
      coeff(p,d, coeff_list);
      dold = d;
      pold = p;
    }
    jmbd = bd;
    jm = (int **) malloc(jmbd*sizeof(int*));
    for (i=0; i<jmbd; i++)
      jm[i] = (int *) malloc(p*sizeof(int));
    if (d == 1) 
    { X = myalloc3(1,n,bd);
      Y = myalloc3(1,m,bd);
      ctr   = 0;
      it[0] = 0;
      for (i=0; i<dim; i++) /* sum over all multiindices jm with |jm| = d */
      { it[0] = it[0]+1; 
        convert(p,d,it,jm[ctr]);
        ptr[ctr] = &coeff_list[i];
        if (ctr < bd-1)
          ctr += 1;
        else
        { multma2vec2(n,p,bd,X[0],S,jm);
          mindec(rc,fov_forward(tag,m,n,bd,x,X[0],y,Y[0]));
          for (k=0; k<bd; k++)
            do 
            { for (j=0; j<m; j++)
                tensor[j][ptr[k]->a] += Y[0][j][k]*ptr[k]->c;
              ptr[k] = ptr[k]->next;
            } while (ptr[k] != NULL);
          if (dim-i < bd)
            bd = dim-i-1;
          ctr = 0;
        }
      }
    }
    else
    { X = myalloc3(n,bd,d);
      Y = myalloc3(m,bd,d);
      ctr = 0; 
      for (i=0; i<d-1; i++)
        it[i] = 1;
      it[d-1] = 0;
      for (i=0; i<dim; i++) /* sum over all multiindices jm with |jm| = d */   
      { it[d-1] = it[d-1]+1;
        for (j=d-2; j>=0; j--)
          it[j] = it[j] + it[j+1]/(p+1);
        for (j=1; j<d; j++) 
          if (it[j] > p) 
            it[j] = it[j-1];
        convert(p,d,it,jm[ctr]);
        ptr[ctr] = &coeff_list[i];
        if (ctr < bd-1)
          ctr += 1;
        else
        { multma3vec2(n,p,d,bd,X,S,jm);
          mindec(rc,hov_forward(tag,m,n,d,bd,x,X,y,Y));
          for (k=0; k<bd; k++)
            do 
            { for (j=0; j<m; j++)
                tensor[j][ptr[k]->a] += Y[j][k][ptr[k]->b-1]*ptr[k]->c;
              ptr[k] = ptr[k]->next;
            } while (ptr[k] != NULL);       
          if (dim-i < bd)
            bd = dim-i-1;
          ctr = 0;
        }
      }
    }
    for (i=0; i<jmbd; i++)
      free((char*) *(jm+i));
    free((char*) jm);
    free((char*) **X); free((char*) *X); free((char*) X);
    free((char*) **Y); free((char*) *Y); free((char*) Y); 
  }
  for(i=0;i<m;i++)
    tensor[i][0] = y[i];
  bd = jmbd;
  free((char*) y); 
  free((char*) it); 
  return rc;
}


/****************************************************************************/
void tensor_value( int d, int m, double *y, double **tensor, int *multi )
{ int i, j, max, ind, add;
  int *im = (int*) malloc(d*sizeof(int));

  max = 0;
  ind = d-1;
  for (i=0; i<d; i++)
  { if (multi[i] > max)
      max = multi[i];
    im[i] = 0; 
  }
  for (i=0; i<d; i++)
  { if (multi[i] == max)  /* olvo 980728 == instead of = */
    { im[ind] = multi[i];
      multi[i] = 0;
      max = 0;
      ind -= 1;
      for (j=0; j<d; j++)
        if (multi[j] > max)
          max = multi[j];
    }
  }
  add = address(d,im);
  for (i=0; i<m; i++)
    y[i] = tensor[i][add];
  free((char*) im); 
}


/****************************************************************************/
/*                                                               THAT'S ALL */
#ifdef __cplusplus
}
#endif

#undef _ADOLC_SRC_
#undef _TAYLORC_C_


























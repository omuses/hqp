/*
   --------------------------------------------------------------
   File gaussexam.C of ADOL-C version 1.8.0       as of Dec/01/98
   --------------------------------------------------------------
   based on gaussexam.C version ADOL-C 1.7

   Example: Gaussian elimination

   Last changes:
     981201  olvo:  new headers
     980821  olvo:  some little changes in output design

   --------------------------------------------------------------
*/

/****************************************************************************/
/*                                                                 INCLUDES */
#include "../../../adolc/adolc.h"
#include <math.h>

/****************************************************************************/
/*                                              active Gaussian elimination */
void gausselim( int n, adoublem& A, adoublev& bv )
{ along i;
  adoublev temp(n);
  adouble r, rj, temps;
  int j, k, ik;
  for (k=0; k<n; k++)
  { for (j=0; j<n; j++)
      fprintf(stdout,"%14.6le ",value(A[k][j]));
    fprintf(stdout,"      %14.6le\n",value(bv[k]));
  }
  fprintf(stdout,"initial state ----------------------\n");
 
/*--------------------------------------------------------------------------*/
  for (k=0; k<n; k++)                                   /* elimination loop */
  { i = k;

/*--------------------------------------------------------------------------*/
    r = fabs(A[k][k]);                                          /* Pivoting */
    for (j=k+1; j<n; j++)
    {  rj = fabs(A[j][k]); /* look for a greater element in the same column */
       condassign(i,(rj >r),j);
       condassign(r,(rj >r),rj);
    }
    fprintf(stdout,"%lf index\n",value(i));

/*--------------------------------------------------------------------------*/
    temp = A[i]; A[i] = A[k]; A[k] = temp;             /* exchange  of rows */
    temps = bv[i]; bv[i] = bv[k]; bv[k] = temps;
    if (!value(A[k][k]))
    { fprintf(stdout," Matrix does not have full rank!\n");
      exit(-1);
    }

    fprintf(stdout,"changed rows: -----------------------\n");
    for (ik=0; ik<n; ik++)
    { for (j=0; j<n; j++)
        fprintf(stdout,"%14.6le ",value(A[ik][j]));
      fprintf(stdout,"      %14.6le\n",value(bv[ik]));
    }

/*--------------------------------------------------------------------------*/
    temps = A[k][k];                                    /* elimination step */
    A[k] = A[k]/temps;
    bv[k] = bv[k]/temps;
    for (j=k+1; j<n; j++)
    { temps = A[j][k];
      A[j] -= temps*A[k];
      bv[j] -= temps*bv[k];
    }

    fprintf(stdout,"step:--------------------------------\n");
    for (ik=0; ik<n; ik++)
    { for(j=0; j<n; j++)
        fprintf(stdout,"%14.6le ",value(A[ik][j]));
      fprintf(stdout,"      %14.6le\n",value(bv[ik]));
    }
  } // endfor elimination loop

  temp = 0.0;
  for(k=n-1; k>=0; k--)
  { temp[k] = (bv[k]-(A[k]*temp))/A[k][k];
    fprintf(stdout,"%14.6le\n",value(temp[k]));
  }
  bv = temp;
  return;
}


/****************************************************************************/
/*                                                             MAIN PROGRAM */
int main() 
{
  int i, j, k, ok = 1;
/*--------------------------------------------------------------------------*/
  short tag = 1;                        /* variables and problem parameters */
  double epsilon = 0.0000000001; // max. allowed difference between results
  int dum            = 1;
  const int max_deg  = 4;        // maximal order of derivation
  const int tayl_num = 2;        // Number of taylor series
  const int size     = 5;
  const int indep    = size*size+size;
  const int depen    = size;
  const int laglength= 2;
  int N              = size*size;

/*--------------------------------------------------------------------------*/
  double* lagras = myalloc(depen);                           /* Lagrangians */
  for (i=0; i<depen; i++)
    lagras[i] = i+1;
  double** lagrav = myalloc(laglength,depen);
  for (j=0; j<laglength; j++)
    for (i=0; i<depen; i++)
      lagrav[j][i] = j+i+1;

/*--------------------------------------------------------------------------*/
  short** nonzero = new short*[laglength];               /* nonzero pattern */
  for (i=0; i<laglength; i++)
    nonzero[i] = new short[indep];

/*--------------------------------------------------------------------------*/
  double**  resultshos = myalloc(indep,max_deg);         /* lots of tensors */
  double*** resultshov = myalloc(laglength,indep,max_deg);
  double*   resultsfos = myalloc(indep);
  double**  resultsfov = myalloc(laglength,indep);
  double*   valuepoint = myalloc(depen);
  double*** arguments  = myalloc(indep,tayl_num,max_deg);
  double**  scalvaluep = myalloc(tayl_num,depen);
  double*** scalargs   = myalloc(tayl_num,indep,max_deg);
  double*** scalres    = myalloc(tayl_num,depen,max_deg);
  double*** taylors    = myalloc(depen,tayl_num,max_deg);

/*--------------------------------------------------------------------------*/
  double yp[size], xp[size*size+size];                  /* passive variable */

/*--------------------------------------------------------------------------*/
  adoublem A(size,size);                                /* active variables */
  adoublev bv(size);         

/*--------------------------------------------------------------------------*/
  trace_on(tag,dum);             /* Taping all calculations with 'adoubles' */
  for (i=0; i<size; i++)
    for (j=0; j<size; j++)
      A[i][j] <<= pow((double)1+j,(double)i); /* indep. vars */
  for (i=0; i<size; i++)
    bv[i] <<= -i-1; /* indep. vars */
  gausselim(size,A,bv);
  bv >>= yp; /* dep. vars */
  trace_off(); 


/*--------------------------------------------------------------------------*/
  /* xp alias A+bv */           /* Initializations for forward calculations */
  for (i=0; i<size; i++) 
    for (j=0; j<size; j++)
      xp[i*size+j] = pow((double)1+j,(double)i); 
  for (i=0; i<size; i++)
    xp[N+i] = -i-1;

  /* directional taylors (input) */
  for (i=0; i<indep; i++)
    for (j=0; j<tayl_num; j++)
      for (k=0; k<max_deg; k++)
      { arguments[i][j][k]=i+j+k+3;
        scalargs[j][i][k]=i+j+k+3;
      } 

/*--------------------------------------------------------------------------*/
  for (j=0; j<tayl_num; j++)                        /* forward calculations */
    hos_forward(tag,depen,indep,max_deg,1,
                xp,scalargs[j],scalvaluep[j],scalres[j]);
  hov_forward(tag,depen,indep,max_deg,tayl_num,
              xp,arguments,valuepoint,taylors);

  // test for correctness +++++++++++++++++++++++++++++++++++++++++++++++
  for (i=0; i<depen; i++)
  { fprintf(stdout,"dependent variable number %d\n",i);
    for (j=0; j<tayl_num; j++)
    { fprintf(stdout,"taylor serie number %d\n",j);
      fprintf(stdout,"hov_f. valuepoint[%d]: %14.6le",i,valuepoint[i]);
      fprintf(stdout,"  =?  hos_f. scalvaluep[%d][%d]: %14.6le",
                     j,i,scalvaluep[j][i]);
      fprintf(stdout,"  =?  yp[%d]: %14.6le\n",i,yp[i]);
      if (fabs(valuepoint[i]-scalvaluep[j][i]) > epsilon)
      { fprintf(stdout,
                "difference is here <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n");
        ok = 0;
      }
      for (k=0; k<max_deg; k++)
      { fprintf(stdout,"hov_f. taylors[%d][%d][%d]: %14.6le",
                       i,j,k,taylors[i][j][k]);
        fprintf(stdout,"  =?  hos_f. scalres[%d][%d][%d]: %14.6le",
                       j,i,k,scalres[j][i][k]);
        fprintf(stdout,"  <- inp. coeff.: %3d\n",((int) scalargs[j][i][k]));
        if (fabs(taylors[i][j][k]-scalres[j][i][k]) > epsilon)
        { fprintf(stdout,
                  "difference is here <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n");
          ok = 0;
        } 
      } 
    } 
  }

/*--------------------------------------------------------------------------*/
                                       /* higher order reverse calculations */
  hos_forward(tag,depen,indep,max_deg,max_deg,
              xp,scalargs[0],scalvaluep[0],scalres[0]);
  fprintf(stdout,
          "reverse sweeps will be done for the first taylor serie only\n");

  hos_reverse(tag,depen,indep,max_deg-1,lagras,resultshos);
  /* fprintf(stdout,"resultshos after hos_reverse :\n");
  for (myi=0; myi<indep; myi++)
  { for (myj=0; myj<max_deg; myj++)
      fprintf(stdout," %14.6le",resultshos[myi][myj]);
    fprintf(stdout,"\n");
  }
  fprintf(stdout,"\n"); */

  hov_reverse(tag,depen,indep,max_deg-1,laglength,lagrav,resultshov,nonzero);
  /* fprintf(stdout,"resultshov after hov_reverse :\n");
  for (myp=0; myp<laglength; myp++)
    for (myi=0; myi<indep; myi++)
    { for (myj=0; myj<max_deg; myj++)
        fprintf(stdout," %14.6le",resultshov[myp][myi][myj]);
      fprintf(stdout,"\n");
    }
  fprintf(stdout,"\n"); */

/*--------------------------------------------------------------------------*/
                                        /* first order reverse calculations */
  hos_forward(tag,depen,indep,max_deg,1,
              xp,scalargs[0],scalvaluep[0],scalres[0]);

  fos_reverse(tag,depen,indep,lagras,resultsfos);
  fov_reverse(tag,depen,indep,laglength,lagrav,resultsfov);

/*--------------------------------------------------------------------------*/
                                                       /* output of results */
  for (i=0; i<laglength; i++)
    if (i == 0)
      for (j=0; j<indep; j++)    
        for (k=0; k<max_deg; k++)
          if (k == 0)
          { fprintf(stdout,"reshov[%d][%d][%d]: %14.6le",
                           i,j,k,resultshov[i][j][k]);
            fprintf(stdout," =? reshos[%d][%d]: %14.6le",
                           j,k,resultshos[j][k]);
            fprintf(stdout," =? resfov[%d][%d]: %14.6le",
                           i,j,resultsfov[i][j]);
            fprintf(stdout," =? resfos[%d]: %14.6le\n",
                           j,resultsfos[j]);
            if ((fabs(resultshov[i][j][k]-resultshos[j][k])  > epsilon) 
             || (fabs(resultshov[i][j][k]- resultsfov[i][j]) > epsilon) 
             || (fabs(resultshov[i][j][k]-resultsfos[j])     > epsilon))
            { fprintf(stdout,
                      "difference is here <<<<<<<<<<<<<<<<<<<"
                      "<<<<<<<<<<<<<<<<\n");
              ok = 0;
            }
          }
          else
          { fprintf(stdout,"reshov[%d][%d[%d]: %14.6le",
                           i,j,k,resultshov[i][j][k]);
            fprintf(stdout," =? reshos[%d][%d]: %14.6le\n",
                           j,k,resultshos[j][k]);
            if (fabs(resultshov[i][j][k]-resultshos[j][k]) > epsilon)
            { fprintf(stdout,
                      "difference is here <<<<<<<<<<<<<<<<<<<"
                      "<<<<<<<<<<<<<<<<\n");
              ok = 0;
            }
          } 
    else
      for (j=0; j<indep; j++)    
        for (k=0; k<max_deg; k++)
          if (k == 0)
          { fprintf(stdout,"reshov[%d][%d[%d]: %14.6le",
                           i,j,k,resultshov[i][j][k]);
            fprintf(stdout," =? resfov[%d][%d]: %14.6le\n",
                           i,j,resultsfov[i][j]);
            if (fabs(resultshov[i][j][k]- resultsfov[i][j]) > epsilon)
            {  fprintf(stdout,
                      "difference is here <<<<<<<<<<<<<<<<<<<"
                      "<<<<<<<<<<<<<<<<\n");
              ok=0;
            } 
          }
          else
            fprintf(stdout,"reshov[%d][%d[%d]: %14.6le\n",
                           i,j,k,resultshov[i][j][k]);

  for (i=0; i<laglength; i++)
    for (j=0; j<indep; j++)
      fprintf(stdout,"nonzero[%d][%d]: %3d\n",i,j,nonzero[i][j]);
  if (!ok)
    fprintf(stdout,"calculation  is not ok<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n"
                   " This message may be caused by very small differences "
                   "(not necessary \nrecognizable in the output)\n");

  return 1;
}   


/*
   --------------------------------------------------------------
   File LUdet.C of ADOL-C version 1.8.5           as of Nov/22/99
   --------------------------------------------------------------

   Example:  * Computation of the determinant of a matrix
               by LU-decomposition of the system matrix without pivoting 
             * application of tapedoc to observe taping of
               the new op_codes for the elementary operations
                  
                     y += x1 * x2;
                     y -= x1 * x2;           

   Last changes:
     990922 olvo    first version 

   --------------------------------------------------------------
*/

/****************************************************************************/
/*                                                                 INCLUDES */
#include "LU.h"              // LU-decomposition  


/****************************************************************************/
/*                                                             MAIN PROGRAM */
int main() 
{ /*------------------------------------------------------------------------*/
  /* variables */
  const int tag   = 1;                       // tape tag
  const int size  = 5;                       // system size
  const int indep = size*size;               // # of indeps
  const int depen = 1;                       // # of deps

  double  A[size][size], a1[size], a2[size], det; // passive variables 
  adouble **AA, *AAp,     Adet;                   // active variables
  double *args  = myalloc1(indep);                // arguments
  double *grad  = myalloc1(indep);                // the gradient
  double **hess = myalloc2(indep,indep);          // the hessian

  int i,j;                  


  /*------------------------------------------------------------------------*/
  /* Info */
  fprintf(stdout,"DETERMINANT by LU-DECOMPOSITION (ADOL-C Example)\n\n");


  /*------------------------------------------------------------------------*/
  /* Allcoation und initialization of the system matrix */
  AA  = new (adouble*)[size];
  AAp = new adouble[size*size];
  for (i=0; i<size; i++)
  { AA[i] = AAp; 
    AAp += size;
  }
  for(i=0; i<size; i++)
  { a1[i] = i*0.25;
    a2[i] = i*0.33;
  }
  for(i=0; i<size; i++) 
  { for(j=0; j<size; j++) 
      A[i][j] = a1[i]*a2[j];
    A[i][i] += i+1;
  }


  /*------------------------------------------------------------------------*/
  /* Taping the computation of the determinant */
  trace_on(tag);  
    /* marking indeps */           
    for(i=0; i<size; i++)
      for(j=0; j<size; j++)
        AA[i][j] <<= (args[i*size+j] = A[i][j]);     // indep. vars 
    /* LU-factorization and computation of determinant */
    LUfact(size,AA);
    Adet = AA[0][0];
    for (i=1; i<size; i++)
      Adet *= AA[i][i];
    /* marking deps */
    Adet >>= det;
  trace_off(); 
  fprintf(stdout," Determinant (original):  %16.4le\n",det);       


  /*------------------------------------------------------------------------*/
  /* Recomputation of determinant */
  function(tag,depen,indep,args,&det);
  fprintf(stdout," Determinant (from tape): %16.4le\n",det);    


  /*------------------------------------------------------------------------*/
  /* Computation of gradient */
  gradient(tag,indep,args,grad);
  fprintf(stdout," Gradient:\n");
  for (i=0; i<size; i++)
  { for (j=0; j<size; j++)
      fprintf(stdout," %14.6le",grad[i*size+j]);
    fprintf(stdout,"\n");
  }

  /*------------------------------------------------------------------------*/
  /* Computation of hessian */
  hessian(tag,indep,args,hess); 
  fprintf(stdout," Part of Hessian:\n");
  for (i=0; i<size; i++)
  { for (j=0; j<size; j++)
      fprintf(stdout," %14.6le",hess[0][i*size+j]);
    fprintf(stdout,"\n");
  }

  /*------------------------------------------------------------------------*/
  /* Tape-documentation */
  tape_doc(tag,depen,indep,args,&det);

  /*------------------------------------------------------------------------*/
  /* Tape statistics */
  int tape_stats[11];
  tapestats(tag,tape_stats);
  fprintf(stdout,"\n    independents   %d\n",tape_stats[0]);
  fprintf(stdout,"    dependents     %d\n",tape_stats[1]);
  fprintf(stdout,"    operations     %d\n",tape_stats[5]);
  fprintf(stdout,"    buffer size    %d\n",tape_stats[4]);
  fprintf(stdout,"    maxlive        %d\n",tape_stats[2]);
  fprintf(stdout,"    valstack size  %d\n\n",tape_stats[3]);
 
  /*------------------------------------------------------------------------*/
  /* That's it */
  return 1;
}









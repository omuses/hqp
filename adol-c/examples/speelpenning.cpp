/*
   --------------------------------------------------------------
   File speelpenning.C of ADOL-C version 1.8.0    as of Dec/01/98
   --------------------------------------------------------------

   Example: Speelpenning's product, described in the manual

   Last changes:
     981201 olvo   last check (new headers)
     980805 andrea first version

   --------------------------------------------------------------
*/

/****************************************************************************/
/*                                                                 INCLUDES */
#include "../adolc/adouble.h"            // use of active doubles and taping
#include "../adolc/drivers/drivers.h"    // use of "Easy to Use" drivers 
                                         // gradient(.) and hessian(.)

#include <iostream>
using namespace std;

#include <malloc.h>
#include <math.h>

/****************************************************************************/
/*                                                             MAIN PROGRAM */
int main() 
{ int n,i,j;
  int tape_stats[11];

  cout << "SPEELPENNINGS PRODUCT (ADOL-C Documented Example)\n\n"; 
  cout << "number of independent variables = ?  \n";
  cin >> n;

  double *xp = new double[n];          
  double  yp = 0.0;
  adouble *x = new adouble[n];         // or: adoublev x(n);
  adouble  y = 1;

  for(i=0; i<n; i++)
    xp[i] = (i+1.0)/(2.0+i);           // some initialization

  trace_on(1);                         // tag = 1, keep = 0 by default
    for(i=0; i<n; i++)
    { x[i] <<= xp[i];                  // or  x <<= xp outside the loop
      y *= x[i]; 
    } // end for
    y >>= yp;
    delete[] x;                        // not needed if x adoublev
  trace_off();

  tapestats(1,tape_stats);             // reading of tape statistics
  cout<<"maxlive "<<tape_stats[2]<<"\n";
  // ..... print other tape stats

  double* g = new double[n];           
  gradient(1,n,xp,g);                  // gradient evaluation

  double** H   = (double**) malloc(n*sizeof(double*));
  for(i=0; i<n; i++)
    H[i] = (double*)malloc((i+1)*sizeof(double)); 
  hessian(1,n,xp,H);                   // H equals (n-1)g since g is
  double errg = 0;                     // homogeneous of degree n-1.
  double errh = 0;
  for(i=0; i<n; i++) 
    errg += fabs(g[i]-yp/xp[i]);       // vanishes analytically.
  for(i=0; i<n; i++) 
  { for(j=0; j<n; j++) 
    { if (i>j)                         // lower half of hessian
        errh += fabs(H[i][j]-g[i]/xp[j]); 
    } // end for 
  } // end for
  cout << yp-1/(1.0+n) << " error in function \n";
  cout << errg <<" error in gradient \n";
  cout << errh <<" consistency check \n";

  return 1;
} // end main


/*
   --------------------------------------------------------------
   File odexam.C of ADOL-C version 1.8.0          as of Dec/01/98
   --------------------------------------------------------------

   Example: Nonlinear ordinary differential equation based on the
            Robertson test problem, described in the manual

   Last changes:
     981201 olvo   last check (new headers)
     981116 olvo   changed tracerhs(..)
     980805 andrea first version

   --------------------------------------------------------------
*/

/****************************************************************************/
/*                                                                 INCLUDES */
#include "../adolc/adouble.h"            // use of active double and taping
#include "../adolc/drivers/odedrivers.h" // use of "Easy To Use" ODE drivers 
#include "../adolc/adalloc.h"            // use of ADOL-C allocation utilities 

// NOTICE: If one wants to include all ADOL-C interfaces without
//         getting in trouble to find out the right header files
//         
//         #include "adolc.h"
//
//         will do the right work. 

#include <iostream>
using namespace std;

/****************************************************************************/
/*                                                          ADOUBLE ROUTINE */
void tracerhs(short int tag, double* py, double* pyprime) 
{ adoublev y(3);                // we left the parameters passive
  adoublev yprime(3);           // and use the vector types

  trace_on(tag);
    y <<= py;                   // initialize and mark independents
    yprime[0] = -sin(y[2]) + 1.0e8*y[2]*(1.0-1.0/y[0]);
    yprime[1] = -10.0*y[0] + 3.0e7*y[2]*(1-y[1]);
    yprime[2] = -yprime[0] - yprime[1];
    yprime >>= pyprime;        // mark and pass dependents
  trace_off(1);                // write tape array onto a file
} // end tracerhs

/****************************************************************************/
/*                                                             MAIN PROGRAM */
int main() 
{ int i,j,deg;
  int n = 3;
  double py[3];
  double pyp[3];

  cout << "MODIFIED ROBERTSON TEST PROBLEM (ADOL-C Documented Example)\n\n";
  cout << "degree of Taylor series =?\n";
  cin >> deg;

  short** nz = new short*[n];
  double  **X;
  double ***Z;
  double ***B;

  X = myalloc2(n,deg+1);
  Z = myalloc3(n,n,deg);
  B = myalloc3(n,n,deg);

  for(i=0; i<n; i++) 
  { py[i]   = (i==0) ? 1.0 : 0.0;       // Initialize the base point
    X[i][0] = py[i];                    // and the Taylor coefficient;
    nz[i]   = new short[n];             // set up sparsity array
  } // end for                        

  tracerhs(1,py,pyp);                   // trace RHS with tag = 1

  forode(1,n,deg,X);                    // compute deg coefficients
  reverse(1,n,n,deg-1,Z,nz);            // U defaults to the identity
  accode(n,deg-1,Z,B,nz);

  cout << "nonzero pattern:\n";
  for(i=0; i<n; i++) 
  { for(j=0; j<n; j++)
      cout << nz[i][j]<<"\t";
    cout <<"\n"; 
  } // end for 
  cout << "\n";
  cout << " 4 = transcend , 3 = rational , 2 = polynomial ,"
       << " 1 = linear , 0 = zero \n";
  cout << " negative number k indicate that entries of all"
       << " B_j with j < -k vanish \n";

  return 1;
} // end main

#define _SFUNC_EXPERIMENTAL_C_
/*
   --------------------------------------------------------------
   File sfunc_experimental.C 
   of ADOL-C version 1.8.0                        as of Dec/01/98
   --------------------------------------------------------------

   Example: function module containing an experimental function

   Each << function module >> contains:
          
     (1) const char* const controlFileName 
     (2) int indepDim; 
 
     (3) void initProblemParameters( void )
     (4) void initIndependents( double* indeps )
     (5) double originalScalarFunction( double* indeps )
     (6) double tapingScalarFunction( int tag, double* indeps )   

   Last changes: 
     981201 olvo new headers
     980805 olvo first version

   --------------------------------------------------------------
*/


/****************************************************************************/
/*                                                                 INCLUDES */
#include "../../../adolc/adolc.h"

#include <time.h>
#include <math.h>


/****************************************************************************/
/*                                                         GLOBAL VARIABLES */

/*--------------------------------------------------------------------------*/
/*                                                        Control file name */
const char* controlFileName = "experimental.ctrl";

/*--------------------------------------------------------------------------*/
/*                                                               Dimensions */
int indepDim;

/*--------------------------------------------------------------------------*/
/*                                       Other problem dependent parameters */
//#define CODE sqrt(indeps[i])
#define CODE sin(indeps[i])
//#define CODE indeps[i]*indeps[i] 


/****************************************************************************/
/*                                                  INIT PROBLEM PARAMETERS */
void initProblemParameters( void )
{ fprintf(stdout,"EXPERIMENTAL EXAMPLE (ADOL-C Example)\n\n");
  if (indepDim <= 0)
  { fprintf(stdout,"    number of independent variables = ? ");
    fscanf(stdin,"%d",&indepDim);
    fprintf(stdout,"\n");
  }
}


/****************************************************************************/
/*                                                        INITIALIZE INDEPs */
void initIndependents( double* indeps )
{ int i;
  for (i=0; i<indepDim; i++)
    indeps[i] = (i+1.0)/(2.0+i);
}


/****************************************************************************/
/*                                                 ORIGINAL SCALAR FUNCTION */

/*--------------------------------------------------------------------------*/
/*                                                    experimental function */
double experimental( int dim, double* indeps ) 
{ int i;
  double y = 1.0;
  for (i=0; i<dim; i++)
    y += CODE; 
  return y;
}

/*--------------------------------------------------------------------------*/
/*                                                   The interface function */
double originalScalarFunction( double* indeps )
{ return experimental(indepDim, indeps);
}


/****************************************************************************/
/*                                                   TAPING SCALAR FUNCTION */

/*--------------------------------------------------------------------------*/
/*                                             active experimental function */
adouble activeExperimental( int dim, adouble* indeps ) 
{ int i;
  adouble y = 1.0;
  for (i=0; i<dim; i++)
    y += CODE; 
  return y;
}

/*--------------------------------------------------------------------------*/
/*                                                   The interface function */
double tapingScalarFunction( int tag, double* indeps )
{ int i;
  trace_on(tag);
  adouble* activeIndeps = new adouble[indepDim];
  adouble* aIP = activeIndeps;
  double*  iP  = indeps; 
  for (i=0; i<indepDim; i++)
     *aIP++ <<= *iP++; 
  adouble ares = activeExperimental(indepDim, activeIndeps);
  double res = 0;
  ares >>= res;
  trace_off();
  return res; 
}

#undef _SFUNC_EXPERIMENTAL_C_








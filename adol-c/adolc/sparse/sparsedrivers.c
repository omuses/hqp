/*---------------------------------------------------------------------------- 
 ADOL-C -- Automatic Differentiation by Overloading in C++
 File:     sparse/sparsedrivers.c
 Revision: $Id: sparsedrivers.c,v 1.1 2004/10/13 14:18:12 e_arnold Exp $
 Contents: All "Easy To Use" C++ interfaces of SPARSE package

 Copyright (c) 2004
               Technical University Dresden
               Department of Mathematics
               Institute of Scientific Computing
  
 This file is part of ADOL-C. This software is provided under the terms of
 the Common Public License. Any use, reproduction, or distribution of the
 software constitutes recipient's acceptance of the terms of this license.
 See the accompanying copy of the Common Public License for more details.

 History: 20040414 kowarz:  adaption to configure - make - make install
          19990308 christo: myalloc1_ushort -> myalloc1_uint
          19990302 christo: new interface of jac_pat(...)
          19981130 olvo:    newly created from driversc.c

----------------------------------------------------------------------------*/
 
#include "../sparse/sparse.h"
#include "../sparse/jacutils.h"

#include <malloc.h>
#include <math.h>

BEGIN_C_DECLS

/*---------------------------------------------------------------------------*/
/*                                                         jacobian pattern  */
/*                                                                           */

int jac_pat( 
      short        tag,        /* tape identification                        */
      int          depen,      /* number of dependent variables              */
      int          indep,      /* number of independent variables            */
      double       *basepoint, /* independant variable values                */
      unsigned int *rb,    
      /* rb[0] = number of blocks of dependent variables
         dependent variable j=0..depen-1 belongs to row block rb[1+j],
         if rb == NULL  each dependent variable will be considered
                        as a block of dependent variables (rb[0]=depen)      */
      unsigned int *cb,
      /* cb[0] = number of blocks of independent variables
         independent variable i=0..indep-1 belongs to column block cb[1+i]
         if cb == NULL  each independent variable will be considered
                        as a block of independent variables (cb[0]=indep)    */
      unsigned int **crs,     
      /* returned compressed row block-index storage 
           crs[ rb[0] ][ non-zero blocks of independent variables per row ]  
         crs[depen. block][0] = non-zero indep. bl. w.r.t current depen. bl.  
         crs[depen. block][1 .. crs[depen. bl.][0]] :                 
           indeces of non-zero blocks of independent variables  with
             respect to the current block of dependent variables             */
      int          *options   
      /* control options                             
                      options[0] : way of bit pattern propagation
                                 0 - automatic detection (default)
                                 1 - forward mode 
                                 2 - reverse mode
                      options[1] : test the computational graph control flow
                                 0 - safe mode (default)
                                 1 - tight mode                              */
              
      )
{
  int             rc= -1, i, j;
  unsigned int    depen_blocks, indep_blocks, *rowbl, *colbl; 
  int             ctrl_options[3];

  if (rb == NULL)
    {
      rb = myalloc1_uint(1+depen);
      rowbl = rb;
      *rowbl++ = depen;
      for (j=0; j<depen; j++)
        *rowbl++ = j;
    }
  if (cb == NULL)
    {
      cb = myalloc1_uint(1+indep);
      colbl = cb;
      *colbl++ = indep;
      for (i=0; i<indep; i++)
        *colbl++ = i;
    }
  
  rowbl = rb;
  depen_blocks = *rowbl++;
  for (j=0; j<depen; j++)
    if ( *rowbl++ >= depen_blocks )
      {
        fprintf(DIAG_OUT,"ADOL-C user error in jac_pat(...) : "
                "bad dependent block index rb[%i]=%i >= rb[0]=%i !\n", 
                j+1, *(rowbl-1), depen_blocks);
        exit(-1);
      }
  colbl = cb;
  indep_blocks = *colbl++;
  for (i=0; i<indep; i++)
    if ( *colbl++ >= indep_blocks )
      {
        fprintf(DIAG_OUT,"ADOL-C user error in jac_pat(...) : "
                "bad independent block index cb[%i]=%i >= cb[0]=%i !\n", 
                i+1, *(colbl-1), indep_blocks);
        exit(-1);
      }

  if (crs == NULL)
    {
      fprintf(DIAG_OUT,"ADOL-C user error in jac_pat(...) : "
                       "parameter crs may not be NULL !\n");
      exit(-1);
    }
  else
    for (i=0; i<depen_blocks; i++)
      crs[i] = NULL;

  if (( options[0] < 0 ) || (options[0] > 2 ))
     options[0] = 0; /* default */
  if (( options[1] < 0 ) || (options[1] > 1 ))
     options[1] = 0; /* default */
  ctrl_options[0] = options[0];
  ctrl_options[1] = options[1];
  ctrl_options[2] = 0; /* no output intern */

  rc = block_pattern( tag, depen, indep,  basepoint,
                      rb, cb, crs, ctrl_options);

  return(rc);
}

/****************************************************************************/

END_C_DECLS

/*----------------------------------------------------------------------------
 ADOL-C -- Automatic Differentiation by Overloading in C++
 File:     sparse/jacutils.c
 Revision: $Id: jacutils.c,v 1.1 2004/10/13 14:18:12 e_arnold Exp $
 Contents: This file contains definitions for the functions prototyped in
           jacutils.h. Each function is an ADOL-C straight C utility.
 
 Copyright (c) 2004
               Technical University Dresden
               Department of Mathematics
               Institute of Scientific Computing
  
 This file is part of ADOL-C. This software is provided under the terms of
 the Common Public License. Any use, reproduction, or distribution of the
 software constitutes recipient's acceptance of the terms of this license.
 See the accompanying copy of the Common Public License for more details.

 History:
   20040414 kowarz:  adaption to configure - make - make install
   20000214 olvo:    call to forward required in tight reverse mode only,
                     in safe reverse mode no basepoint is available!
                     (By the way, the return values rc are never checked!)
   20000214 olvo:    fixed bug according to Andreas Waechter from CMU
   19990308 christo: myalloc1_ushort -> myalloc1_uint
   19990308 christo: number of blocks : unsigned short -> unsigned int
   19990308 christo: bit patterns : unsigned int -> unsigned long int
   19990302 christo: new interface of jac_pat(...)
   19981203 christo: non-transponded matrices for reverse 
   19981125 christo: changed options rowblocks, colblocks
   19981115 christo: non-contiguous blocks
   19981101 christo: changed options, default=safe=0

----------------------------------------------------------------------------*/


/****************************************************************************/
/*                                                                  DEFINES */

/* define SHOWBITPATTERN to print the Jacobian bit pattern in reverse mode  */
#undef SHOWBITPATTERN


/****************************************************************************/
/*                                                       NECESSARY INCLUDES */
#include "../sparse/jacutils.h"
#include "../interfaces.h"
#include "../adalloc.h"

#include <malloc.h>

BEGIN_C_DECLS

/*****************************************************************************/
/*                                                    JACOBIAN BLOCK PATTERN */

/* ------------------------------------------------------------------------- */
int block_pattern(
   short  tag,        /* tape identification                        */
   int    depen,      /* number of dependent variables              */
   int    indep,      /* number of independent variables            */
   double *basepoint, /* independant variable values                */
   unsigned int *rowblocks,
      /* rb[0] = number of blocks of dependent variables
                 dependent variable j=0..depen-1 belongs to block rb[1+j]   */
   unsigned int *columnblocks,
      /* cb[0] = number of blocks of independent variables
                 independent variable i=0..indep-1 belongs to block cb[1+i] */
   unsigned int **crs,
      /* compressed block row storage [rb[0]][1 + indep.bl. per row]
         crs[depen. bl.][0] = non-zero indep. bl. w.r.t depen. bl.
         crs[depen. bl.][1 .. crs[depen. bl.][0]] :
         indeces of non-zero blocks of independent variables with
         respect to the current block of dependent variables                */
   int *options       /* control options                            */
      /* options[0] : way of bit pattern propagation
                      0 - automatic detection (default)
                      1 - forward mode 
                      2 - reverse mode
         options[1] : test the computational graph control flow
                      0 - safe variant (default)
                      1 - tight variant                
         options[2] : output ( 0 in 1 in 2 in 3 )
                      0 - no output (default)
                      1 - sparsity (block) patterns output
                      2 - mode specification
                      3 - seed/Jacobian matrix output as unsigned longs     */
) {
  unsigned int   depen_blocks, indep_blocks, *rowbl, *colbl;

  int                rc= 3, outp;
  char               forward_mode, tight_mode, each_depvar_is_block;
  int                i, ii, j, jj, k, k_old, bits_per_long, 
                     i_blocks_per_strip, d_blocks_per_strip;
  int                this_strip_i_bl_idx, next_strip_i_bl_idx, 
                     this_strip_d_bl_idx, next_strip_d_bl_idx;
  int                stripmined_calls, strip_idx;
  int                p_stripmine, q_stripmine, p_ind_bl_bp, q_dep_bl_bp, 
                     i_bl_idx, d_bl_idx;
  unsigned long int  value1, v;
  unsigned long int  **seed=NULL, *s, **jac_bit_pat=NULL, *jac;
  unsigned int       *rb;
  unsigned char      *indep_blocks_flags=NULL, *i_b_flags;
  double             *valuepoint=NULL;


  depen_blocks = rowblocks[0];    /*number of blocks of dependent variables  */
  indep_blocks = columnblocks[0]; /*number of blocks of independent variables*/
  rowbl = rowblocks + 1;
  colbl = columnblocks + 1;

  if ( options[0] == 0 ) {
    if ( depen >= indep ) 
      options[0] = 1; /* forward */
    else
      options[0] = 2; /* reverse */
  }

  if ( options[0] == 1 )
    forward_mode = 1;
  else
    forward_mode = 0;

  if ( options[1] == 1 )
    tight_mode = 1;
  else
    tight_mode = 0;
  
  outp = options[2];

  if ( ! forward_mode )
    valuepoint = myalloc1(depen);

  /* bit pattern parameters */
  
  /* number of bits in an unsigned long int variable */
  bits_per_long = 8 * sizeof(unsigned long int);
  /* olvo 20000214 nl: inserted explicit cast to unsigned long int */
  value1 =  (unsigned long int) 1 << (bits_per_long - 1); /* 10000....0 */

  /* =================================================== forward propagation */
  if ( forward_mode )
    {

      if (( tight_mode ) && ( basepoint == NULL ))
        {
          fprintf(DIAG_OUT, "ADOL-C error in jac_pat(...) :  supply basepoint x for tight mode.\n");
          exit(-1);
        }
      else
        if (outp > 1)
	  {
	    printf("jac_pat started in forward, ");
	    if ( tight_mode )
	      printf("tight mode !\n\n");
	    else
	      printf("safe mode !\n\n");
	  }
      
      /* indep partial derivatives for the whole Jacobian */
      
      /* number of unsigned longs to store the whole seed / Jacobian matrice */
      p_ind_bl_bp = indep_blocks / bits_per_long 
                    + ( (indep_blocks % bits_per_long) != 0 );

      /* number of unsigned longs to store the seed / Jacobian strips */
      if ( p_ind_bl_bp <= PQ_STRIPMINE_MAX )
      {
        p_stripmine = p_ind_bl_bp;
        stripmined_calls = 1;
      }
      else
        {
          p_stripmine = PQ_STRIPMINE_MAX;
          stripmined_calls = p_ind_bl_bp / PQ_STRIPMINE_MAX
                             + ( (p_ind_bl_bp % PQ_STRIPMINE_MAX) != 0 );
        }
      
      /* number of independent blocks per seed / Jacobian strip */
      i_blocks_per_strip = p_stripmine * bits_per_long;

      /* allocate memory --------------------------------------------------- */

      if ( ! (indep_blocks_flags = (unsigned char*)
              calloc(i_blocks_per_strip, sizeof(char)) ) )
        {
          fprintf(DIAG_OUT, "ADOL-C error, "__FILE__
                  ":%i : \njac_pat(...) unable to allocate %i bytes !\n",
                  __LINE__, i_blocks_per_strip*sizeof(char));
          exit(-1);
        }

      seed        = myalloc2_ulong(indep, p_stripmine);
      jac_bit_pat = myalloc2_ulong(depen, p_stripmine);
          
      if ( depen == depen_blocks )
	{
	  each_depvar_is_block = 1;
	  rb = rowbl;
	  for (j=0; j<depen; j++)
	    if ( *rb++ != j )
	      each_depvar_is_block = 0;
	}
      else
	each_depvar_is_block = 0;

      /* strip-mining : repeated forward calls ----------------------------- */

      for (strip_idx = 0; strip_idx < stripmined_calls; strip_idx++)
        {
          /* build a partition of the seed matrix (indep x indep_blocks) --- */
          /* (indep x i_blocks_per_strip) as a bit pattern                   */
	  s = seed[0]; 
	  for (i=0; i<indep; i++) 
	    for (ii=0; ii<p_stripmine; ii++) /* 2 loops if short -> int !!! */
	      *s++ = 0; /* set old seed matrix to 0 */

          this_strip_i_bl_idx = strip_idx * i_blocks_per_strip;
	  next_strip_i_bl_idx = (strip_idx+1) * i_blocks_per_strip;
          if ( next_strip_i_bl_idx > indep_blocks )
            next_strip_i_bl_idx = indep_blocks;
          v = value1; /* 10000....0 */
          
          for (i=0; i<indep; i++)
	    if ( (this_strip_i_bl_idx <= colbl[i]) 
		&& (colbl[i] < next_strip_i_bl_idx) )
	      {
		/*printf("indep var.  %i  -> block  %i  -> strip  %i\n", 
                i, colbl[i], strip_idx);*/
		/* block colbl[i] belongs to this strip */
		ii = (colbl[i] - this_strip_i_bl_idx) /  bits_per_long;
		seed[i][ii] = v >> ((colbl[i] - this_strip_i_bl_idx) %  bits_per_long);
	      }

          if (outp > 2)
            {
              printf("seed matrix (indep x indep_blocks_per_strip) as \
unsigned long ints : \n");
              for (i=0; i<indep; i++)
                {
                  for (ii=0; ii<p_stripmine; ii++)
                    printf("%10lu  ",seed[i][ii]);
                  printf("\n");
                }
              printf("\n");
            }
      
          /* bit pattern propagation by forward ---------------------------- */

          if ( tight_mode )
            rc = int_forward_tight( tag, depen, indep, p_stripmine,
                                    basepoint, seed, valuepoint, jac_bit_pat);
          else
            rc = int_forward_safe ( tag, depen, indep, p_stripmine, 
                                    seed, jac_bit_pat);

          if (outp > 2)
            {
              printf("Jacobian (depen x indep_blocks_per_strip) as unsigned \
long ints: \n");
              for (j=0; j<depen; j++)
                {
                  for (ii=0; ii<p_stripmine; ii++)
                    printf("%10lu  ",jac_bit_pat[j][ii]);
                  printf("\n");
                }
              printf("\n");
            }
    
          /* extract (block) pattern from bit patterns --------------------- */

          for (d_bl_idx = 0; d_bl_idx < depen_blocks; d_bl_idx++) 
            {
	      if ( each_depvar_is_block )
		{
		  j = d_bl_idx; /* dependent block index */
		  ii = -1; v = 0;
                    
		  jac = jac_bit_pat[j];
		  i_b_flags = indep_blocks_flags;
		  for (i_bl_idx = 0; i_bl_idx < i_blocks_per_strip; i_bl_idx++)
		    {
		      if ( !v )
			{
			  v =  value1; /* 10000....0 */
			  ii++;
			}		      
		      if ( v & jac[ii] )
			*i_b_flags = 1;
		      i_b_flags++;
		      
		      v = v >> 1;
		    }
		}
	      else
		for (j=0; j<depen; j++)
		  if ( rowbl[j] == d_bl_idx ) /* dependent block index */
		    {
		      ii = -1; v = 0;
		      
		      jac = jac_bit_pat[j];
		      i_b_flags = indep_blocks_flags;
		      for (i_bl_idx = 0; i_bl_idx < i_blocks_per_strip; i_bl_idx++)
			{
			  if ( !v )
			    {
			      v =  value1; /* 10000....0 */
			      ii++;
			    } 
			  if ( v & jac[ii] )
			    *i_b_flags = 1;
			  i_b_flags++;

			  v = v >> 1;
			}
		    }
              
              if ( strip_idx == 0 )
                k_old = 0;
              else
                k_old = crs[d_bl_idx][0];
	      k = 0;
	      i_b_flags = indep_blocks_flags;
              for (i = 0; i < i_blocks_per_strip; i++)
                k += *i_b_flags++;

	      if ((k > 0 ) || ( strip_idx == 0 ))
		{
		  if ( ! (crs[d_bl_idx] = (unsigned int*)realloc(crs[d_bl_idx], 
                (k_old+k+1)*sizeof(unsigned int))) )
		    {
		      fprintf(DIAG_OUT, "ADOL-C error, "__FILE__
			      ":%i : \njac_pat(...) unable to allocate %i bytes !\n",
			      __LINE__, (k_old+k+1)*sizeof(unsigned int));
		      exit(-1);
		    }
		  if ( strip_idx == 0 )
		    crs[d_bl_idx][0]  = 0;
		  if ( k > 0 )
		    {
		      k = crs[d_bl_idx][0] + 1;
		      i_b_flags = indep_blocks_flags;
		      for (i = 0; i < i_blocks_per_strip; i++)
			{
			  if ( *i_b_flags )
			    {
			      crs[d_bl_idx][k++] = this_strip_i_bl_idx + i;
			      *i_b_flags = 0;
			    }
			  i_b_flags++;
			}
		      /* current/total number of non-zero blocks of indep. vars. */
		      crs[d_bl_idx][0] = k - 1;
		    }
		}
    
            }
 
        } /* strip_idx */

    } /* forward */


  /* =================================================== reverse propagation */
  else
    {
      if (outp > 1)
	{
	  printf("jac_pat started in reverse, ");
	  if ( tight_mode )
	      printf("tight mode !\n\n");
	    else
	      printf("safe mode !\n\n");
	}

      /* depen weight vectors for the whole Jacobian */
      
      /* number of unsigned longs to store the whole seed / Jacobian matrice */
      q_dep_bl_bp = depen_blocks / bits_per_long 
                    + ( (depen_blocks % bits_per_long) != 0 );
      
       /* number of unsigned longs to store the seed / Jacobian strips */
      if ( q_dep_bl_bp <= PQ_STRIPMINE_MAX )
      {
        q_stripmine = q_dep_bl_bp;
        stripmined_calls = 1;
      }
      else
        {
          q_stripmine = PQ_STRIPMINE_MAX;
          stripmined_calls = q_dep_bl_bp / PQ_STRIPMINE_MAX
                             + ( (q_dep_bl_bp % PQ_STRIPMINE_MAX) != 0 );
        }
      
      /* number of dependent blocks per seed / Jacobian strip */
      d_blocks_per_strip = q_stripmine * bits_per_long;

      /* allocate memory --------------------------------------------------- */
      if ( ! (indep_blocks_flags = (unsigned char*)calloc(indep_blocks, 
              sizeof(unsigned char)) ) )
        {
          fprintf(DIAG_OUT, "ADOL-C error, "__FILE__
                  ":%i : \njac_pat(...) unable to allocate %i bytes !\n",
                  __LINE__, indep_blocks*sizeof(unsigned char));
          exit(-1);
        }

      seed        = myalloc2_ulong(q_stripmine, depen);
      jac_bit_pat = myalloc2_ulong(q_stripmine, indep);


      /* olvo 20000214: call to forward required in tight mode only,
         in safe mode no basepoint available! */
      if ( tight_mode )
      { if ( basepoint == NULL )
        { fprintf(DIAG_OUT, "ADOL-C error in jac_pat(..) :  ");
          fprintf(DIAG_OUT, "no basepoint x for tight mode supplied.\n");
          exit(-1);
        }

        rc = zos_forward(tag, depen, indep, 1, basepoint, valuepoint);
      }

      /* strip-mining : repeated reverse calls ----------------------------- */

      for (strip_idx = 0; strip_idx < stripmined_calls; strip_idx++)
        {
          /* build a partition of the seed matrix (depen_blocks x depen)     */
          /* (d_blocks_per_strip x depen) as a bit pattern                   */
	  s = seed[0];
          for (jj=0; jj<q_stripmine; jj++) /* 2 loops if short -> int !!! */
            for (j=0; j<depen; j++)
	      *s++ = 0; /* set old seed matrix to 0 */
          
          this_strip_d_bl_idx = strip_idx * d_blocks_per_strip;
	  next_strip_d_bl_idx = (strip_idx+1) * d_blocks_per_strip;
	  if ( next_strip_d_bl_idx > depen_blocks )
            next_strip_d_bl_idx = depen_blocks;
          v = value1; /* 10000....0 */

          for (j=0; j<depen; j++)
	    if ( (this_strip_d_bl_idx <= rowbl[j]) 
		&& (rowbl[j] < next_strip_d_bl_idx) )
	      {
		/*printf("depen var.  %i  -> block  %i  -> strip  %i\n",j, rowbl[j], 
        strip_idx);*/
		/* block rowbl[j] belongs to this strip */
      jj = (rowbl[j] - this_strip_d_bl_idx) /  bits_per_long;
      seed[jj][j] = v >> ((rowbl[j] - this_strip_d_bl_idx) % bits_per_long);
              }

          if (outp > 2)
            {
              printf("seed matrix (depen_blocks_per_strip x depen) ");
              printf("as unsigned long ints: \n");
              for (jj=0; jj<q_stripmine; jj++)
                {
                  for (j=0;j<depen;j++)
                    printf("%10lu  ",seed[jj][j]);
                  printf("\n");
                }
              printf("\n");
            }

          /* bit pattern propagation by reverse ---------------------------- */

          if ( tight_mode )
            rc = int_reverse_tight( tag, depen, indep, q_stripmine,
                                    seed, jac_bit_pat);
          else
            rc = int_reverse_safe ( tag, depen, indep, q_stripmine,
                                    seed, jac_bit_pat);

          if (outp > 2)
            {
              printf("Jacobian (depen_blocks_per_strip x indep) ");
              printf("as unsigned long ints: \n");
              for (jj=0; jj<q_stripmine; jj++)
                {
                  for (i=0;i<indep;i++)
                    printf("%10lu  ",jac_bit_pat[jj][i]);
                  printf("\n");
                }
              printf("\n");
            }

          /* extract (block) pattern from bit patterns --------------------- */
          
#ifdef SHOWBITPATTERN
          if (outp > 2)
            printf("Jacobian (depen_blocks_per_strip x indep) ");
            printf("as a bit pattern: \n");
#endif

          jj = -1;
          v = 0; 
          for (d_bl_idx = this_strip_d_bl_idx; 
               d_bl_idx < next_strip_d_bl_idx; d_bl_idx++)
            {
              if ( !v )
                {
                  v =  value1; /* 10000....0 */
                  jj++;
                }
              jac = jac_bit_pat[jj];
              for (i=0; i<indep; i++)
                {
                  if ( v & *jac++ )
                    {
                      indep_blocks_flags[colbl[i]] = 1;
#ifdef SHOWBITPATTERN
                      if (outp > 2)
                        printf("X ");
#endif
                    }
#ifdef SHOWBITPATTERN
                  else
                    if (outp > 2)
                      printf("%c ",183);
#endif
                }
              
              v = v >> 1;

#ifdef SHOWBITPATTERN              
              if (outp > 2)
                printf("\n");
#endif          
              k=0;
	      i_b_flags = indep_blocks_flags;
              for (i=0; i<indep_blocks; i++)
                k += *i_b_flags++;
              
              if ( ! (crs[d_bl_idx] = (unsigned int*)malloc((k+1)*sizeof(unsigned int))) )
                {
                  fprintf(DIAG_OUT, "ADOL-C error, "__FILE__
                          ":%i : \njac_pat(...) unable to allocate %i bytes !\n",
                          __LINE__, (k+1)*sizeof(unsigned int));
                  exit(-1);
                }
              crs[d_bl_idx][0] = k; /* number of non-zero indep. blocks */
              k=1;
	      i_b_flags = indep_blocks_flags;
              for (i=0; i<indep_blocks; i++)
		{
		  if ( *i_b_flags )
		    {
		      crs[d_bl_idx][k++] = i;
		      *i_b_flags = 0;
		    }
		  i_b_flags++;
		}
            }

        } /* strip_idx */

    } /* reverse */

  if ( ! forward_mode )
    {
      free((char*)valuepoint); valuepoint=NULL;
    }
  free((char*)*seed); free((char*)seed); seed=NULL;
  free((char*)*jac_bit_pat); free((char*)jac_bit_pat); jac_bit_pat=NULL;
  free((char*)indep_blocks_flags); indep_blocks_flags=NULL;

  if (outp > 0)
    {
      printf("\nJacobian Block Pattern :\n");
      for (j=0; j<depen_blocks; j++)
        {
          printf("dep. bl[%i], %i nzbl :", j, crs[j][0]);
          for (i=1; i<=crs[j][0]; i++)
            printf("%i ", crs[j][i]);
          printf("\n");
        }
      printf("\n");
    }

  return(rc);
}


/*****************************************************************************/
/*                                               MEMORY MANAGEMENT UTILITIES */

/* ------------------------------------------------------------------------- */
unsigned int *myalloc1_uint(int m)
{
  unsigned int *A = (unsigned int*)malloc(m*sizeof(unsigned int));
  if (A == NULL){
    fprintf(DIAG_OUT, "ADOL-C error, "__FILE__
            ":%i : \nmyalloc1_ushort cannot allocate %i bytes\n",
            __LINE__, m*sizeof(unsigned int));
    exit (-1);
  } /* endif */
  return A;
}


/* ------------------------------------------------------------------------- */
unsigned long int *myalloc1_ulong(int m)
{
  unsigned long int *A = (unsigned long int*)  calloc(m,sizeof(unsigned long int));
  if (A == NULL){
    fprintf(DIAG_OUT, "ADOL-C error, "__FILE__
            ":%i : \nmyalloc1_ulong cannot allocate %i bytes\n",
            __LINE__, m*sizeof(unsigned long int));
    exit (-1);
  } /* endif */
  return A;
}


/* ------------------------------------------------------------------------- */
unsigned long int **myalloc2_ulong(int m,int n)
{
  unsigned long int *Adum = (unsigned long int*)  calloc(m*n,sizeof(unsigned long int));
  unsigned long int **A   = (unsigned long int**) calloc(m,sizeof(unsigned long int*));
  int i;
  if (Adum == NULL){
    fprintf(DIAG_OUT, "ADOL-C error, "__FILE__
            ":%i : \nmyalloc2_ulong cannot allocate %i bytes\n",
            __LINE__, m*n*sizeof(unsigned long int));
    exit (-1);
  } /* endif */
  if (A == NULL){
    fprintf(DIAG_OUT, "ADOL-C error, "__FILE__
            ":%i : \nmyalloc2_ulong cannot allocate %i bytes\n",
            __LINE__, m*sizeof(unsigned long int*));
    exit (-1);
  } /* endif */
  for(i=0;i<m;i++)
    {
      A[i] = Adum;
      Adum += n;
    }
  return A;

  /* To deallocate an array set up by   A = myalloc2_ulong(m,n)   */
  /*    use  free((char*)*A); free((char*)A);  in that order      */

 }


/****************************************************************************/
/*                                                               THAT'S ALL */

END_C_DECLS

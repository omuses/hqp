/*
   --------------------------------------------------------------
   File jacpatexam.C  of ADOL-C version 1.8.2     as of Mar/10/98
   --------------------------------------------------------------

   Example: Jacobian (Block) Pattern Example of ADOL-C 
  
   Last changed:
    990308 christo: number of blocks : unsigned short -> unsigned int    
    990302 christo: new interface of jac_pat(...)
    981201 olvo:    new headers
    981125 christo: new call to jac_pat(..)
    981113 christo: new example
  
   -------------------------------------------------------------- 
*/

/****************************************************************************/
/*                                                                 INCLUDES */
#include "../../../adolc/adolc.h"
#include "../clock/myclock.h"

#include <string.h>
#include <math.h>

#include <iostream>
using namespace std;

/****************************************************************************/
/*                                                                  DEFINES */

#define TAG 11


/****************************************************************************/
/*                                                     EVALUATION FUNCTIONS */

/*--------------------------------------------------------------------------*/
const unsigned int N = 5, M = 6;

void eval_small(short tag, double *xp, double *yp)
{
  unsigned int i,j;

  trace_on(tag); 

  adouble *x,*y;
  x = new adouble[N];
  y = new adouble[M];
  for (i=0;i<N;i++)
    x[i] <<= xp[i];

  int PD1B = __LINE__;

  y[0] = pow(x[0],1) + pow(x[1],2) + pow(x[2],3);
  y[1] = x[0]        * x[1]        / x[2];
  y[2] =                                          asin(x[3]*0.1); 
  y[3] = sqrt(x[0])  + sqrt(x[1])  + sqrt(x[2]) + sqrt(x[3])      + sqrt(x[4]);
  y[4] =                                          log(x[3])       + exp(x[4]);
  y[5] =                                          cos(x[3])       - sin(x[4]);

  int PD1E = __LINE__;

  for (j=0;j<M;j++)
    y[j] >>= yp[j];

  delete []x;
  delete []y;

  trace_off(1); // force a numbered tape file to be written 

  cout << "\nproblem definition in  "<<__FILE__<<",  lines  "<<PD1B<<" - "<<PD1E<<"\n";  
}


/*--------------------------------------------------------------------------*/
const unsigned int NM = 961; // PQ_STRIPMINE_MAX * 8*sizeof(unsigned long int) + 1

void eval_arrow_like_matrix(short tag, double *xp, double *yp)
{
  unsigned int i,j;
  
  trace_on(tag); 

  adouble *x,*y;
  x = new adouble[NM];
  y = new adouble[NM];
  for (i=0;i<NM;i++)
    x[i] <<= xp[i];

  int PD2B = __LINE__;

  for (i=0;i<NM;i++)
    {
      /* dense diagonal and dense last column*/
      y[i] = cos(x[i]) + sin(x[NM-1]);
    }
  for (i=0;i<NM;i++)
    /* dense last row */
    y[NM-1] += sin(x[i]);

  int PD2E = __LINE__;

  for (j=0;j<NM;j++)
    y[j] >>= yp[j];

  delete []x;
  delete []y;

  trace_off(1); // force a numbered tape file to be written 

  cout << "\nproblem definition in  "<<__FILE__<<",  lines  "<<PD2B<<" - "<<PD2E<<"\n";
}


/****************************************************************************/
/*                                                             MAIN PROGRAM */
int main(void) 
{
  short  tag;
  int    ret_c = -1, choice;
  int    oper,buf_size,maxlive,deaths;
  unsigned int depen, indep;
  int    tape_stats[11] ;
  unsigned int    i, j, minnz, maxnz, nz, nzref, nz_rel;
  double z1, z2, t0, t1, t2, t3, t4, t5=0.0, t6=0.0; 
  char   mode[8], outp, full_jac;
  int    precision, width;
  
/*--------------------------------------------------------------------------*/
/*            variables needed for the Jacobian (block) pattern exploration */

  unsigned int  depenbl=0;   // number of blocks of dependent variables
  unsigned int  indepbl=0;   // number of blocks of independent variables 
  unsigned int  *rowbl=NULL; // starting index of each block of 
                             // dependent variables
  unsigned int  *colbl=NULL; // starting index of each block of 
                             // independent variables
  unsigned int  **jac_block_pat=NULL; // compressed block row storage
  double        *base, *value;
  double        basepoint ;
  int           ctrl[2];

  cout << "----------------------------------------------------------------\n";
  cout << "\n                Jacobian (Block) Pattern Example\n\n";
  cout << "Tape identification tag ( [-4..-1] for standart examples ) :  ?\b" ;
  cin  >> choice ;

  cout << "\n\nOutput Jacobian (block) pattern? (y/n)  ?\b";
  cin >> outp;

  cout << "\n\nCompare with the full Jacobian calculation? (y/n)  ?\b";
  cin >> full_jac;
  if (( full_jac == 'y' ) || ( full_jac == 'Y'))
    full_jac = 1;
  else
    full_jac = 0;

  cout << "----------------------------------------------------------------\n";

/*--------------------------------------------------------------------------*/

  
  if ( choice < 0 ) // Take function in the "eval(...)" routines -------------
    {
      if ( choice > -4 )
        {
          base = new double[N];
          for (i=0;i<N;i++)
            base[i] = i+1;
          
          value = new double[M];
          
          tag = TAG;
          eval_small(tag, base, value);
          
          cout << "\n\nCreated ADOL-C tape with identification tag " << tag << ".\n\n";
          cout.flush();
      
          if ( choice == -1 )//-----------------------------------------------
            {
              int BD1B = __LINE__;

              // variables in contiguous blocks
              depenbl = 4;
              rowbl = new unsigned int[1+M];
              rowbl[0] = depenbl;
              rowbl[1+0] = 0;
              rowbl[1+1] = 0;
              rowbl[1+2] = 1;
              rowbl[1+3] = 2;
              rowbl[1+4] = 3;
              rowbl[1+5] = 3;
              
              indepbl = 3;
              colbl = new unsigned int[1+N];
              colbl[0] = indepbl;
              colbl[1+0] = 0;
              colbl[1+1] = 0;
              colbl[1+2] = 0;
              colbl[1+3] = 1;
              colbl[1+4] = 2;

              int BD1E = __LINE__;
              cout << "variables (block) definition in  "<<__FILE__<<",  lines  "<<BD1B<<" - "<<BD1E<<"\n";
            }
          
          if ( choice == -2 )//-----------------------------------------------
            {
              int BD2B = __LINE__;

              // variables in non-contiguous blocks
              depenbl = 3;
              rowbl = new unsigned int[1+M];
              rowbl[0] = depenbl;
              rowbl[1+0] = 0;
              rowbl[1+1] = 1;
              rowbl[1+2] = 1;
              rowbl[1+3] = 2;
              rowbl[1+4] = 1;
              rowbl[1+5] = 0;
              
              indepbl = 3;
              colbl = new unsigned int[1+N];
              colbl[0] = indepbl;
              colbl[1+0] = 0;
              colbl[1+1] = 1;
              colbl[1+2] = 2;
              colbl[1+3] = 1;
              colbl[1+4] = 0;

              int BD2E = __LINE__;
              cout << "variables (block) definition in  "<<__FILE__<<",  lines  "<<BD2B<<" - "<<BD2E<<"\n";
            }
          
          if ( choice == -3 )//-----------------------------------------------
            {
              int BD3B = __LINE__;

              // each variable as a separate block  
              depenbl = M;
              rowbl = new unsigned int[1+M];
              rowbl[0] = depenbl;
              for (j=0;j<M;j++)
                rowbl[j+1] = j; 

              indepbl = N;
              colbl = new unsigned int[1+N];
              colbl[0] = indepbl;
              for (i=0;i<N;i++)
                colbl[i+1] = i;

              int BD3E = __LINE__;
              cout << "variables (block) definition in  "<<__FILE__<<",  lines  "<<BD3B<<" - "<<BD3E<<"\n";
            }
        }
      else // ( choice == -4 ) -----------------------------------------------
        {
          int BD4B = __LINE__;

          base = new double[NM];
          for (i=0;i<NM;i++)
            base[i] = i;
          
          value = new double[NM];
          
          tag = TAG;
          eval_arrow_like_matrix(tag, base, value);
          
          cout << "\n\nCreated ADOL-C tape with identification tag " << tag << ".\n\n";
          cout.flush();

          // each variable as a separate block
          depenbl = 0;          
          rowbl = NULL;

          indepbl = 0;
          colbl = NULL;          

          int BD4E = __LINE__;
          cout << "variables (block) definition in  "<<__FILE__<<",  lines  "<<BD4B<<" - "<<BD4E<<"\n";
        }

      tapestats(tag,tape_stats);              /* Reading of tape statistics */
      indep    = tape_stats[0];
      depen    = tape_stats[1];
      buf_size = tape_stats[4];
      oper     = tape_stats[5];
      deaths   = tape_stats[3];
      maxlive  = tape_stats[2];
  
      cout<<"\nTape " << tag << ": \n";
      cout<<"  independents   "<<indep<<"\n";
      cout<<"  dependents     "<<depen<<"\n";
      cout<<"  operations     "<<oper<<"\n";
      cout<<"  buffer size    "<<buf_size<<"\n";
      cout<<"  maxlive        "<<maxlive<<"\n";
      cout<<"  valstack size  "<<deaths<<"\n";
      cout<<"\n";
      cout<<"dependent blocks   "<<depenbl<<"\n";
      cout<<"independent blocks "<<indepbl<<"\n\n";
      
      if (( depenbl == 0 ) || ( depenbl == depen ))
        cout << "=> each dependent variable as a block of dependent variable(s)\n";
      if (( indepbl == 0 ) || ( indepbl == indep ))
        cout << "=> each independent variable as a block of independent variable(s)\n";
    }

  else // ( choice >= 0 ) : Take a written tape ------------------------------
    {
      tag = choice;

      cout << "\nproblem definition in  tape "<<tag<<"\n";

      tapestats(tag,tape_stats) ;
      indep    = tape_stats[0] ;
      depen    = tape_stats[1] ;
      buf_size = tape_stats[4];
      oper     = tape_stats[5];
      deaths   = tape_stats[3];
      maxlive  = tape_stats[2];
      
      cout<<"\nTape " << tag << ": \n";
      cout<<"  independents   "<<indep<<"\n";
      cout<<"  dependents     "<<depen<<"\n";
      cout<<"  operations     "<<oper<<"\n";
      cout<<"  buffer size    "<<buf_size<<"\n";
      cout<<"  maxlive        "<<maxlive<<"\n";
      cout<<"  valstack size  "<<deaths<<"\n";

      cout << "\n\nbasepoint[0.."<< indep <<"] = ?\b";
      cin  >> basepoint;

      base  = new double[indep];
      value = new double[depen];

      for (i=0;i<indep;i++)
        base[i] = basepoint;

#define BDTB __LINE__
      // each variable as a separate block
      depenbl = 0;
      rowbl = NULL;

      indepbl = 0;
      colbl = NULL;
#define BDTE __LINE__
      cout << "\nvariables (block) definition in  "<<__FILE__<<",  lines  "<<BDTB<<" - "<<BDTE<<"\n";

      cout<<"\n";
      cout<<"dependent blocks   "<<depenbl<<"\n";
      cout<<"independent blocks "<<indepbl<<"\n\n";

      if (( depenbl == 0 ) || ( depenbl == depen ))
        cout << "=> each dependent variable as a block of dependent variable(s)\n";
      if (( indepbl == 0 ) || ( indepbl == indep ))
        cout << "=> each independent variable as a block of independent variable(s)\n";
    }

  //tape_doc(tag,depen,indep,base,value); // write a tape into a tex-file


/*==========================================================================*/
/*Jacobian (block) pattern automatic: if depen >= indep forward else reverse*/

  if ( depen >= indep )
    strcpy(mode,"forward");
  else
    strcpy(mode,"reverse");
  cout << "\nJacobian Block Pattern by automatic mode = " << mode << ", safe ...\n";

  if ( rowbl != NULL )
    depenbl = rowbl[0];
  else
    depenbl = depen;
  if ( colbl != NULL )
    indepbl = colbl[0];
  else
    indepbl = indep;

  jac_block_pat = new unsigned int* [depenbl];
  
  ctrl[0] = 0; // automatic mode choice
  ctrl[1] = 0; // safe

  z1 = myclock() ;
  ret_c = jac_pat( tag, depen, indep,  base, 
                   rowbl, colbl, jac_block_pat, ctrl);
  z2 = myclock() ;

  if (( outp == 'y' ) || ( outp == 'Y'))
    {
      for (i=0;i<depenbl;i++) 
        {
          cout <<"depen. blocks["<< i <<"], "<< jac_block_pat[i][0] <<" non-zero indep. blocks :\n";
          for (j=1;j<=jac_block_pat[i][0];j++)
            cout << jac_block_pat[i][j] <<"  ";
          cout <<"\n";
        }
      cout.flush();
    }

  t0 = z2 - z1;   

  nz = 0; minnz = indepbl; maxnz = 0;
  for (i=0;i<depenbl;i++) 
    {
      nz += jac_block_pat[i][0]; 
      if ( jac_block_pat[i][0] < minnz )
        minnz = jac_block_pat[i][0];
      if ( jac_block_pat[i][0] > maxnz )
        maxnz = jac_block_pat[i][0];  
    }
  nz_rel = (int) ceil (100*nz / ((double)depenbl*indepbl));
  cout << nz << " non-zero Jacobian elements (blocks) of total " << depenbl*indepbl << " elements (blocks) <= " << nz_rel << "%\n";
  cout << "min " << minnz << " non-zeros per row;    max " << maxnz << " non-zeros per row;\n"; 
  nzref = nz;

  for (i=0;i<depenbl;i++) 
    delete[] jac_block_pat[i];
  delete[] jac_block_pat;jac_block_pat=NULL;


/*--------------------------------------------------------------------------*/
/*                               Jacobian (block) pattern by forward, tight */

  cout << "\nJacobian Block Pattern by forward, tight ...\n";

  if ( rowbl != NULL )
    depenbl = rowbl[0];
  else
    depenbl = depen;
  if ( colbl != NULL )
    indepbl = colbl[0];
  else
    indepbl = indep;

  jac_block_pat = new unsigned int* [depenbl];

  ctrl[0] = 1; // forward
  ctrl[1] = 1; // tight

  z1 = myclock() ;
  ret_c = jac_pat( tag, depen, indep,  base, 
                   rowbl, colbl, jac_block_pat, ctrl);
  z2 = myclock() ;
  
  if (( outp == 'y' ) || ( outp == 'Y'))
    {
      for (i=0;i<depenbl;i++) 
        {
          cout <<"depen. blocks["<< i <<"], "<< jac_block_pat[i][0] <<" non-zero indep. blocks :\n";
          for (j=1;j<=jac_block_pat[i][0];j++)
            cout << jac_block_pat[i][j] <<"  ";
          cout <<"\n";
        }
      cout.flush();
    }

  t1 = z2 - z1; 

  nz = 0; minnz = indepbl; maxnz = 0;
  for (i=0;i<depenbl;i++) 
    {
      nz += jac_block_pat[i][0]; 
      if ( jac_block_pat[i][0] < minnz )
        minnz = jac_block_pat[i][0];
      if ( jac_block_pat[i][0] > maxnz )
        maxnz = jac_block_pat[i][0];  
    }
  nz_rel = (int) ceil (100*nz / ((double)depenbl*indepbl));
  cout << nz << " non-zero Jacobian elements (blocks) of total " << depenbl*indepbl << " elements (blocks) <= " << nz_rel << "%\n";
  cout << "min " << minnz << " non-zeros per row;    max " << maxnz << " non-zeros per row;\n"; 
  if ( nz != nzref )
    cout << "\n\n!!! This method found a different number of non-zeros !!!\n\n"; 

  for (i=0;i<depenbl;i++) 
    delete[] jac_block_pat[i];
  delete[] jac_block_pat;jac_block_pat=NULL;


/*--------------------------------------------------------------------------*/
/*                                Jacobian (block) pattern by forward, safe */
  
  cout << "\nJacobian Block Pattern by forward, safe ...\n";

  if ( rowbl != NULL )
    depenbl = rowbl[0];
  else
    depenbl = depen;
  if ( colbl != NULL )
    indepbl = colbl[0];
  else
    indepbl = indep;

  jac_block_pat = new unsigned int* [depenbl];

  ctrl[0] = 1; // forward
  ctrl[1] = 0; // safe

  z1 = myclock() ;
  ret_c = jac_pat( tag, depen, indep,  base, 
                   rowbl, colbl, jac_block_pat, ctrl);
  z2 = myclock() ;

  if (( outp == 'y' ) || ( outp == 'Y'))
    {
      for (i=0;i<depenbl;i++) 
        {
          cout <<"depen. blocks["<< i <<"], "<< jac_block_pat[i][0] <<" non-zero indep. blocks :\n";
          for (j=1;j<=jac_block_pat[i][0];j++)
            cout << jac_block_pat[i][j] <<"  ";
          cout <<"\n";
        }
      cout.flush();
    }

  t2 = z2 - z1; 

  nz = 0; minnz = indepbl; maxnz = 0;
  for (i=0;i<depenbl;i++) 
    {
      nz += jac_block_pat[i][0]; 
      if ( jac_block_pat[i][0] < minnz )
        minnz = jac_block_pat[i][0];
      if ( jac_block_pat[i][0] > maxnz )
        maxnz = jac_block_pat[i][0];  
    }
  nz_rel = (int) ceil (100*nz / ((double)depenbl*indepbl));
  cout << nz << " non-zero Jacobian elements (blocks) of total " << depenbl*indepbl << " elements (blocks) <= " << nz_rel << "%\n";
  cout << "min " << minnz << " non-zeros per row;    max " << maxnz << " non-zeros per row;\n"; 
  if ( nz != nzref )
    cout << "\n\n!!! This method found a different number of non-zeros !!!\n\n"; 
  
  for (i=0;i<depenbl;i++) 
    delete[] jac_block_pat[i];
  delete[] jac_block_pat;jac_block_pat=NULL;


/*--------------------------------------------------------------------------*/
/*                               Jacobian (block) pattern by reverse, tight */

  cout << "\nJacobian Block Pattern by reverse, tight ...\n";

  if ( rowbl != NULL )
    depenbl = rowbl[0];
  else
    depenbl = depen;
  if ( colbl != NULL )
    indepbl = colbl[0];
  else
    indepbl = indep;

  jac_block_pat = new unsigned int* [depenbl];

  ctrl[0] = 2; // reverse
  ctrl[1] = 1; // tight

  z1 = myclock() ;
  ret_c = jac_pat( tag, depen, indep,  base, 
                   rowbl, colbl, jac_block_pat, ctrl);
  z2 = myclock() ;

  if (( outp == 'y' ) || ( outp == 'Y'))
    {
      for (i=0;i<depenbl;i++) 
        {
          cout <<"depen. blocks["<< i <<"], "<< jac_block_pat[i][0] <<" non-zero indep. blocks :\n";
          for (j=1;j<=jac_block_pat[i][0];j++)
            cout << jac_block_pat[i][j] <<"  ";
          cout <<"\n";
        }
      cout.flush();
    }

  t3 = z2 - z1;

  nz = 0; minnz = indepbl; maxnz = 0;
  for (i=0;i<depenbl;i++) 
    {
      nz += jac_block_pat[i][0]; 
      if ( jac_block_pat[i][0] < minnz )
        minnz = jac_block_pat[i][0];
      if ( jac_block_pat[i][0] > maxnz )
        maxnz = jac_block_pat[i][0];  
    }
  nz_rel = (int) ceil (100*nz / ((double)depenbl*indepbl));
  cout << nz << " non-zero Jacobian elements (blocks) of total " << depenbl*indepbl << " elements (blocks) <= " << nz_rel << "%\n";
  cout << "min " << minnz << " non-zeros per row;    max " << maxnz << " non-zeros per row;\n"; 
  if ( nz != nzref )
    cout << "\n\n!!! This method found a different number of non-zeros !!!\n\n"; 

  for (i=0;i<depenbl;i++) 
    delete[] jac_block_pat[i];
  delete[] jac_block_pat;jac_block_pat=NULL;


/*--------------------------------------------------------------------------*/
/*                                Jacobian (block) pattern by reverse, safe */

  cout << "\nJacobian Block Pattern by reverse, safe ...\n";

  if ( rowbl != NULL )
    depenbl = rowbl[0];
  else
    depenbl = depen;
  if ( colbl != NULL )
    indepbl = colbl[0];
  else
    indepbl = indep;

  jac_block_pat = new unsigned int* [depenbl];

  ctrl[0] = 2; // reverse
  ctrl[1] = 0; // safe

  z1 = myclock() ;
  ret_c = jac_pat( tag, depen, indep,  base, 
                   rowbl, colbl, jac_block_pat, ctrl);
  z2 = myclock() ;

  if (( outp == 'y' ) || ( outp == 'Y'))
    {
      for (i=0;i<depenbl;i++) 
        {
          cout <<"depen. blocks["<< i <<"], "<< jac_block_pat[i][0] <<" non-zero indep. blocks :\n";
          for (j=1;j<=jac_block_pat[i][0];j++)
            cout << jac_block_pat[i][j] <<"  ";
          cout <<"\n";
        }
      cout.flush();
    }

  t4 = z2 - z1;

  nz = 0; minnz = indepbl; maxnz = 0;
  for (i=0;i<depenbl;i++) 
    {
      nz += jac_block_pat[i][0]; 
      if ( jac_block_pat[i][0] < minnz )
        minnz = jac_block_pat[i][0];
      if ( jac_block_pat[i][0] > maxnz )
        maxnz = jac_block_pat[i][0];  
    }
  nz_rel = (int) ceil (100*nz / ((double)depenbl*indepbl));
  cout << nz << " non-zero Jacobian elements (blocks) of total " << depenbl*indepbl << " elements (blocks) <= " << nz_rel << "%\n";
  cout << "min " << minnz << " non-zeros per row;    max " << maxnz << " non-zeros per row;\n"; 
  if ( nz != nzref )
    cout << "\n\n!!! This method found a different number of non-zeros !!!\n\n"; 

  for (i=0;i<depenbl;i++) 
    delete[] jac_block_pat[i];
  delete[] jac_block_pat;jac_block_pat=NULL;


  /* full Jacobian evaluation -----------------------------------------------*/

  if ( full_jac )
    {
      /*---------------------------------------------------------------------*/
      /*                        variables needed for the evaluation routines */

      double **Jac = new double*[depen];
      for (i=0;i<depen;i++)
        Jac[i] = new double[indep];
      double **I = new double*[indep];
      for (i=0;i<indep;i++)
        {
          I[i] = new double[indep];
          for (j=0;j<indep;j++)
            I[i][j] = 0.0;
          I[i][i] = 1.0;
        }


      /*---------------------------------------------------------------------*/
      /*                full Jacobian evaluation by forward, no strip-mining */

      cout << "\nFull Jacobian evaluation by forward(..), no \"strip-mining\" ...\n";
      
      z1 = myclock() ;
      
      forward( tag, depen, indep, indep, base, I, value, Jac);
      
      z2 = myclock() ;
      
      t5 = z2 - z1;
      

      /*---------------------------------------------------------------------*/
      /*    full Jacobian evaluation by the jacobian driver, no strip-mining */

      cout << "\nFull Jacobian evaluation by the jacobian driver, no \"strip-mining\" ...\n";
      
      z1 = myclock() ;
      
      jacobian( tag, depen, indep, base, Jac);

      z2 = myclock() ;

      t6 = z2 - z1;

    }


/*--------------------------------------------------------------------------*/
                                                       /* output of timings */
  width = 8;
  precision = 2;

  cout.setf(ios::fixed,ios::floatfield);
  cout.setf(ios::right,ios::adjustfield);
  cout.precision(precision);

  cout << "\n\n----------------------------------------------------------------\n";
  cout << "\n\nTime to explore the Jacobian (Block) Pattern by :\n\n";
  cout << " automatic mode = " << mode << ", safe      :  ";
  cout.width(width);
  cout << t0 << " sec.\n";
  cout << " forward, tight                      :  ";
  cout.width(width);
  cout << t1 << " sec.\n";
  cout << " forward, safe                       :  ";
  cout.width(width);
  cout << t2 << " sec.\n";
  cout << " reverse, tight                      :  ";
  cout.width(width);
  cout << t3 << " sec.\n";
  cout << " reverse, safe                       :  ";
  cout.width(width);
  cout << t4 << " sec.\n\n";
  if ( full_jac )
    {
      cout << " full Jacobian evaluation, forward   :  ";
      cout.width(width);
      cout << t5 << " sec.\n";
      cout << " full Jacobian evaluation, jacobian  :  ";
      cout.width(width);
      cout << t6 << " sec.\n";
    }

  if ( ! ( t0 && t1 && t2 && t3 && t4 ) )
    cout <<"\n! Zero timing due to low problem dimension.\n";

  cout <<"\nOK, done.\n";
  cout.flush();

  return 1;  
}

      **********************RKSUITE read.me file *************************
                   ****************************************


                  RKSUITE Release 1.0  November 1991

                             written by

          R.W. Brankin (*),  I. Gladwell(**),  and  L.F. Shampine (**)

 (*)  Numerical Algorithms Group Ltd.
      Wilkinson House
      Jordan Hill Road
      Oxford OX2 8DR
      U.K.
      email: richard@nag.co.uk
             na.brankin@na-net.ornl.gov
      International phone: + 44 865 511245
      International fax: + 44 865 310139

 (**) Department of Mathematics
      Southern Methodist University
      Dallas, Texas 75275
      U.S.A.
      email: h5nr1001@vm.cis.smu.edu
      U.S. phone: (214) 692-2542
      U.S. fax: (214) 692-4138


RKSUITE is a suite of codes based on Runge-Kutta methods for the numerical
solution of the initial value problem for a first order system of ordinary
differential equations.  It is the result of decades of research and
development of such methods and software by the authors.  It supersedes
some very widely used codes written by the authors and their coauthors,
namely, the RKF45 code that is available in several books and its descendant
DDERKF in the SLATEC library, and D02PAF and the associated codes in the NAG
Fortran library.

RKSUITE is being made available free of charge to the scientific community as
a public service.  It is expected that anyone making substantial use of the
software will acknowledge this use, and in particular, give a proper citation
in any publications resulting from this use.  A suitable reference is:

           R.W. Brankin, I. Gladwell, and L.F. Shampine, RKSUITE: a suite
           of Runge-Kutta codes for the initial value problem for ODEs,
           Softreport 92-S1, Department of Mathematics, Southern Methodist 
           University, Dallas, Texas, U.S.A, 1992.
         
The authors have tested the codes on a variety of problems, computers, and
compilers.  The codes are believed to perform correctly on problems for
which they were designed.  Of course, errors are possible in a software 
project of this size.  The authors assume no responsibility for the 
consequences of errors resulting from the use of this free software.  They
would greatly appreciate notification of unsatisfactory performance that 
might indicate the presence of a bug.  Constructive criticism would also be
much appreciated.

Future releases are planned that will add capabilities to the suite and 
correct any errors that might be discovered.  
                      
                      ============================================
                      ============================================
                      YOU SHOULD READ THE DOCUMENTATION CAREFULLY.
                      ============================================
                      ============================================


                      ===================================================
                      ===================================================
                      THE EASIEST WAY TO SOLVE A PROBLEM IS OFTEN TO EDIT 
                      ONE OF THE TEMPLATES PROVIDED IN THIS DIRECTORY.
                      ===================================================
                      ===================================================


                             Installation Details

All machine-dependent aspects of the suite have been isolated in the 
subroutine ENVIRN found in the rksuite.for file.  Some environmental
parameters must be specified in this subroutine.  The values in the
distribution version are those appropriate to the IEEE arithmetic 
standard.  They must be altered, if necessary, to values appropriate to 
the computing system you are using before calling the codes of the suite.
If the IEEE arithmetic standard values are not appropriate for your
system, appropriate values can be often be obtained by calling routines
named in the Comments of ENVIRN.

         ================================================================
         ================================================================
         TO MAKE SURE THAT YOU SPECIFY THESE MACHINE AND INSTALLATION 
         DEPENDENT QUANTITIES PROPERLY, WHEN THE DISTRIBUTION VERSION IS
         CALLED IT WRITES A MESSAGE ABOUT THE QUANTITIES TO THE STANDARD 
         OUTPUT CHANNEL THEN TERMINATES THE RUN.  THE VALUES PROVIDED IN 
         THE DISTRIBUTION VERSION SHOULD BE ALTERED, IF NECESSARY, THEN
         THE "WRITE" AND "STOP" STATEMENTS MUST BE COMMENTED OUT.
         ================================================================
         ================================================================

The distribution version of rksuite.for is in DOUBLE PRECISION. A REAL 
version is also available.  When solving ordinary differential equations 
on many popular computers, it is advisable to use DOUBLE PRECISION. However,
if DOUBLE PRECISION provides more than about 20 significant figures, the 
REAL version will usually be satisfactory, provided that the accuracy 
required of the solution is meaningful in the REAL machine precision. 

The values of the coefficients of the Runge-Kutta methods are supplied in
the subroutine CONST as DOUBLE PRECISION constants. Some of the coefficients
are given to 30 significant figures.  It is likely that your compiler will 
round these constants correctly to the nearest representable machine number.
If possible, you should request your compiler to round the constants rather
than truncate them.  Your compiler might warn you that it has shortened
the representation of the constants. The warning does not imply anything is
wrong, but you might wish to take action to avoid receiving such messages
every time you compile RKSUITE. 


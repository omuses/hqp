#define _TAYUTILC_CPP_
#define _ADOLC_SRC_
/*
   --------------------------------------------------------------
   File tayutilC.C of ADOL-C version 1.8.0        as of Nov/30/98
   --------------------------------------------------------------
   The provided class clean_up makes sure the once the
   program leaves, any temporary taylor file is deleted.
   Last changes: 
     981130 olvo: newly created from utils.C

   --------------------------------------------------------------
*/


/****************************************************************************/
/*                                                                 INCLUDES */
#include "dvlparms.h"
#include "usrparms.h"
#include "tayutil.h"

#include <stdio.h>


/****************************************************************************/
/*                                                            CLASS CLEANUP */
/* Added class clean-up, so that when the program leaves, it will clean
   up the temporary file */

/*--------------------------------------------------------------------------*/
/* class definition ( constructor & destructor only) */
class cleanup {
  int valid;
public:
  cleanup();
  ~cleanup();
};

cleanup::cleanup()
{ valid = 0;
}

cleanup::~cleanup()
{ if (taylor_access())
  { close_taylor();
    remove(FNAME3); /*     Complies with ANSI standard */ 
    /*   unlink(FNAME3);   works on some UNIX systems */ 
  }
}


/*--------------------------------------------------------------------------*/
/* one static instance that does all work */
static cleanup at_end;


/****************************************************************************/
/*                                                               THAT'S ALL */

#undef _ADOLC_SRC_
#undef _TAYUTILC_CPP_


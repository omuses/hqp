/*---------------------------------------------------------------------------- 
 ADOL-C -- Automatic Differentiation by Overloading in C++
 File:     tayutilc.cpp
 Revision: $Id: tayutilc.cpp,v 1.1 2004/10/13 14:18:12 e_arnold Exp $
 Contents: The provided class clean_up makes sure the once the
           program leaves, any temporary taylor file is deleted.

 Copyright (c) 2004
               Technical University Dresden
               Department of Mathematics
               Institute of Scientific Computing
  
 This file is part of ADOL-C. This software is provided under the terms of
 the Common Public License. Any use, reproduction, or distribution of the
 software constitutes recipient's acceptance of the terms of this license.
 See the accompanying copy of the Common Public License for more details.

 History:
          20040423 kowarz: adapted to configure - make - make install
          19981130 olvo:   newly created from utils.C

----------------------------------------------------------------------------*/

/****************************************************************************/
/*                                                                 INCLUDES */
#include "tayutil.h"
#include "tayutil_p.h"

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

inline cleanup::~cleanup()
{ if (taylor_access())
  { 
    close_taylor();
    remove(FNAME3); /*     Complies with ANSI standard */ 
    /*   unlink(FNAME3);   works on some UNIX systems */ 
  }
}


/*--------------------------------------------------------------------------*/
/* one static instance that does all work */
static cleanup at_end;


/****************************************************************************/
/*                                                               THAT'S ALL */

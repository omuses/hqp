##############################################################################
## Makefile.am -- Process this file with automake to produce Makefile.in
## Revision: $Id: Makefile.am,v 1.1 2004/10/13 15:05:26 e_arnold Exp $
##
## Copyright (C) 2003
##               Technical University Dresden
##               Department of Mathematics
##               Institute of Scientific Computing
##############################################################################

AUTOMAKE_OPTIONS     = gnu

EXTRA_DIST           = TODO BUGS

if DOCEXA
SUBDIRS              = adolc examples doc
else
if ADDEXA
SUBDIRS              = adolc examples doc
else
SUBDIRS              = adolc doc
endif
endif

MAINTAINERCLEANFILES = Makefile.in aclocal.m4 configure adolc/config-h.in \
                       adolc/stamp-h.in *~

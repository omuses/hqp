# Generated automatically from makefile.in by configure.
#
# Makefile for Meschach via autoconf
#
# Copyright (C) David Stewart & Zbigniew Leyk 1993
#
# $Id: makefile,v 1.4 2002/12/09 10:57:47 e_arnold Exp $
#

include ../makedefs

#OBJ_SUFFIX = .o
#MES_CC = gcc

#MES_CFLAGS = -O2 -fPIC -Wall -pedantic
#MES_DEFS = -DU_INT_DEF=1

#MES_CFLAGS = -O2 -fPIC -Wall -ansi -pedantic
#MES_DEFS = 

O = $(OBJ_SUFFIX)

srcdir = .
VPATH = .

CC = $(MES_CC)

DEFS = $(MES_DEFS)
LIBS =  -lm
RANLIB = ranlib


CFLAGS = $(MES_CFLAGS)


.c$O:
	$(CC) -c $(CFLAGS) $(DEFS) $<

SHELL = /bin/sh
MES_PAK = mesch12b
TAR = tar
SHAR = stree -u
ZIP = zip -r -l
FLIST = FILELIST

###############################

LIST1 = copy$O err$O matrixio$O memory$O vecop$O matop$O pxop$O \
	submat$O init$O otherio$O machine$O matlab$O ivecop$O version$O \
	meminfo$O memstat$O
LIST2 = lufactor$O bkpfacto$O chfactor$O qrfactor$O solve$O hsehldr$O \
	givens$O update$O norm$O hessen$O symmeig$O schur$O svd$O fft$O \
	mfunc$O bdfactor$O
LIST3 = sparse$O sprow$O sparseio$O \
	iter0$O splufctr$O \
	itersym$O iternsym$O \
	spchfctr$O spbkp$O spswap$O 
ZLIST1 = zmachine$O zcopy$O zmatio$O zmemory$O zvecop$O zmatop$O znorm$O \
	 zfunc$O 
ZLIST2 = zlufctr$O zsolve$O zmatlab$O zhsehldr$O zqrfctr$O \
         zgivens$O  zhessen$O zschur$O
LISTHQP = $(LIST1) $(LIST2) sparse$O sprow$O sparseio$O iter0$O splufctr$O \
	addon_hqp$O spaddon_hqp$O addon2_hqp$O

# they are no longer supported
# if you use them add oldpart to all and sparse
OLDLIST = conjgrad$O lanczos$O arnoldi$O

ALL_LISTS = $(LIST1) $(LIST2) $(LIST3) $(ZLIST1) $(ZLIST2) $(OLDLIST)

HBASE = err.h meminfo.h machine.h matrix.h

HLIST = $(HBASE) iter.h matlab.h matrix2.h  oldnames.h sparse.h \
	sparse2.h  zmatrix.h zmatrix2.h

TORTURE = torture$O sptort$O ztorture$O memtort$O itertort$O \
	 mfuntort$O iotort$O

OTHERS = dmacheps.c extras.c fmacheps.c maxint.c  makefile.in \
	 README configure configure.in machine.h.in copyright \
	 tutorial.c tutadv.c rk4.dat ls.dat makefile $(FLIST)


# Different configurations
# default: Hqp configuration
hqp:	$(LISTHQP)
$(LISTHQP): $(HBASE) matrix2.h sparse.h sparse2.h 

# the dependencies **between** the parts are for dmake
all:  part1 part2 part3 zpart1 zpart2 
part2: part1
part3: part2
basic: part1 part2
sparse: part1 part2 part3 
zpart2: zpart1
complex: part1 part2 zpart1 zpart2


$(LIST1): $(HBASE)
part1: $(LIST1)
	ar ru meschach.a $(LIST1); $(RANLIB) meschach.a
	$(RANLIB) meschach.a

$(LIST2): $(HBASE) matrix2.h
part2: $(LIST2)
	ar ru meschach.a $(LIST2); $(RANLIB) meschach.a
	$(RANLIB) meschach.a

$(LIST3): $(HBASE) sparse.h sparse2.h
part3: $(LIST3)
	ar ru meschach.a $(LIST3); $(RANLIB) meschach.a
	$(RANLIB) meschach.a

$(ZLIST1): $(HBASE) zmatrix.h
zpart1: $(ZLIST1)
	ar ru meschach.a $(ZLIST1); $(RANLIB) meschach.a
	$(RANLIB) meschach.a

$(ZLIST2): $(HBASE) zmatrix.h zmatrix2.h 
zpart2: $(ZLIST2)
	ar ru meschach.a $(ZLIST2); $(RANLIB) meschach.a
	$(RANLIB) meschach.a

$(OLDLIST): $(HBASE) sparse.h sparse2.h 
oldpart: $(OLDLIST)
	ar ru meschach.a $(OLDLIST); $(RANLIB) meschach.a
	$(RANLIB) meschach.a

#######################################

tar:
	- /bin/rm -f $(MES_PAK).tar
	chmod 644 `echo $(ALL_LISTS) | sed -e 's/\$O/.c/g'` \
	$(OTHERS) $(HLIST)  `echo $(TORTURE) | sed -e 's/\$O/.c/g'` 
	chmod 755 configure
	$(MAKE) list
	$(TAR) cvf $(MES_PAK).tar \
	 `echo $(ALL_LISTS) | sed -e 's/\$O/.c/g'` \
	$(HLIST)  $(OTHERS) \
	`echo $(TORTURE) | sed -e 's/\$O/.c/g'` \
	MACHINES DOC

# use this only for PC machines	
msdos-zip:
	- /bin/rm -f $(MES_PAK).zip
	chmod 644 `echo $(ALL_LISTS) | sed -e 's/\$O/.c/g'` \
	$(OTHERS) $(HLIST) `echo $(TORTURE) | sed -e 's/\$O/.c/g'` 
	chmod 755 configure
	$(MAKE) list
	$(ZIP)  $(MES_PAK).zip \
	 `echo $(ALL_LISTS) | sed -e 's/\$O/.c/g'` \
	$(HLIST)  $(OTHERS) `echo $(TORTURE) | sed -e 's/\$O/.c/g'` \
	MACHINES DOC


fullshar:
	- /bin/rm -f $(MES_PAK).shar;
	chmod 644 `echo $(ALL_LISTS) | sed -e 's/\$O/.c/g'` \
	$(OTHERS) $(HLIST) `echo $(TORTURE) | sed -e 's/\$O/.c/g'` 
	chmod 755 configure
	$(MAKE) list
	$(SHAR) `echo $(ALL_LISTS) | sed -e 's/\$O/.c/g'` \
	$(HLIST)  $(OTHERS) `echo $(TORTURE) | sed -e 's/\$O/.c/g'` \
	MACHINES DOC > $(MES_PAK).shar

shar:
	- /bin/rm -f meschach1.shar meschach2.shar meschach3.shar \
	meschach4.shar oldmeschach.shar meschach0.shar 
	chmod 644 `echo $(ALL_LISTS) | sed -e 's/\$O/.c/g'` \
	$(OTHERS) $(HLIST) `echo $(TORTURE) | sed -e 's/\$O/.c/g'` 
	chmod 755 configure	
	$(MAKE) list
	$(SHAR) `echo $(LIST1) | sed -e 's/\$O/.c/g'` > meschach1.shar
	$(SHAR) `echo $(LIST2) | sed -e 's/\$O/.c/g'` > meschach2.shar
	$(SHAR) `echo $(LIST3) | sed -e 's/\$O/.c/g'` > meschach3.shar	
	$(SHAR) `echo $(ZLIST1) | sed -e 's/\$O/.c/g'` \
	  `echo $(ZLIST2) | sed -e 's/\$O/.c/g'` > meschach4.shar
	$(SHAR) `echo $(OLDLIST) | sed -e 's/\$O/.c/g'` > oldmeschach.shar
	$(SHAR) $(OTHERS) `echo $(TORTURE) | sed -e 's/\$O/.c/g'` \
	  $(HLIST)  DOC MACHINES  > meschach0.shar

list:
	/bin/rm -f $(FLIST)
	ls -lR `echo $(ALL_LISTS) | sed -e 's/\$O/.c/g'` \
	`echo $(TORTURE) | sed -e 's/\$O/.c/g'` \
	$(HLIST) $(OTHERS) MACHINES DOC \
	|awk '/^$$/ {print};/^[-d]/ {printf("%s %s   %10d %s %s %s %s\n", \
	 $$1,$$2,$$5,$$6,$$7,$$8,$$9)}; /^[^-d]/ {print}' \
	> $(FLIST)



clean:
	/bin/rm -f *$O core asx5213a.mat iotort.dat 

cleanup:
	/bin/rm -f *$O core asx5213a.mat iotort.dat *.a 

realclean:
	/bin/rm -f *$O core asx5213a.mat iotort.dat *.a
	/bin/rm -f torture sptort ztorture memtort itertort mfuntort iotort
#	/bin/rm -f makefile machine.h config.status maxint macheps

alltorture: torture sptort ztorture memtort itertort mfuntort iotort

torture:torture$O meschach.a
	$(CC) $(CFLAGS) $(DEFS) -o torture torture$O \
	meschach.a $(LIBS)
sptort:sptort$O meschach.a
	$(CC) $(CFLAGS) $(DEFS) -o sptort sptort$O \
	meschach.a $(LIBS)
memtort: memtort$O meschach.a
	$(CC) $(CFLAGS) $(DEFS) -o memtort memtort$O \
	meschach.a $(LIBS)
ztorture:ztorture$O meschach.a
	$(CC) $(CFLAGS) $(DEFS) -o ztorture ztorture$O \
	meschach.a $(LIBS)
itertort: itertort$O meschach.a
	$(CC) $(CFLAGS) $(DEFS) -o itertort itertort$O \
	meschach.a $(LIBS)

iotort: iotort$O meschach.a
	$(CC) $(CFLAGS) $(DEFS) -o iotort iotort$O \
	meschach.a $(LIBS)
mfuntort: mfuntort$O meschach.a
	$(CC) $(CFLAGS) $(DEFS) -o mfuntort mfuntort$O \
	meschach.a $(LIBS)
tstmove: tstmove$O meschach.a
	$(CC) $(CFLAGS) $(DEFS) -o tstmove tstmove$O \
	meschach.a $(LIBS)
tstpxvec: tstpxvec$O meschach.a
	$(CC) $(CFLAGS) $(DEFS) -o tstpxvec tstpxvec$O \
	meschach.a $(LIBS)


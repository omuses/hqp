#**************************************************************************
#
#  --------------------------------------------------------------
#  makefile (INS) of ADOL-C version 1.8.7       as of Mar/10/2000
#  --------------------------------------------------------------
#  makefile for generating all other makefiles and total clean up
#
#**************************************************************************

COMPONENTS = makefile.SRC.head makefile.SRC.tail


gnuinstall: $(COMPONENTS) gnu_comp
	cat makefile.SRC.head gnu_comp makefile.SRC.tail > ../SRC/makefile
	@echo installation of makefiles finished


aixinstall: $(COMPONENTS) aix_comp
	cat makefile.SRC.head aix_comp makefile.SRC.tail > ../SRC/makefile
	@echo installation of makefiles finished


suninstall: $(COMPONENTS) sun_comp
	cat makefile.SRC.head sun_comp makefile.SRC.tail > ../SRC/makefile
	@echo installation of makefiles finished


decinstall: $(COMPONENTS) dec_comp
	cat makefile.SRC.head dec_comp makefile.SRC.tail > ../SRC/makefile
	@echo installation of makefiles finished


sgiinstall: $(COMPONENTS) sgi_comp
	cat makefile.SRC.head sgi_comp makefile.SRC.tail > ../SRC/makefile
	@echo installation of makefiles finished


xxxinstall: $(COMPONENTS) xxx_comp
	cat makefile.SRC.head xxx_comp makefile.SRC.tail > ../SRC/makefile
	@echo installation of makefiles finished


uninstall:
	( cd ../SRC; make cleanall)
	-rm ../SRC/makefile
	-rm *~
	@echo uninstallation finished










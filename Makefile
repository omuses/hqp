include makedefs

all:
	$(MAKE) -f Makefile.hqp
	$(MAKE) omuses

ASRCS = any
omuses:
	cd omu; $(MAKE); cd ..
	mv omu/$(LIB_PREFIX)omu.* lib/
	if test -f lib/omu.dll; then cp lib/omu.dll odc/; fi
	cd odc; $(MAKE); cd ..
	$(MAKE) test

test:
	@if test -n "$(findstring OMU_WITH_ADOLC,$(OMU_DEFS))"; then \
	  echo $(findstr OMU_WITH_ADOLC,$(OMU_DEFS)) \
	  echo "Testing Crane example..."; \
	  cd odc; ./run Crane; cd ..; \
	fi
	@if test -n "yes"; then \
	  echo "Testing external S-function..."; \
	  cd odc; ./run dic_mex_sfunction_est; cd ..; \
	fi

doc::
	cd doc; doxygen; cd ..
#	cd doc/latex; $(MAKE) refman.pdf; mv refman.pdf ..; cd ../..

clean:
	cd odc; $(MAKE) clean; cd ..
	rm -f hxi/*~
	cd malloc; $(MAKE) clean; cd ..
	cd omu; $(MAKE) clean; cd ..
	$(MAKE) -f Makefile.hqp clean
	cd hqp_docp; $(MAKE) clean; cd ..
	rm -f hqp_cute/*~
	rm -f doc/*~

distclean: clean
	rm -f makedefs makedirs odc/Makefile odc/mex.tcl hqp_docp/Makefile
	rm -f lib/pkgIndex.tcl
	rm -rf doc/Doxyfile doc/latex doc/refman.pdf
	rm -rf *.cache

LIB_DIR_ROOT = $(INSTALL_PREFIX)/lib
INC_DIR_ROOT = $(INSTALL_PREFIX)/include
INC_DIR = $(INC_DIR_ROOT)/hqp
install::
	$(MAKE) -f Makefile.hqp install
	@PWD=`pwd`
# install libomu
	@if test ! -d $(LIB_DIR_ROOT); then mkdir $(LIB_DIR_ROOT); fi
	$(INSTALL) lib/$(LIB_PREFIX)omu$(LIB_SUFFIX) \
	  $(LIB_DIR_ROOT)/$(LIB_PREFIX)omu-$(VERSION)$(LIB_SUFFIX)
	rm -f $(LIB_DIR_ROOT)/$(LIB_PREFIX)omu$(LIB_SUFFIX)
	cd "$(LIB_DIR_ROOT)"; \
	ln -s $(LIB_PREFIX)omu-$(VERSION)$(LIB_SUFFIX) \
	  $(LIB_PREFIX)omu$(LIB_SUFFIX); \
	cd "$(PWD)"
# make directory structure for includes
	@if test ! -d $(INC_DIR_ROOT); then mkdir $(INC_DIR_ROOT); fi
	@if test ! -d $(INC_DIR)-$(VERSION); then \
	  mkdir $(INC_DIR)-$(VERSION); fi
# install include files
	for f in omu/*.h hxi/*.h; do \
	  $(INSTALL_DATA) $$f $(INC_DIR)-$(VERSION)/; done
# complete directory structure for includes
	rm -rf $(INC_DIR_ROOT)/hqp
	cd "$(INC_DIR_ROOT)"; \
	ln -s hqp-$(VERSION) hqp; \
	cd "$(PWD)"

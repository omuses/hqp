include makedefs

all:
	$(MAKE) -f Makefile.hqp
	if test ! "$(CXX)" = "cl -nologo"; then $(MAKE) omuses; fi

omuses:
	cd adol-c/INS; $(MAKE) xxxinstall; cd ../..
	cd adol-c/SRC; $(MAKE); cd ../..
	cd omu; $(MAKE); cd ..
	mv omu/$(LIB_PREFIX)omu$(LIB_SUFFIX) lib/
	cd odc; $(MAKE); ./run Crane; cd ..

clean:
	cd odc; $(MAKE) clean; cd ..
	cd malloc; $(MAKE) clean; cd ..
	cd omu; $(MAKE) clean; cd ..
	cd adol-c/SRC; $(MAKE) cleanall; rm -f makefile; cd ../..
	$(MAKE) -f Makefile.hqp clean
	cd hqp_docp; $(MAKE) clean; cd ..
	rm -f hqp_cute/*~

distclean: clean
	rm -f makedefs makedirs odc/Makefile hqp_docp/Makefile

LIB_DIR = $(INSTALL_PREFIX)/lib
INC_DIR_ROOT = $(INSTALL_PREFIX)/include
INC_DIR = $(INC_DIR_ROOT)/hqp
install:
	@if test ! -d $(LIB_DIR); then mkdirhier $(LIB_DIR); fi
	$(INSTALL) lib/libhqp$(LIB_SUFFIX) lib/libomu$(LIB_SUFFIX) $(LIB_DIR)
	@if test ! -d $(INC_DIR_ROOT); then mkdir $(INC_DIR_ROOT); fi
	@if test ! -d $(INC_DIR); then mkdir $(INC_DIR); fi
	for f in meschach/*.h iftcl/*.h hqp/*.h adol-c/SRC/*.h omu/*.h; do \
	  $(INSTALL_DATA) $$f $(INC_DIR); done
	@if test ! -d $(INC_DIR)/DRIVERS; then mkdir $(INC_DIR)/DRIVERS; fi
	for f in adol-c/SRC/DRIVERS/*.h; do \
	  $(INSTALL_DATA) $$f $(INC_DIR)/DRIVERS; done
	@if test ! -d $(INC_DIR)/SPARSE; then mkdir $(INC_DIR)/SPARSE; fi
	for f in adol-c/SRC/SPARSE/*.h; do \
	  $(INSTALL_DATA) $$f $(INC_DIR)/SPARSE; done
	@if test ! -d $(INC_DIR)/TAPEDOC; then mkdir $(INC_DIR)/TAPEDOC; fi
	for f in adol-c/SRC/TAPEDOC/*.h; do \
	  $(INSTALL_DATA) $$f $(INC_DIR)/TAPEDOC; done

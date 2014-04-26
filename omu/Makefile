# Makefile for omu
# ("makedepend" may be used for depencies)

include ../makedirs
include ../makedefs

O = $(OBJ_SUFFIX)

CCDEFS = -I. $(MES_INCDIR) $(TCL_INCDIR) $(IF_INCDIR) \
          $(HQP_INCDIR) $(ADOL_INCDIR) $(HXI_INCDIR) $(MEX_INCDIR) \
          $(OMU_DEFS)

ASRCS =	Omu_IntOdeTs.C

SRCS =	Omu_Init.C \
	Hqp_Omuses.C \
	Omu_Variables.C \
	Omu_Vars.C \
	Omu_Dependents.C \
	Omu_Deps.C \
	Omu_Program.C \
	Omu_Integrator.C \
	Omu_IntODE.C \
	Omu_IntEuler.C \
	Omu_IntRK4.C \
	Omu_IntDopri5.C \
	Omu_IntIMP.C \
	Omu_IntGRK4.C \
	Omu_IntSDIRK.C \
        $(OMU_INTRKSUITE_C) \
        $(OMU_INTDASPK_C) \
	adoublev.C \
	Omu_Model.C \
	Prg_SFunctionOpt.C \
	Prg_SFunctionEst.C \
	Hxi_SFunction.C \
	Hxi_SimStruct.C \
	$(MEX_SRCS)

OBJS = $(SRCS:.C=$O) $(RKSUITE_O) fmi$O sfun_fmu$O

LIBRARY = omu

lib: $(OBJS) $(DASPK) Makefile
	rm -f $(LIB_PREFIX)$(LIBRARY)$(MKLIB_SUFFIX)
	$(LD)$(LIB_PREFIX)$(LIBRARY)$(MKLIB_SUFFIX) $(OBJS) \
	$(DASPK_OBJS) $(OMU_MACH_OBJS)
	$(RANLIB) $(LIB_PREFIX)$(LIBRARY)$(LIB_SUFFIX)

fmi.c: ../hxi/fmi.tcl
	../hqp/tpc -ts 80 -o fmi.c ../hxi/fmi.tcl

daspk:
	PWD=`pwd`; cd ../daspk; $(MAKE); cd "$(PWD)"

$(RKSUITE_O): $(RKSUITE_O:$O=.f)
	PWD=`pwd`; cd "`dirname $(RKSUITE_O)`"; \
	$(FORTRAN_COMP) $(FORTRAN_FLAGS) -c rksuite.f; cd "$(PWD)"

sfun_fmu$O: ../hxi/sfun_fmu.c
	$(CC) -c $(CFLAGS) $(CCDEFS) ../hxi/sfun_fmu.c

Hxi_MEX_SFunction$O: ../hxi/Hxi_MEX_SFunction.C
	$(CXX) -c $(CXXFLAGS) $(CCDEFS) ../hxi/Hxi_MEX_SFunction.C

Hxi_SFunction$O: ../hxi/Hxi_SFunction.C
	$(CXX) -c $(CXXFLAGS) $(CCDEFS) ../hxi/Hxi_SFunction.C

Hxi_SimStruct$O: ../hxi/Hxi_SimStruct.C
	$(CXX) -c $(CXXFLAGS) $(CCDEFS) ../hxi/Hxi_SimStruct.C

.C$O:
	$(CXX) -c $(CXXFLAGS) $(CCDEFS) $<

.c$O:
	$(CC) -c $(CFLAGS) $(CCDEFS) $<

clean:
	rm -f fmi.c *.o *.obj *.[pi]db $(LIB_PREFIX)omu.* *~
	PWD=`pwd`; cd ../daspk; $(MAKE) clean; cd "$(PWD)"
	PWD=`pwd`; cd ../rksuite; rm -f *.o *.obj *~; cd "$(PWD)"

depend:
	makedepend -- $(CCFLAGS) -- $(SRCS)

# DO NOT DELETE THIS LINE -- make depend depends on it.

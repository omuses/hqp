# Makefile for libhqp
# ("makedepend" is used for depencies)

include ../makedirs
include ../makedefs

O = $(OBJ_SUFFIX)

CCDEFS = $(MES_INCDIR) $(TCL_INCDIR) $(IF_INCDIR) $(ASCEND_INCDIR) \
	 $(LPSOLVE_INCDIR) $(HQP_INCDIR) $(HQP_DEFS)

SRCS = 	Hqp_Init.C \
	sprcm.C \
        spBKP.C \
        matBKP.C \
        bdBKP.C \
	Hqp_impl.C \
	Hqp_DynLoad.C \
	Hqp_Program.C \
	Hqp_Solver.C \
	Hqp_Client.C \
	Hqp_IpsFranke.C \
	Hqp_IpsMehrotra.C \
	Hqp_IpMatrix.C \
	Hqp_IpSpBKP.C \
	Hqp_IpRedSpBKP.C \
	Hqp_IpLQDOCP.C \
	Hqp_IpPardiso.C \
	t_mesch.C \
	Hqp_IpSpSC.C \
	meschext_hl.C \
	Hqp_SqpSolver.C \
	Hqp_SqpPowell.C \
	Hqp_SqpSchittkowski.C \
	Hqp_HL.C \
	Hqp_HL_Gerschgorin.C \
	Hqp_HL_DScale.C \
	Hqp_HL_BFGS.C \
	Hqp_HL_SparseBFGS.C \
	Hqp_SqpProgram.C \
	Hqp_Docp.C \
        Hqp_MipSolver.C \
        $(LPSOLVE_SRCS) \
	$(CUTE_SRCS) \
	$(ASCEND_SRCS)

OBJS = $(SRCS:.C=$O) hqp_solve$O
# add own qsort to HQP sources for comparison on different platforms
#OBJS = $(SRCS:.C=$O) qsort$O hqp_solve$O

LIBRARY = hqp

all: lib

lib: $(OBJS) Makefile
	rm -f $(LIB_PREFIX)$(LIBRARY)$(MKLIB_SUFFIX)
	$(LD)$(LIB_PREFIX)$(LIBRARY)$(MKLIB_SUFFIX) $(OBJS) \
          ../meschach/*$O ../iftcl/*$O $(HQP_MACH_OBJS)
	$(RANLIB) $(LIB_PREFIX)$(LIBRARY)$(LIB_SUFFIX)

hqp_solve.c: tpc hqp_solve.tcl
	./tpc -ts 80 -o hqp_solve.c hqp_solve.tcl

tpc:
	$(CC) -o tpc tpc.c

.C$O:
	$(CXX) -c $(CXXFLAGS) $(CCDEFS) $<

.c$O:
	$(CC) -c $(CFLAGS) $(CCDEFS) $<

depend:
	makedepend -- $(CCDEFS) -- $(SRCS)

clean:
	rm -f tpc$(EXE_SUFFIX) hqp_solve.c *.o *.obj *.[pi]db \
	  $(LIB_PREFIX)hqp.* *~

# DO NOT DELETE THIS LINE -- make depend depends on it.

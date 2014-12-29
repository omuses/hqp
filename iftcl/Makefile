# Makefile for iftcl
# rf, 6/20/94
# ("makedepend" is used for depencies)

include ../makedefs
include ../makedirs

O = $(OBJ_SUFFIX)

CCDEFS = -I. $(HQP_INCDIR) $(MES_INCDIR) $(TCL_INCDIR) $(IF_DEFS)

SRCS = 	If.C \
	If_List.C \
	If_ListElement.C \
	If_Element.C \
	If_Procedure.C \
	If_Int.C \
	If_IntVec.C \
	If_Bool.C \
	If_Real.C \
	If_RealVec.C \
	If_RealMat.C \
	If_String.C \
	If_Class.C

OBJS = $(SRCS:.C=$O)

all: $(OBJS) Makefile

.C$O:
	$(CXX) -c $(CXXFLAGS) $(CCDEFS) $<

iftest: $(OBJS) iftest$O Makefile
	$(CXX) -o iftest $(OBJS) iftest$O ../meschach/*$O ../hqp/sprcm$O \
	$(LDFLAGS_START) $(TCL_LIBDIR) $(TCL_STUB_LIB) $(TCL_LIBS) -lm

depend:
	makedepend -- $(CCFLAGS) -- $(SRCS)

clean:
	rm -f *.o *.obj *.[pi]db *.ilk *.lib *.exp *.exe *~

# DO NOT DELETE THIS LINE -- make depend depends on it.

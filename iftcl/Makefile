# Makefile for iftcl
# rf, 6/20/94
# ("makedepend" is used for depencies)

include ../makedefs
include ../makedirs

O = $(OBJ_SUFFIX)

CCDEFS = $(HQP_INCDIR) $(MES_INCDIR) $(TCL_INCDIR) $(IF_DEFS)

SRCS = 	If.C \
	If_List.C \
	If_ListElement.C \
	If_Element.C \
	If_Command.C \
	If_Proc.C \
	If_Variable.C \
	If_Int.C \
	If_Bool.C \
	If_Real.C \
	If_RealVec.C \
	If_RealMat.C \
	If_IntVec.C \
	If_Class.C

OBJS = $(SRCS:.C=$O)

all: $(OBJS) Makefile

.C$O:
	$(CXX) -c $(CXXFLAGS) $(CCDEFS) $<

iftest: $(OBJS) iftest$O Makefile
	$(CC) -o iftest $(CCFLAGS) $(OBJS) iftest$O \
        ../meschach/*$O ../hqp/sprcm$O -ltcl8.0 -lm -ldl

depend:
	makedepend -- $(CCFLAGS) -- $(SRCS)

clean:
	rm -f *.o *.obj *~

# DO NOT DELETE THIS LINE -- make depend depends on it.

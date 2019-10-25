export PETSC_ARCH:=arch-darwin-c-real
include ${PETSC_DIR}/lib/petsc/conf/variables
#include ${PETSC_DIR}/lib/petsc/conf/rules
include ${SLEPC_DIR}/lib/slepc/conf/slepc_common


# the following provides CCPPFLAGS, CLINKER, CC_INCLUDES, SLEPC_INCLUDE, SLEPC_LIB
#include $(SLEPC_DIR)/conf/slepc_common

INCPATHS = -I. -I../common -I../atmosphere -I/usr/local/include -I/usr/include ${PETSC_CC_INCLUDES} ${SLEPC_CC_INCLUDES}

all: Modess

.PHONY: clean

# link	
Modess: Modess_main.o SolveModNB.o ProcessOptionsNB.o Modess_lib.o
	${CXX_LINKER} -o $@ $^  -L/usr/lib ../../lib/libatmosphere.a ../../lib/libcommon.a  ${CXX_LINKER_FLAGS} ${SLEPC_LIB} ${PETSC_LIB} -lgsl -lgslcblas -lm -lfftw3 
	cp $@ ../../bin
	
# compile 
%.o: %.cpp
	${CXX} ${INCPATHS} -fpic -c -Wall ${CXX_FLAGS} -Wwrite-strings -Wno-strict-aliasing -Wno-unknown-pragmas -o $@ $<

clean::
	-$(RM) -rf Modess_main.o SolveModNB.o ProcessOptionsNB.o Modess_lib.o Modess

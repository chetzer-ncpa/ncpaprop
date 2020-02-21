#SLEPC_DIR=/home/claus/dev/ncpaprop/src/extern/slepc-3.2-p3
#export PETSC_VERSION:=3.2-p5
#export SLEPC_VERSION:=3.2-p3
#export PETSC_DIR:=/home/claus/dev/ncpaprop/src/extern/petsc-3.2-p5
#export SLEPC_DIR:=/home/claus/dev/ncpaprop/src/extern/slepc-3.2-p3
#export PETSC_ARCH:=arch-linux-gnu-c-debug
include ${PETSC_DIR}/lib/petsc/conf/variables
include ${SLEPC_DIR}/lib/slepc/conf/slepc_common

# the following provides CCPPFLAGS, CLINKER, CC_INCLUDES, SLEPC_INCLUDE, SLEPC_LIB
#include $(SLEPC_DIR)/conf/slepc_common

#INCPATHS=-I. -I../common -I../atmosphere -I/usr/local/include -I/usr/include ${CCPPFLAGS}
INCPATHS=-I. -I../common -I../atmosphere -I/usr/local/include -I/usr/include ${PETSC_CC_INCLUDES} ${SLEPC_CC_INCLUDES}
OBJS=wnlrt.o ProcessOptionsNRT.o nonlinearRay.o linearRay3DStrat.o
TARGET=wnlrt

#CXX=mpicxx
#CXX_LINKER=mpicxx

all: $(TARGET)

.PHONY: clean

# link	
$(TARGET): $(OBJS)
	${CXX_LINKER} -o $@ $^  -L/usr/lib ../../lib/libatmosphere.a ../../lib/libcommon.a -lgsl -lgslcblas -lm -lfftw3 
	cp $@ ../../bin
	
# compile 
%.o: %.cpp
	${CXX} -fpic -c -Wall ${INCPATHS} -o $@ $<

clean::
	-$(RM) -rf $(OBJS) $(TARGET)

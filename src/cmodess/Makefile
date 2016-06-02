#SLEPC_DIR=/home/claus/dev/ncpaprop/src/extern/slepc-3.2-p3
export PETSC_VERSION:=3.2-p5
export SLEPC_VERSION:=3.2-p3
export PETSC_DIR:=/home/claus/dev/ncpaprop/src/extern/petsc-3.2-p5
export SLEPC_DIR:=/home/claus/dev/ncpaprop/src/extern/slepc-3.2-p3
export PETSC_ARCH:=arch-linux-gnu-c-complex

# the following provides CCPPFLAGS, CLINKER, CC_INCLUDES, SLEPC_INCLUDE, SLEPC_LIB
include $(SLEPC_DIR)/conf/slepc_common

INCPATHS=-I. -I../common -I../atmosphere -I/usr/local/include -I/usr/include ${CCPPFLAGS}
OBJS=CModess_main.o SolveCModNB.o ProcessOptionsNB.o CModess_lib.o
TARGET=CModess

all: $(TARGET)

.PHONY: clean

# link	
$(TARGET): $(OBJS)
	g++ -o $@ $^  -L/usr/lib ../../lib/libatmosphere.a ../../lib/libcommon.a $(SLEPC_LIB) -lgsl -lgslcblas -lm -lfftw3 
	cp $@ ../../bin
	
# compile 
%.o: %.cpp
	g++ -fpic -c -ggdb -Wall ${INCPATHS} -Wwrite-strings -Wno-strict-aliasing -Wno-unknown-pragmas -o $@ $<

clean::
	-$(RM) -rf $(OBJS) $(TARGET)

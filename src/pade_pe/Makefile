#SLEPC_DIR=/home/claus/dev/ncpaprop/src/extern/slepc-3.2-p3
export PETSC_VERSION:=3.2-p5
export SLEPC_VERSION:=3.2-p3
export PETSC_DIR:=/home/claus/dev/ncpaprop/src/extern/petsc-3.2-p5
export SLEPC_DIR:=/home/claus/dev/ncpaprop/src/extern/slepc-3.2-p3
export PETSC_ARCH:=arch-linux-gnu-c-debug

# the following provides CCPPFLAGS, CLINKER, CC_INCLUDES, SLEPC_INCLUDE, SLEPC_LIB
include $(SLEPC_DIR)/conf/slepc_common

OBJS=binaryreader.o atmlib.o ProcessOptionsPE.o pe_main.o
TARGET=pape
INCPATHS=-I. -I../common -I../atmosphere -I/usr/include ${CCPPFLAGS}

all: $(TARGET)

.PHONY: clean

# link	
$(TARGET): $(OBJS)
	g++ -o $@ $^  -L/usr/lib ../../lib/libatmosphere.a ../../lib/libcommon.a $(SLEPC_LIB) -lgsl -lgslcblas -lm -lfftw3 
	cp $@ ../../bin
	
# compile 
%.o: %.cpp
	g++ -fpic -c -ggdb -Wall -I. -I../common -I../atmosphere -I/usr/local/include -o $@ $<

clean::
	-$(RM) -rf $(OBJS) $(TARGET)

LIBNAMES=atmosphere common
STATICLIBS=$(foreach f,$(LIBNAMES),../../lib/lib$(f).a)

GSL_LIBNAMES=gsl gslcblas
GSL_LIBS=$(foreach d,$(GSL_LIBNAMES),-l$(d))

# Set up the build order
SUBDIRS= common atmosphere

#COMPLEX_VERSION = 1

#ifeq (${COMPLEX_VERSION}, 0)
#PETSC_ARCH=arch-linux2-c-debug
#PETSC_DIR=/home/doru/lib/petsc-3.2-p5
#SLEPC_DIR=/home/doru/lib/slepc-3.2-p3
#endif

# use the following dirs for PETSC with complex numbers
#ifeq (${COMPLEX_VERSION}, 1)
#PETSC_ARCH=arch-linux2-c-debug-complex
#PETSC_DIR=/home/doru/lib/petsc-3.3-p3
#SLEPC_DIR=/home/doru/lib/slepc-3.3-p2
#endif

# the following provides CCPPFLAGS, CLINKER, CC_INCLUDES, SLEPC_INCLUDE, SLEPC_LIB
#include ${PETSC_DIR}/conf/variables
#include ${PETSC_DIR}/conf/rules

include ${SLEPC_DIR}/conf/slepc_common

PETC_SLEPC_INCLUDES = ${CCPPFLAGS}

export CXX=$(PETSC_DIR)/$(PETSC_ARCH_COMPLEX)/bin/mpic++

INCPATHS=$(foreach d,$(SUBDIRS),-I../$(d)) -I/usr/include ${PETC_SLEPC_INCLUDES}

export CXXFLAGS=-fpic -c -g -Wall $(INCPATHS) -Wwrite-strings -Wno-strict-aliasing -Wno-unknown-pragmas
export LDFLAGS=-L/usr/lib

SOURCES=$(wildcard *.cpp) 
OBJECTS=$(SOURCES:.cpp=.o)

.PHONY: clean

# link	
ModessRD2WCM: ModessRDCM_main.o SolveModNBRDCM.o ProcessOptionsNBRDCM.o ModessRDCM_lib.o 
	${CXX} $(LDFLAGS) -o $@ $^  $(STATICLIBS) $(GSL_LIBS) ${SLEPC_LIB}
	cp $@ ../../bin
	
# compile 
%.o: %.cpp
	${CXX} $(CXXFLAGS) -o $@ $< 		
	
# cleanup 
clean::
	@echo SOURCES = ${SOURCES}
	@echo OBJECTS = ${OBJECTS}
	rm -rf ${OBJECTS} ModessRD2WCM

	
	

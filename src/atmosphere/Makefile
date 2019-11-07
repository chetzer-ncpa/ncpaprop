export PETSC_ARCH:=arch-darwin-c-real
include ${PETSC_DIR}/lib/petsc/conf/variables
#include ${PETSC_DIR}/lib/petsc/conf/rules
include ${SLEPC_DIR}/lib/slepc/conf/slepc_common

SOURCES=AtmosphericProfile.cpp AtmosphericSpecification.cpp JetProfile.cpp ProfileGroup.cpp SampledProfile.cpp Slice.cpp Sounding.cpp
OBJECTS=$(SOURCES:.cpp=.o)

# Set default values for compiler and flags
TARGET=libatmosphere.a
LINKCOMMAND=ar curv 

.PHONY: clean tidy atmosphere

all: atmosphere

atmosphere: libatmosphere.a

libatmosphere.a: $(OBJECTS)
	ar curv $@ $^
	cp $@ ../../lib

%.o: %.cpp
	$(CXX) -fpic -c -ggdb -Wall -I. -I../common -I../atmosphere -I/usr/local/include -o $@ $<

clean::
	rm -rf *.o $(TARGET) $(OBJECTS)

tidy:
	rm -rf *.o


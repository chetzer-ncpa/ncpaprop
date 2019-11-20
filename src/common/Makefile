export PETSC_ARCH:=arch-darwin-c-real
include ${PETSC_DIR}/lib/petsc/conf/variables
#include ${PETSC_DIR}/lib/petsc/conf/rules
include ${SLEPC_DIR}/lib/slepc/conf/slepc_common

SOURCES=anyoption.cpp binaryreader.cpp geographic.cpp util.cpp units.cpp
OBJS=$(SOURCES:.cpp=.o)

# Set default values for compiler and flags
TARGET=libcommon.a
LINKCOMMAND=ar curv 
STANDARD=-std=c++14
DEBUG=-ggdb
WARNING_FLAGS=-Wall
CXXFLAGS=-fpic -c $(DEBUG) $(STANDARD) $(WARNING_FLAGS)
INCLUDEPATH=-I. -I../common -I../atmosphere -I/usr/local/include
#CXX=g++

.PHONY: clean tidy common

all: common

common: libcommon.a

libcommon.a: $(OBJS)
	$(LINKCOMMAND) $@ $^
	cp $@ ../../lib

%.o: %.cpp
	${CXX} $(CXXFLAGS) $(INCLUDEPATH) -o $@ $<

clean::
	rm -rf *.o $(TARGET) $(OBJS)

tidy:
	rm -rf *.o


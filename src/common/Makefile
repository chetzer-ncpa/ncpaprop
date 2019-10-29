SOURCES=anyoption.cpp binaryreader.cpp geographic.cpp util.cpp
OBJECTS=$(SOURCES:.cpp=.o)

# Set default values for compiler and flags
TARGET=libcommon.a
LINKCOMMAND=ar curv 
CXX=g++

.PHONY: clean tidy common

all: common

common: libcommon.a

libcommon.a: $(OBJECTS)
	${LINKCOMMAND} $@ $^
	cp $@ ../../lib

%.o: %.cpp
	${CXX} -fpic -c -ggdb -Wall -I. -I../common -I../atmosphere -I/usr/local/include -o $@ $<

clean:
	rm -rf *.o $(TARGET) $(OBJECTS)

tidy:
	rm -rf *.o


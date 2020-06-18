include /home/claus/code/devel/ncpaprop/extern/petsc/lib/petsc/conf/variables
include /home/claus/code/devel/ncpaprop/extern/slepc/slepc-3.12.2/lib/slepc/conf/slepc_common


# the following provides CCPPFLAGS, CLINKER, CC_INCLUDES, SLEPC_INCLUDE, SLEPC_LIB
#include $(SLEPC_DIR)/conf/slepc_common

INCPATHS = -I. -I../common -I../atmosphere -I/usr/local/include -I/usr/include ${PETSC_CC_INCLUDES} ${SLEPC_CC_INCLUDES}
SOURCES=binaryreader.cpp geographic.cpp util.cpp units.cpp parameterset.cpp
OBJS=$(SOURCES:.cpp=.o)
TARGET=libcommon.a


all: $(TARGET)

.PHONY: clean

# link	
$(TARGET): $(OBJS)
	ar curv $@ $^
	cp $@ ../../lib
	
# compile 
%.o: %.cpp
	${CXX} ${INCPATHS} -fpic -c -Wall ${CXX_FLAGS} -Wwrite-strings -Wno-strict-aliasing -Wno-unknown-pragmas -Wno-unused-but-set-variable -Wno-maybe-uninitialized -Wno-unused-result -o $@ $<

clean::
	-$(RM) -rf $(OBJS) $(TARGET)

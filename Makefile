# Variable check
# CXXFLAGS = -fpic -c -Wall
# INCLUDEFLAGS = -I. -I../common -I../atmosphere -I/usr/local/include -I/usr/include
# WARNINGFLAGS = -Wwrite-strings -Wno-strict-aliasing -Wno-unknown-pragmas -Wno-unused-but-set-variable -Wno-maybe-uninitialized -Wno-unused-result
# STATICLIBS = ../../lib/libatmosphere.a ../../lib/libcommon.a
# LIBS = -lgsl -lgslcblas -lm -lfftw3
# LDFLAGS =  -L/usr/lib
# PETSC_DIR = /home/claus/code/devel/ncpaprop/extern/petsc
# SLEPC_DIR = /home/claus/code/devel/ncpaprop/extern/slepc/slepc-3.12.2
# PETSC_ARCH_REAL = arch-linux-c-real-debug
# PETSC_ARCH_COMPLEX = arch-linux-c-complex-debug
# PETSC_INCLUDE_FILE_GENERIC = /home/claus/code/devel/ncpaprop/extern/petsc/lib/petsc/conf/variables
# SLEPC_INCLUDE_FILE_GENERIC = /home/claus/code/devel/ncpaprop/extern/slepc/slepc-3.12.2/lib/slepc/conf/slepc_common

#package-parts=common atmosphere modess modbb modess_rd_1wcm pade_pe wmod cmodess cmodbb wnlrt tdpape
package-parts=common atmosphere modess 


# Build everything
all: $(package-parts)

# build general utility libraries
common:
	$(MAKE) -C src/common all PETSC_ARCH=arch-linux-c-real-debug PETSC_DIR=/home/claus/code/devel/ncpaprop/extern/petsc SLEPC_DIR=/home/claus/code/devel/ncpaprop/extern/slepc/slepc-3.12.2

# build atmospheric utility libraries
atmosphere:
	$(MAKE) -C src/atmosphere all PETSC_ARCH=arch-linux-c-real-debug PETSC_DIR=/home/claus/code/devel/ncpaprop/extern/petsc SLEPC_DIR=/home/claus/code/devel/ncpaprop/extern/slepc/slepc-3.12.2

# ray tracing routines
raytrace:
	$(MAKE) -C src/raytrace all PETSC_ARCH=arch-linux-c-real-debug PETSC_DIR=/home/claus/code/devel/ncpaprop/extern/petsc SLEPC_DIR=/home/claus/code/devel/ncpaprop/extern/slepc/slepc-3.12.2

# Normal modes, single frequency, effective sound speed approximation
modess:
	$(MAKE) -C src/modess all PETSC_ARCH=arch-linux-c-real-debug PETSC_DIR=/home/claus/code/devel/ncpaprop/extern/petsc SLEPC_DIR=/home/claus/code/devel/ncpaprop/extern/slepc/slepc-3.12.2

# Normal modes, broadband
modbb:
	$(MAKE) -C src/modbb all PETSC_ARCH=arch-linux-c-real-debug PETSC_DIR=/home/claus/code/devel/ncpaprop/extern/petsc SLEPC_DIR=/home/claus/code/devel/ncpaprop/extern/slepc/slepc-3.12.2

# Normal modes, range-dependent, one-way coupled modes
modess_rd_1wcm:
	$(MAKE) -C src/modess_rd_1wcm all PETSC_ARCH=arch-linux-c-real-debug PETSC_DIR=/home/claus/code/devel/ncpaprop/extern/petsc SLEPC_DIR=/home/claus/code/devel/ncpaprop/extern/slepc/slepc-3.12.2

# High-angle high-mach parabolic equation
pade_pe:
	$(MAKE) -C src/pade_pe all PETSC_ARCH=arch-linux-c-real-debug PETSC_DIR=/home/claus/code/devel/ncpaprop/extern/petsc SLEPC_DIR=/home/claus/code/devel/ncpaprop/extern/slepc/slepc-3.12.2

# Wide-angle high-mach number normal modes, single frequency
wmod:
	$(MAKE) -C src/wmod all PETSC_ARCH=arch-linux-c-real-debug PETSC_DIR=/home/claus/code/devel/ncpaprop/extern/petsc SLEPC_DIR=/home/claus/code/devel/ncpaprop/extern/slepc/slepc-3.12.2
	
# Complex normal modes, single frequency, effective sound speed
cmodess:
	$(MAKE) -C src/cmodess all PETSC_ARCH=arch-linux-c-complex-debug PETSC_DIR=/home/claus/code/devel/ncpaprop/extern/petsc SLEPC_DIR=/home/claus/code/devel/ncpaprop/extern/slepc/slepc-3.12.2

# Complex normal modes, broadband
cmodbb:
	$(MAKE) -C src/cmodbb all PETSC_ARCH=arch-linux-c-complex-debug PETSC_DIR=/home/claus/code/devel/ncpaprop/extern/petsc SLEPC_DIR=/home/claus/code/devel/ncpaprop/extern/slepc/slepc-3.12.2

wnlrt:
	$(MAKE) -C src/wnlrt all PETSC_ARCH=arch-linux-c-real-debug PETSC_DIR=/home/claus/code/devel/ncpaprop/extern/petsc SLEPC_DIR=/home/claus/code/devel/ncpaprop/extern/slepc/slepc-3.12.2

tdpape:
	$(MAKE) -C src/tdpape all PETSC_ARCH=arch-linux-c-real-debug PETSC_DIR=/home/claus/code/devel/ncpaprop/extern/petsc SLEPC_DIR=/home/claus/code/devel/ncpaprop/extern/slepc/slepc-3.12.2
	
#modess_rd_2wcm:
#	$(MAKE) -C src/modess_rd_2wcm ModessRD2WCM PETSC_ARCH=arch-linux-c-complex-debug

# clean everything up
clean:
	-$(MAKE) -C src/common clean  PETSC_ARCH=arch-linux-c-real-debug PETSC_DIR=/home/claus/code/devel/ncpaprop/extern/petsc SLEPC_DIR=/home/claus/code/devel/ncpaprop/extern/slepc/slepc-3.12.2
	-$(MAKE) -C src/atmosphere clean  PETSC_ARCH=arch-linux-c-real-debug PETSC_DIR=/home/claus/code/devel/ncpaprop/extern/petsc SLEPC_DIR=/home/claus/code/devel/ncpaprop/extern/slepc/slepc-3.12.2
	-$(MAKE) -C src/raytrace clean  PETSC_ARCH=arch-linux-c-real-debug PETSC_DIR=/home/claus/code/devel/ncpaprop/extern/petsc SLEPC_DIR=/home/claus/code/devel/ncpaprop/extern/slepc/slepc-3.12.2
	-$(MAKE) -C src/modess clean  PETSC_ARCH=arch-linux-c-real-debug PETSC_DIR=/home/claus/code/devel/ncpaprop/extern/petsc SLEPC_DIR=/home/claus/code/devel/ncpaprop/extern/slepc/slepc-3.12.2
	-$(MAKE) -C src/modbb clean  PETSC_ARCH=arch-linux-c-real-debug PETSC_DIR=/home/claus/code/devel/ncpaprop/extern/petsc SLEPC_DIR=/home/claus/code/devel/ncpaprop/extern/slepc/slepc-3.12.2
	-$(MAKE) -C src/modess_rd_1wcm clean  PETSC_ARCH=arch-linux-c-real-debug PETSC_DIR=/home/claus/code/devel/ncpaprop/extern/petsc SLEPC_DIR=/home/claus/code/devel/ncpaprop/extern/slepc/slepc-3.12.2
	-$(MAKE) -C src/pade_pe clean  PETSC_ARCH=arch-linux-c-real-debug PETSC_DIR=/home/claus/code/devel/ncpaprop/extern/petsc SLEPC_DIR=/home/claus/code/devel/ncpaprop/extern/slepc/slepc-3.12.2
	-$(MAKE) -C src/wmod clean  PETSC_ARCH=arch-linux-c-real-debug PETSC_DIR=/home/claus/code/devel/ncpaprop/extern/petsc SLEPC_DIR=/home/claus/code/devel/ncpaprop/extern/slepc/slepc-3.12.2
	-$(MAKE) -C src/cmodess clean  PETSC_ARCH=arch-linux-c-complex-debug PETSC_DIR=/home/claus/code/devel/ncpaprop/extern/petsc SLEPC_DIR=/home/claus/code/devel/ncpaprop/extern/slepc/slepc-3.12.2
	-$(MAKE) -C src/cmodbb clean  PETSC_ARCH=arch-linux-c-complex-debug PETSC_DIR=/home/claus/code/devel/ncpaprop/extern/petsc SLEPC_DIR=/home/claus/code/devel/ncpaprop/extern/slepc/slepc-3.12.2
	-$(MAKE) -C src/tdpape clean  PETSC_ARCH=arch-linux-c-real-debug PETSC_DIR=/home/claus/code/devel/ncpaprop/extern/petsc SLEPC_DIR=/home/claus/code/devel/ncpaprop/extern/slepc/slepc-3.12.2
	-$(MAKE) -C src/wnlrt clean  PETSC_ARCH=arch-linux-c-real-debug PETSC_DIR=/home/claus/code/devel/ncpaprop/extern/petsc SLEPC_DIR=/home/claus/code/devel/ncpaprop/extern/slepc/slepc-3.12.2
#	-$(MAKE) -C src/modess_rd_2wcm clean  PETSC_ARCH=arch-linux-c-complex-debug
	-rm bin/CModBB  bin/CModess  bin/ModBB  bin/Modess  bin/ModessRD1WCM  bin/pape  bin/raytrace.2d  bin/raytrace.3d  bin/WMod lib/libatmosphere.a lib/libcommon.a bin/tdpape bin/wnlrt

# Test all executables to make sure their output matches the dev machine output
test:
	$(MAKE) -C test test

# Run the tests that take forever
testlong:
	$(MAKE) -C test testlong

.PHONY:  common atmosphere raytrace modess modbb modess_rd_1wcm pade_pe wmod cmodess cmodbb modess_rd_2wcm clean test

# Variable check
# CXXFLAGS = -fpic -c -ggdb -Wall
# LIBS = -lgsl -lgslcblas -lm -lfftw3 
# LDFLAGS =  -L/usr/lib
# FC = gfortran
# PETSC_VERSION = 3.2-p5
# SLEPC_VERSION = 3.2-p3
# PETSC_DIR = /home/claus/dev/ncpaprop/src/extern/petsc-3.2-p5
# SLEPC_DIR = /home/claus/dev/ncpaprop/src/extern/slepc-3.2-p3
# PETSC_OS = linux-gnu
PETSC_ARCH_REAL = arch-linux-c-real
PETSC_ARCH_COMPLEX = arch-linux-c-complex
# PACKAGE_NAME = NCPA Propagation Modeling Suite
# PACKAGE_TARNAME = ncpaprop
# PACKAGE_VERSION = 1.3.17
distdir=ncpaprop-1.3.17
tarball=ncpaprop-1.3.17.tar.gz
package-parts=common atmosphere raytrace modess modbb modess_rd_1wcm pade_pe wmod cmodess cmodbb wnlrt tdpape


# Build everything
all: $(package-parts)

dist: $(tarball)

# Build the distribution tarball
$(tarball): $(distdir)
	tar cvzhf $(tarball) $(distdir)
	rm -rf $(distdir)

$(distdir): FORCE
	mkdir -p $(distdir)/bin $(distdir)/lib $(distdir)/src $(distdir)/samples $(distdir)/test/results  $(distdir)/test/CModBB $(distdir)/test/ModBB $(distdir)/test/pape $(distdir)/docs
	cp configure Makefile.in config.guess config.sub install-sh changelog.txt $(distdir)/
	cp docs/*.pdf $(distdir)/docs/
	rsync -auv samples/* $(distdir)/samples/
	rsync -auv test/CModess test/col* test/compare_* test/Makefile.in test/Modess test/ModessRD1WCM test/raytrace.?d test/WMod test/run_*.bash $(distdir)/test/
#	rsync -auv test/CModBB test/ModBB test/pape $(distdir)/test/

	mkdir -p $(distdir)/src/extern
	cp src/extern/Makefile.in src/extern/petsc-3.2-p5.tar.gz src/extern/slepc-3.2-p3.tar.gz $(distdir)/src/extern/

	mkdir -p $(distdir)/src/atmosphere
	cp src/atmosphere/*.h src/atmosphere/*.cpp src/atmosphere/Makefile.in $(distdir)/src/atmosphere/

	mkdir -p $(distdir)/src/cmodbb
	cp src/cmodbb/*.h src/cmodbb/*.cpp src/cmodbb/Makefile.in $(distdir)/src/cmodbb/

	mkdir -p $(distdir)/src/cmodess
	cp src/cmodess/*.h src/cmodess/*.cpp src/cmodess/Makefile.in $(distdir)/src/cmodess/

	mkdir -p $(distdir)/src/common
	cp src/common/*.h src/common/*.cpp src/common/Makefile.in $(distdir)/src/common/


	mkdir -p $(distdir)/src/modbb
	cp src/modbb/*.h src/modbb/*.cpp src/modbb/Makefile.in $(distdir)/src/modbb/

	mkdir -p $(distdir)/src/modess
	cp src/modess/*.h src/modess/*.cpp src/modess/Makefile.in $(distdir)/src/modess/

	mkdir -p $(distdir)/src/modess_rd_1wcm
	cp src/modess_rd_1wcm/*.h src/modess_rd_1wcm/*.cpp src/modess_rd_1wcm/Makefile.in $(distdir)/src/modess_rd_1wcm/

	mkdir -p $(distdir)/src/pade_pe
	cp src/pade_pe/*.h src/pade_pe/*.cpp src/pade_pe/Makefile.in $(distdir)/src/pade_pe/

	mkdir -p $(distdir)/src/raytrace
	cp src/raytrace/*.h src/raytrace/*.cpp src/raytrace/Makefile.in $(distdir)/src/raytrace/

	mkdir -p $(distdir)/src/wmod
	cp src/wmod/*.h src/wmod/*.cpp src/wmod/Makefile.in $(distdir)/src/wmod/

	mkdir -p $(distdir)/src/wnlrt
	cp src/wnlrt/*.h src/wnlrt/*.cpp src/wnlrt/Makefile.in $(distdir)/src/wnlrt

	mkdir -p $(distdir)/src/tdpape
	cp src/tdpape/*.h src/tdpape/*.cpp src/tdpape/Makefile.in $(distdir)/src/tdpape


FORCE:
	-rm $(tarball)
	-rm -rf $(distdir)

# build external libraries (PETSc, SLEPc)
#extern:
#	$(MAKE) -C src/extern extern

# build general utility libraries
common:
	$(MAKE) -C src/common common PETSC_ARCH=$(PETSC_ARCH_REAL)

# build atmospheric utility libraries
atmosphere:
	$(MAKE) -C src/atmosphere atmosphere PETSC_ARCH=$(PETSC_ARCH_REAL)

# ray tracing routines
raytrace:
	$(MAKE) -C src/raytrace raytrace PETSC_ARCH=$(PETSC_ARCH_REAL)

# Normal modes, single frequency, effective sound speed approximation
modess:
	$(MAKE) -C src/modess Modess PETSC_ARCH=$(PETSC_ARCH_REAL)

# Normal modes, broadband
modbb:
	$(MAKE) -C src/modbb ModBB PETSC_ARCH=$(PETSC_ARCH_REAL)

# Normal modes, range-dependent, one-way coupled modes
modess_rd_1wcm:
	$(MAKE) -C src/modess_rd_1wcm ModessRD1WCM PETSC_ARCH=$(PETSC_ARCH_REAL)

# High-angle high-mach parabolic equation
pade_pe:
	$(MAKE) -C src/pade_pe pape PETSC_ARCH=$(PETSC_ARCH_REAL)

# Wide-angle high-mach number normal modes, single frequency
wmod:
	$(MAKE) -C src/wmod WMod PETSC_ARCH=$(PETSC_ARCH_REAL)
	
# Complex normal modes, single frequency, effective sound speed
cmodess:
	$(MAKE) -C src/cmodess CModess PETSC_ARCH=$(PETSC_ARCH_COMPLEX)

# Complex normal modes, broadband
cmodbb:
	$(MAKE) -C src/cmodbb CModBB PETSC_ARCH=$(PETSC_ARCH_COMPLEX)

wnlrt:
	$(MAKE) -C src/wnlrt wnlrt PETSC_ARCH=$(PETSC_ARCH_REAL)

tdpape:
	$(MAKE) -C src/tdpape tdpape PETSC_ARCH=$(PETSC_ARCH_REAL)
	
#modess_rd_2wcm:
#	$(MAKE) -C src/modess_rd_2wcm ModessRD2WCM PETSC_ARCH=$(PETSC_ARCH_COMPLEX)

# clean everything up, including external libraries
clean:
	-$(MAKE) -C src/common clean  PETSC_ARCH=$(PETSC_ARCH_REAL)
	-$(MAKE) -C src/atmosphere clean  PETSC_ARCH=$(PETSC_ARCH_REAL)
	-$(MAKE) -C src/raytrace clean  PETSC_ARCH=$(PETSC_ARCH_REAL)
	-$(MAKE) -C src/modess clean  PETSC_ARCH=$(PETSC_ARCH_REAL)
	-$(MAKE) -C src/modbb clean  PETSC_ARCH=$(PETSC_ARCH_REAL)
	-$(MAKE) -C src/modess_rd_1wcm clean  PETSC_ARCH=$(PETSC_ARCH_REAL)
	-$(MAKE) -C src/pade_pe clean  PETSC_ARCH=$(PETSC_ARCH_REAL)
	-$(MAKE) -C src/wmod clean  PETSC_ARCH=$(PETSC_ARCH_REAL)
	-$(MAKE) -C src/cmodess clean  PETSC_ARCH=$(PETSC_ARCH_COMPLEX)
	-$(MAKE) -C src/cmodbb clean  PETSC_ARCH=$(PETSC_ARCH_COMPLEX)
	-$(MAKE) -C src/wnlrt clean PETSC_ARCH=$(PETSC_ARCH_REAL)
	-$(MAKE) -C src/tdpape clean PETSC_ARCH=$(PETSC_ARCH_REAL)
#	-$(MAKE) -C src/modess_rd_2wcm clean  PETSC_ARCH=$(PETSC_ARCH_COMPLEX)
	#-$(MAKE) -C src/extern clean  PETSC_ARCH=$(PETSC_ARCH_REAL)
	-rm bin/wnlrt bin/tdpape bin/CModBB  bin/CModess  bin/ModBB  bin/Modess  bin/ModessRD1WCM  bin/pape  bin/raytrace.2d  bin/raytrace.3d  bin/WMod lib/libatmosphere.a lib/libcommon.a

# Test all executables to make sure their output matches the dev machine output
test:
	$(MAKE) -C test test

# Run the tests that take forever
testlong:
	$(MAKE) -C test testlong

.PHONY: extern common atmosphere raytrace modess modbb modess_rd_1wcm pade_pe wmod cmodess cmodbb modess_rd_2wcm clean FORCE dist test

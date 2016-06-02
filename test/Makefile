# Variable check
# CXXFLAGS = -fpic -c -ggdb -Wall
# LIBS = -lgsl -lgslcblas -lm -lfftw3 
# FC = gfortran
# PETSC_VERSION = 3.2-p5
# SLEPC_VERSION = 3.2-p3
# PETSC_DIR = /home/claus/dev/ncpaprop/src/extern/petsc-3.2-p5
# SLEPC_DIR = /home/claus/dev/ncpaprop/src/extern/slepc-3.2-p3
# PETSC_OS = linux-gnu
# PETSC_ARCH_REAL = arch-linux-gnu-c-debug
# PETSC_ARCH_COMPLEX = arch-linux-gnu-c-complex
# PACKAGE_NAME = NCPA Propagation Modeling Suite
# PACKAGE_TARNAME = ncpaprop
# PACKAGE_VERSION = 1.3.17
# BASHPATH = /bin/bash
# PERLPATH = /usr/bin/perl

# Build everything
all: test

test: cleanlog calculate compare

testlong: calculate.long compare.long

.PHONY: test all testclean calculate compare calculate.raytrace.2d calculate.raytrace.3d calculate.modess calculate.cmodess calculate.modessrd1wcm calculate.wmod calculate.pape calculate.modbb calculate.cmodbb compare.raytrace.2d compare.raytrace.3d compare.modess compare.cmodess compare.modessrd1wcm compare.wmod compare.pape compare.modbb compare.cmodbb testlong compare.long calculate.long cleanlog

calculate: calculate.raytrace.2d calculate.raytrace.3d calculate.modess calculate.cmodess calculate.modessrd1wcm calculate.wmod

calculate.long: calculate.pape calculate.modbb calculate.cmodbb

compare: compare.raytrace.2d compare.raytrace.3d compare.modess compare.cmodess compare.modessrd1wcm compare.wmod

compare.long: compare.pape compare.modbb compare.cmodbb

calculate.raytrace.2d:
	@echo ""
	@echo "*** Running 2D Ray Theory Routines ***"
	@echo "<<raytrace.2d>>" > ./testlog.txt
	/bin/bash run_raytrace.2d_test.bash >> ./testlog.txt

calculate.raytrace.3d:
	@echo ""
	@echo "*** Running 3D Ray Theory Routines ***"
	@echo "<<raytrace.3d>>" >> ./testlog.txt
	/bin/bash run_raytrace.3d_test.bash >> ./testlog.txt

calculate.modess:
	@echo "" 
	@echo "*** Running Effective-Sound-Speed Modal Routines ***" 
	@echo "<<Modess>>" >> ./testlog.txt 
	/bin/bash run_Modess_test.bash >> ./testlog.txt 

calculate.cmodess:
	@echo ""
	@echo "*** Running Complex Effective-Sound-Speed Modal Routines ***"
	@echo "<<CModess>>" >> ./testlog.txt
	/bin/bash run_CModess_test.bash >> ./testlog.txt 

calculate.modessrd1wcm:
	@echo ""
	@echo "*** Running One-Way Coupled Modal Routines ***"
	@echo "<<ModessRD1WCM>>" >> ./testlog.txt
	/bin/bash run_ModessRD1WCM_test.bash >> ./testlog.txt 

calculate.wmod:
	@echo ""
	@echo "*** Running Wide-Angle Modal Routines ***"
	@echo "<<WMod>>" >> ./testlog.txt
	/bin/bash run_WMod_test.bash >> ./testlog.txt 

calculate.pape:
	@echo ""
	@echo "*** Running PE Routine - this may take a few minutes ***"
	@echo "<<pape>>" >> ./testlog.txt
	/bin/bash run_pape_test.bash >> ./testlog.txt

calculate.modbb:
	@echo ""
	@echo "*** Running Broadband Modal Routines - you may want to go get some coffee ***"
	@echo "<<ModBB>>" >> ./testlog.txt
	/bin/bash run_ModBB_test.bash >> ./testlog.txt

calculate.cmodbb:
#	@echo ""
#	@echo "*** Running Complex Broadband Modal Routines - you may want to go have dinner ***"
#	@echo "<<CModBB>>" >> ./testlog.txt
#	/bin/bash run_CModBB_test.bash >> ./testlog.txt

compare.raytrace.2d:
	@echo ""
	@echo "*** Checking 2D Ray Theory Calculation Results ***"
	/usr/bin/perl compare_raytrace.2d_test

compare.raytrace.3d:
	@echo ""
	@echo "*** Checking 3D Ray Theory Calculation Results ***"
	/usr/bin/perl compare_raytrace.3d_test

compare.modess:
	@echo ""
	@echo "*** Checking Effective-Sound-Speed Modal Calculation Results ***"
	/bin/bash compare_Modess_test

compare.cmodess:
	@echo ""
	@echo "*** Checking Complex Effective-Sound-Speed Modal Calculation Results ***"
	/bin/bash compare_CModess_test

compare.modessrd1wcm:
	@echo ""
	@echo "*** Checking One-Way Coupled Modal Calculation Results ***"
	/bin/bash compare_ModessRD1WCM_test

compare.wmod:
	@echo ""
	@echo "*** Checking Wide-Angle Modal Calculation Results ***"
	/bin/bash compare_WMod_test

compare.pape:
#	@echo ""
#	@echo "*** Checking PE Calculation Results ***"
#	/bin/bash compare_pape_test

compare.modbb:
#	@echo ""
#	@echo "*** Checking Broadmand Modal Calculation Results ***"
#	/bin/bash compare_ModBB_test

compare.cmodbb:
#	@echo ""
#	@echo "*** Checking Complex Broadband Modal Calculation Results ***"
#	/bin/bash compare_CModBB_test

#	$(MAKE) testclean

testclean:
	@rm -rf results
	@mkdir results

cleanlog:
	-rm testlog.txt

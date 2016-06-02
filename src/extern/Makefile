PETSC_VERSION=3.2-p5
PETSC_OS=linux-gnu
PETSC_DIR=/home/claus/dev/ncpaprop/src/extern/petsc-3.2-p5/
SLEPC_DIR=/home/claus/dev/ncpaprop/src/extern/slepc-3.2-p3/


extern: petsc slepc

petsc: .petsc-real .petsc-complex

slepc: .slepc-real .slepc-complex

.petsc-real: .petsc-untar
	cd petsc-3.2-p5 ; \
	./configure PETSC_ARCH=arch-linux-gnu-c-debug --with-cc=gcc --with-fc=gfortran --download-f-blas-lapack --with-shared-libraries --with-mpi=0 ; \
	cd ..
	$(MAKE) -C $(PETSC_DIR) PETSC_DIR=$(PETSC_DIR) PETSC_ARCH=arch-linux-gnu-c-debug all
	$(MAKE) -C $(PETSC_DIR) PETSC_DIR=$(PETSC_DIR) PETSC_ARCH=arch-linux-gnu-c-debug test
	touch .petsc-real

.petsc-complex: .petsc-untar
	cd petsc-3.2-p5 ; \
	./configure PETSC_ARCH=arch-linux-gnu-c-complex --with-cc=gcc --with-fc=gfortran --download-f-blas-lapack --with-cxx=g++ --with-scalar-type=complex --with-clanguage=cxx --with-fortran-kernels=generic --with-shared-libraries --with-mpi=0 ; \
	cd ..
	$(MAKE) -C $(PETSC_DIR) PETSC_DIR=$(PETSC_DIR) PETSC_ARCH=arch-linux-gnu-c-complex all
	$(MAKE) -C $(PETSC_DIR) PETSC_DIR=$(PETSC_DIR) PETSC_ARCH=arch-linux-gnu-c-complex test
	touch .petsc-complex

.slepc-real: .slepc-untar
	cd slepc-3.2-p3 ; \
	export PETSC_ARCH=arch-linux-gnu-c-debug ; \
	export PETSC_DIR=$(PETSC_DIR) ; \
	./configure
	$(MAKE) -C $(SLEPC_DIR) PETSC_ARCH=arch-linux-gnu-c-debug SLEPC_DIR=$(SLEPC_DIR) PETSC_DIR=$(PETSC_DIR)
	$(MAKE) -C $(SLEPC_DIR) PETSC_ARCH=arch-linux-gnu-c-debug SLEPC_DIR=$(SLEPC_DIR) PETSC_DIR=$(PETSC_DIR) test
	touch .slepc-real

.slepc-complex: .slepc-untar
	cd slepc-3.2-p3 ; \
	export PETSC_ARCH=arch-linux-gnu-c-complex ; \
	export PETSC_DIR=$(PETSC_DIR) ; \
	./configure
	$(MAKE) -C $(SLEPC_DIR) PETSC_ARCH=arch-linux-gnu-c-complex SLEPC_DIR=$(SLEPC_DIR) PETSC_DIR=$(PETSC_DIR)
	$(MAKE) -C $(SLEPC_DIR) PETSC_ARCH=arch-linux-gnu-c-complex SLEPC_DIR=$(SLEPC_DIR) PETSC_DIR=$(PETSC_DIR) test
	touch .slepc-complex



.petsc-untar: 
	-rm -rf petsc-3.2-p5 ; \
	tar xvzf petsc-3.2-p5.tar.gz ; \
	touch .petsc-untar

.slepc-untar:
	-rm -rf slepc-3.2-p3 ; \
	tar xvzf slepc-3.2-p3.tar.gz ; \
	touch .slepc-untar

clean:
	-rm -rf petsc-3.2-p5 slepc-3.2-p3 .petsc-untar .slepc-untar .petsc-real .petsc-complex .slepc-real .slepc-complex
	

.PHONY: petsc slepc clean

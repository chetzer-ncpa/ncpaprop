# ncpaprop

# Installation

1. Run ./configure with appropriate parameters.  Examples:

To link to an existing PETSc/SLEPc installation:

	./configure PETSC_DIR=/code/petsc SLEPC_DIR=/code/slepc PETSC_ARCH_REAL=arch-linux-c-real PETSC_ARCH_COMPLEX=arch-linux-c-complex --with-autodependencies

To download and install PETSc and SLEPc locally to the ncpaprop installation:

	./configure --with-localpetsc --with-autodependencies

See the manual for detailed information on additional parameters.

2. Run 

	make
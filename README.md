# ncpaprop

# Installation

The automake installation system is under revision and will not currently work.  To
install manually:

1. Ensure the PETSc and SLEPc packages are installed and the environment is correctly
set up (i.e. set the PETSC_DIR and SLEPC_DIR variables).  Both the real and complex
versions of PETSc and SLEPc should be built.  An example process to do this can be
found in ./install_petsc_slepc.bash.  You may edit the variables in the header of this
file and run it, but be sure to read it carefully first to understand what it is
doing.  The script is meant as a guide and error checking is minimal, so it may fail
partway through without notice.  If directories are created matching your selected
architecture names in $PETSC_DIR and $SLEPC_DIR, the build likely succeeded.  You can
double-check this by running (using architecture names arch-linux-c-real and
arch-linux-c-complex):

  `
  cd $PETSC_DIR
  make PETSC_ARCH=arch-linux-c-real check
  make PETSC_ARCH=arch-linuc-c-complex check
  cd $SLEPC_DIR
  make PETSC_ARCH=arch-linux-c-real check
  make PETSC_ARCH=arch-linuc-c-complex check
  `

  Installation problems with PETSc and SLEPc are outside the scope of this document to
  solve, and the user is directed to the documentation for those packages.

2. Edit the master Makefile in this directory and set the PETSC_ARCH_REAL and
PETSC_ARCH_COMPLEX variables to the values selected when they were built.  These
correspond to directory names in $PETSC_DIR.

3.  Run 'make'.

If build is successful, executables will be placed in bin/.  Refer to the documentation
for instructions on example model runs.

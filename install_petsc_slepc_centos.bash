#!/bin/bash

###### IMPORTANT IMPORTANT IMPORTANT IMPORTANT IMPORTANT #####
#
# This script uses commands that successfully installed and configured the PETSc
# and SLEPc libraries on the ncpaprop development system.  It is NOT meant to be used
# straight out of the box, and has been intentionally disabled to prevent this.  In
# order to make use of it, the environmental variables in the section designated as 
# "Start Customizing Here" must be set appropriately for your system.  This script is 
# provided as guidance for installation and should work properly once all variables are
# properly set, but make sure you read and understand the entire thing before you run it.
# The safest thing is to use this as guidance but run the commands manually in case of 
# unexpected error.  No warranty is provided or implied by the provision of this script.
#
# GSL, FFTW3, git, make, and g++ are prerequisites of ncpaprop.  You can use your normal 
# package manager to obtain these.
#
# You will need to have git installed to obtain PETSc.  FBLAS-LAPACK and MPICH are also
# required but can be obtained as part of the installation process if desired.  We recommend
# using your system package manager to install these beforehand.  You will need wget to
# download SLEPc, or change the command to use your command-line downloader of choice.
#
# Note: in CentOS, openmpi is currently disabled because of difficulty getting PETSc
# to play nicely with the system install of openmpi.
#
# Note: in CentOS, the prerequisites can be obtained with:
#   yum install -y blas-devel lapack-devel gsl-devel fftw-devel gcc-c++ git make
#
###### IMPORTANT IMPORTANT IMPORTANT IMPORTANT IMPORTANT #####

##### Start Customizing Here #####

# Install locations: Where the PETSc libraries should be installed.  These will be
# created if they do not already exist
PETSC_INSTALL_BASE=${HOME}/code/petsc
SLEPC_INSTALL_BASE=${HOME}/code/slepc

# Current SLEPc version, which can be checked and updated from http://slepc.upv.es/download/
SLEPC_VERSION="3.12.0"

# Operating system type (i.e. linux, darwin, etc), used to generate the unique PETSc architecture strings.
# Leave this alone to use the default system value, or customize as you wish
OS_TYPE=$(uname -s | awk '{print tolower($0)}')

# Set these to "true" if you want the PETSc installer to go get the appropriate packages for you and install
# them locally to the PETSc instance
INSTALL_FBLAS_LAPACK="false"
INSTALL_MPICH="false"

# set this to "true" to install the debug libraries, "false" to install the optimized libraries
USE_DEBUG="false"

# Set this to "true" to activate the script
I_HAVE_FINISHED_SETTING_UP="false"

##### Stop Making Changes Here #####

# exit early if the reader hasn't done everything they need to do
if [ "$I_HAVE_FINISHED_SETTING_UP" != "true" ] ; then
	echo "Script is disabled.  Go back and read all of the instructions again."
	exit
fi

# create install directory for PETSc
if [[ -d "$PETSC_INSTALL_BASE" ]] ; then
	echo "$PETSC_INSTALL_BASE already exists"
else
	echo "Creating $PETSC_INSTALL_BASE"
	mkdir -p $PETSC_INSTALL_BASE
	if [[ ! -d "$PETSC_INSTALL_BASE" ]] ; then
		echo "Error creating $PETSC_INSTALL_BASE"
		exit
	fi
fi
cd $PETSC_INSTALL_BASE
export PETSC_DIR=$PETSC_INSTALL_BASE

# download
if [[ ! -f "configure" ]] ; then
	echo "No PETSc configure script found, cloning the repository from git..."
	git clone -b maint https://gitlab.com/petsc/petsc.git .
else
	echo "Repository found, using existing repository"
fi

# set up configure commands
extra_config="--with-fc=0"
if [ "$INSTALL_FBLAS_LAPACK" = "true" ] ; then
	extra_config="$extra_config --download-fblaslapack"
fi
if [ "$INSTALL_MPICH" = "true" ] ; then
	extra_config="$extra_config --download-mpich"
else
	extra_config="$extra_config --with-mpi=0"
fi

# first, configure and make for the real version
if [ "$USE_DEBUG" = "true" ] ; then
	debug_str="-debug"
	debug_config="--with-debugging=1"
else
	debug_str=""
	debug_config="--with-debugging=0"
fi
export PETSC_ARCH="arch-${OS_TYPE}-c-real${debug_str}"

# see if it's already been made
if [[ -f "${PETSC_DIR}/${PETSC_ARCH}/include/petscconf.h" ]] ; then
	echo "PETSc architecture ${PETSC_ARCH} appears to be already built."
else
	./configure PETSC_ARCH=${PETSC_ARCH} PETSC_DIR=${PETSC_DIR} ${debug_config} ${extra_config}
	make PETSC_DIR=${PETSC_DIR} PETSC_ARCH=${PETSC_ARCH} all
fi

# now for the complex version
export PETSC_ARCH="arch-${OS_TYPE}-c-complex${debug_str}"
if [[ -f "${PETSC_DIR}/${PETSC_ARCH}/include/petscconf.h" ]] ; then
        echo "PETSc architecture ${PETSC_ARCH} appears to be already built."
else
	./configure PETSC_ARCH=${PETSC_ARCH} PETSC_DIR=${PETSC_DIR} ${debug_config} ${extra_config} --with-scalar-type=complex
	make PETSC_DIR=${PETSC_DIR} PETSC_ARCH=${PETSC_ARCH} all
fi

# create install directory for SLEPc
if [[ -d "$SLEPC_INSTALL_BASE" ]] ; then
	echo "$SLEPC_INSTALL_BASE already exists"
else
	echo "Creating $SLEPC_INSTALL_BASE"
	mkdir -p $SLEPC_INSTALL_BASE
	if [[ ! -d "$SLEPC_INSTALL_BASE" ]] ; then
		echo "Error creating $SLEPC_INSTALL_BASE"
		exit
	fi
fi
cd $SLEPC_INSTALL_BASE
export SLEPC_DIR="${SLEPC_INSTALL_BASE}/slepc-${SLEPC_VERSION}"
if [[ ! -d "$SLEPC_DIR" ]] ; then
	echo "SLEPc directory ${SLEPC_DIR} does not exist yet, downloading..."
	wget "http://slepc.upv.es/download/distrib/slepc-${SLEPC_VERSION}.tar.gz"
	tar xvzf "slepc-${SLEPC_VERSION}.tar.gz"
fi
cd $SLEPC_DIR

# configure and make real version
export PETSC_ARCH="arch-${OS_TYPE}-c-real${debug_str}"
./configure
make

# configure and make complex version
export PETSC_ARCH="arch-${OS_TYPE}-c-complex${debug_str}"
./configure
make

# provide info
echo "Installation complete"
#echo "When configuring ncpaprop for compilation, use the flags --with-petsc-real=arch-${OS_TYPE}-c-real${debug_str} and --with-petsc-complex=arch-${OS_TYPE}-c-complex${debug_str}"
echo "To make ncpaprop in one step, first edit the Makefile and set:"
echo "PETSC_ARCH_REAL = arch-${OS_TYPE}-c-real${debug_str}"
echo "PETSC_ARCH_COMPLEX = arch-${OS_TYPE}-c-complex${debug_str}"
echo "Then, type 'make all'"
echo ""
echo "Add to ~/.bashrc (or other startup script):"
echo "export PETSC_DIR=${PETSC_DIR}"
echo "export SLEPC_DIR=${SLEPC_DIR}"
echo ""
echo "Note: This installation skipped the 'make check' phase of compilation.  If you wish to do this, execute:"
echo "cd ${PETSC_DIR}"
echo "make PETSC_DIR=${PETSC_DIR} PETSC_ARCH=arch-${OS_TYPE}-c-real${debug_str} check"
echo "make PETSC_DIR=${PETSC_DIR} PETSC_ARCH=arch-${OS_TYPE}-c-complex${debug_str} check"
echo "cd ${SLEPC_DIR}"
echo "make SLEPC_DIR=${SLEPC_DIR} PETSC_DIR=${PETSC_DIR} PETSC_ARCH=arch-${OS_TYPE}-c-real${debug_str} check"
echo "make SLEPC_DIR=${SLEPC_DIR} PETSC_DIR=${PETSC_DIR} PETSC_ARCH=arch-${OS_TYPE}-c-complex${debug_str} check"
echo "If this is your first time installing PETSc and SLEPc, we highly recommend"
echo "performing this step before proceeding"

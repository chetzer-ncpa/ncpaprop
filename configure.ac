#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([NCPA Propagation Modeling Suite], [1.4.0], [claus@olemiss.edu], ncpaprop)
#AM_INIT_AUTOMAKE

AC_CONFIG_SRCDIR([src/atmosphere/Atmosphere.h])
#AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])

# Look for bash and perl.  If not found, check with the user, cause
# that's unusual for a Linux-like system
AC_PATH_PROG([BASHPATH],[bash],)
AS_IF([test "x$BASHPATH" = "x"],[
	AC_MSG_FAILURE(["No bash interpreter found on your path!  Please check this."])
])
AC_PATH_PROG([PERLPATH],[perl],)
AS_IF([test "x$PERLPATH" = "x"],[
	AC_MSG_FAILURE(["No perl interpreter found on your path!  Please check this."])
])

# Check to see if auto-prereqs are enabled
#AC_ARG_ENABLE([auto-dependencies],
#	AS_HELP_STRING([--enable-auto-dependencies],
#		[Attempt to satisfy missing dependencies automatically]
#	)
#)
#AC_ARG_ENABLE([library-guess],
#	AS_HELP_STRING([--disable-library-guess],
#		[Disable guesswork for Fink-installed library path]
#	)
#)

AC_PROG_RANLIB

AC_CANONICAL_HOST
# check host OS for support
case $host_os in
	darwin* )
		AC_MSG_NOTICE([Detected MacOS operating system])

		# for MacOS, need to know if port or fink is installed
		AS_IF([test "x$PACKMAN" = "x"],[
			AC_CHECK_PROGS([PACKMAN],[brew port fink],)
		])
		AC_MSG_NOTICE([Detected package manager = ${PACKMAN}])
		AC_SUBST([PETSC_USE_SHARED_LIBRARIES],[--with-shared-libraries=0])
		;;
	linux* )
		AC_MSG_NOTICE([Detected Linux operating system])

		# For Linux, need to know if apt-get or yum is installed
		AS_IF([test "x$PACKMAN" = "x"],[
			AC_CHECK_PROGS([PACKMAN],[apt-get yum],)
		])
		AC_MSG_NOTICE([Detected package manager = ${PACKMAN}])
		AC_SUBST([PETSC_USE_SHARED_LIBRARIES],[--with-shared-libraries])
		;;
	*)
		AC_MSG_ERROR([Your platform is not currently supported])
		;;
esac
#AC_MSG_NOTICE([enable-auto-dependencies = ${enable_auto_dependencies}])
#AC_MSG_NOTICE([enable-library-guess = ${enable_library_guess}])

# check to see if auto-dependencies are enabled and we have a package manager
#AS_IF([test "x$enable_auto_dependencies" = "xyes"],[
#	AS_CASE(["x$PACKMAN"],
#		["x"],[
#			AX_NCPA_MESSAGES_PACKMAN_NOT_FOUND
#			AC_MSG_FAILURE([No package manager found])
#		]
#	)
#	AC_MSG_NOTICE([Using ${PACKMAN} as package manager])
#])
#AS_IF([test "x$PACKMAN" = "xfink"],[
#	fullfink=$(which fink)
#	AC_MSG_NOTICE([Full path of fink is $fullfink])
#	finkbase=$(echo $fullfink | perl -e 'while(<>) { print $1 if m#(.+)/bin/fink$#; }')
#	AC_MSG_NOTICE([Using fink at $finkbase])
#])

# Checks for programs.  First, the g++ C++ compiler
#AC_CHECK_PROGS([CXX_TEST],[g++],)
#AS_IF([test "x${CXX_TEST}" = "x"], [],[
#		AX_NCPA_MESSAGES_GXX_NOT_FOUND
#	])
#])

AC_PROG_CXX
AC_PROG_CC
AC_PROG_CPP
dnl AC_PROG_INSTALL

# Checks for libraries.
AC_CHECK_LIB([fftw3], [fftw_plan_dft_1d],,[
	# Are we going to try to auto-correct?
	#AS_IF([test "x$enable_auto_dependencies" = "xyes"],[
	#	AS_CASE([$host_os],
	#		[darwin*],[
	#			mkdir prereq
	#			cd prereq
	#			curl -O http://www.fftw.org/fftw-3.3.4.tar.gz
	#			tar xvzf fftw-3.3.4.tar.gz
	#			cd fftw-3.3.4
	#			./configure
	#			make
	#			sudo make install
	#			cd ../..
	#			rm -r prereq
	#		],
	#		[linux*],[
	#			AS_CASE([$PACKMAN],
	#				[*apt-get],[
	#					sudo $PACKMAN -y install libfftw3-dev
	#				],
	#				[*yum],[
	#					sudo $PACKMAN install fftw-devel
	#				],[
	#					AX_NCPA_MESSAGES_PACKMAN_NOT_FOUND
	#					AC_MSG_FAILURE([No package manager found])
	#			])
	#	])
	#	unset ac_cv_lib_fftw3_fftw_plan_dft_1d
	#	AC_CHECK_LIB([fftw3], [fftw_plan_dft_1d],,[AX_NCPA_MESSAGES_LFFTW3_NOT_FOUND])
	#],[
		AX_NCPA_MESSAGES_LFFTW3_NOT_FOUND
	#])
])

AC_CHECK_LIB([m],[cos],,[AC_MSG_ERROR([Missing M library])])
AC_CHECK_LIB([gslcblas],[cblas_dgemm],,[
	# Are we going to try to auto-correct?
	#AS_IF([test "x$enable_auto_dependencies" = "xyes"],[
	#	AS_CASE([$host_os],
	#		[darwin*],[
	#			mkdir prereq
	#			cd prereq
	#			curl -O ftp://ftp.gnu.org/gnu/gsl/gsl-1.16.tar.gz
	#			tar xvzf gsl-1.16.tar.gz
	#			cd gsl-1.16
	#			./configure --disable-shared --disable-dependency-tracking
	#			make
	#			sudo make install
	#			cd ../..
	#			rm -r prereq
	#		],
	#		[linux*],[
	#			AS_CASE([$PACKMAN],
	#				[*apt-get],[
	#					sudo $PACKMAN -y install libgsl0-dev
	#				],
	#				[*yum],[
	#					sudo $PACKMAN install gsl-devel
	#				],[
    #                                           AX_NCPA_MESSAGES_PACKMAN_NOT_FOUND
    #                                          AC_MSG_FAILURE([No package manager found])
    #                         ])
	#		],[
	#			AC_MSG_ERROR([OS not supported])
	#	])
	#	unset ac_cv_lib_gslcblas_cblas_dgemm
	#	AC_CHECK_LIB([gslcblas],[cblas_dgemm],,[AX_NCPA_MESSAGES_LGSL_NOT_FOUND])
	#],[
		AX_NCPA_MESSAGES_LGSL_NOT_FOUND
	#])
])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm],,[
	# Are we going to try to auto-correct?
	#AS_IF([test "x$enable_auto_dependencies" = "xyes"],[
	#	AS_CASE([$host_os],
	#		[darwin*],[
	#			mkdir prereq
	#			cd prereq
	#			curl -O ftp://ftp.gnu.org/gnu/gsl/gsl-1.16.tar.gz
	#			tar xvzf gsl-1.16.tar.gz
	#			cd gsl-1.16
	#			./configure --disable-shared --disable-dependency-tracking
	#			make
	#			sudo make install
	#			cd ..
	#			rm -r prereq
	#		],
	#		[linux*],[
	#			AS_CASE([$PACKMAN],
	#				[*apt-get],[
	#					sudo $PACKMAN -y install libgsl0-dev
	#				],
	#				[*yum],[
	#					sudo $PACKMAN install gsl-devel
	#				],[
    #                                            AX_NCPA_MESSAGES_PACKMAN_NOT_FOUND
    #                                            AC_MSG_FAILURE([No package manager found])
    #                            ])
	#		],[
	#			AC_MSG_ERROR([OS not supported])
	#	])
	#	unset ac_cv_lib_gsl_gsl_blas_dgemm
	#	AC_CHECK_LIB([gsl],[gsl_blas_dgemm],,[AX_NCPA_MESSAGES_LGSL_NOT_FOUND])
	#],[
		AX_NCPA_MESSAGES_LGSL_NOT_FOUND
	#])
])

# If MacPorts is used and the selected g++ version is 4.8, try to guess the
# library path.  For some reason it doesn't get set properly when linking
# things to the PETSc libraries
#AS_IF([test "x$PACKMAN" = "xport"],[
#	AS_VAR_SET(GXX_VERSION_SHORT,$(echo ${GXX_VERSION} | perl -e '$v=<>; if ($v=~/^(\d+\.\d+)/) { print $1; }'))
#	AS_IF([test "x$GXX_VERSION_SHORT" = "x4.8"],[
#		LDFLAGS="${LDFLAGS} -L/opt/local/lib/gcc48"
#	])
#])


# Checks for header files.
AC_CHECK_HEADERS([stdlib.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FUNC_STRTOD
AC_CHECK_FUNCS([floor memset mkdir modf pow sqrt strpbrk])

# Checks for environmental variables
AS_IF([test "x$PETSC_DIR" = "x"],[
	AC_MSG_FAILURE([No value provided for PETSC_DIR!])
])
AS_IF([test "x$SLEPC_DIR" = "x"],[
	AC_MSG_FAILURE([No value provided for PETSC_DIR!])
])
AS_IF([test "x$PETSC_ARCH_REAL" = "x"],[
	AC_MSG_FAILURE([No value provided for PETSC_DIR!])
])
AS_IF([test "x$PETSC_ARCH_COMPLEX" = "x"],[
	AC_MSG_FAILURE([No value provided for PETSC_DIR!])
])

# Check the PETSc and SLEPc installations
# First, check that the directories exist
AC_CHECK_FILE( $PETSC_DIR,
	[
		AC_MSG_NOTICE([PETSc directory $PETSC_DIR OK.])
	],
	[
		AC_MSG_FAILURE([PETSc directory $PETSC_DIR not found.  Check your PETSC_DIR variable.])
	]
)
AC_CHECK_FILE( $SLEPC_DIR,
	[
		AC_MSG_NOTICE([SLEPc directory $SLEPC_DIR OK])
	],
	[
		AC_MSG_FAILURE([SLEPc directory $SLEPC_DIR not found.  Check your SLEPC_DIR variable.])
	],
)
AS_IF([test ! -d $PETSC_DIR],[
	AC_MSG_FAILURE([Specified PETSc directory ${PETSC_DIR} is not a directory!])
])
AS_IF([test ! -d $SLEPC_DIR],[
	AC_MSG_FAILURE([Specified PETSc directory ${PETSC_DIR} is not a directory!])
])

# Now see if the generic include files are present
PETSC_INCLUDE_FILE_GENERIC="${PETSC_DIR}/lib/petsc/conf/variables"
SLEPC_INCLUDE_FILE_GENERIC="${SLEPC_DIR}/lib/slepc/conf/slepc_common"
AC_CHECK_FILE( $PETSC_INCLUDE_FILE_GENERIC,
	[
		AC_MSG_NOTICE([PETSc generic include file OK])
	],
	[
		AC_MSG_FAILURE([PETSc generic include file $PETSC_INCLUDE_FILE_GENERIC not found.  Your PETSc build may be incomplete."])
	])
AC_CHECK_FILE( $SLEPC_INCLUDE_FILE_GENERIC,
	[
		AC_MSG_NOTICE([SLEPc generic include file OK])
	],
	[
		AC_MSG_FAILURE([SLEPc generic include file $SLEPC_INCLUDE_FILE_GENERIC not found.  Your SLEPc build may be incomplete.])
	])

# Finally check the architecture-specific files
AC_CHECK_FILE( ${PETSC_DIR}/${PETSC_ARCH_REAL}/lib/petsc/conf/petscvariables,
	[
		AC_MSG_NOTICE([PETSc real architecture-specific include file OK])
	],
	[
		AC_MSG_FAILURE([PETSc real architecture-specific include file ${PETSC_DIR}/${PETSC_ARCH_REAL}/lib/petsc/conf/petscvariables not found.  Check your PETSC_ARCH_REAL variable and your PETSc installation.])
	]
	)
AC_CHECK_FILE( ${PETSC_DIR}/${PETSC_ARCH_COMPLEX}/lib/petsc/conf/petscvariables,
	[
		AC_MSG_NOTICE([PETSc complex architecture-specific include file OK"])
	],
	[
		AC_MSG_FAILURE([PETSc complex architecture-specific include file ${PETSC_DIR}/${PETSC_ARCH_COMPLEX}/lib/petsc/conf/petscvariables not found.  Check your PETSC_ARCH_COMPLEX variable and your PETSc installation.])
	]
	)
AC_CHECK_FILE( ${SLEPC_DIR}/${PETSC_ARCH_REAL}/lib/slepc/conf/slepcvariables,
	[
		AC_MSG_NOTICE([SLEPc real architecture-specific include file OK])
	],
	[
		AC_MSG_FAILURE([SLEPc real architecture-specific include file ${SLEPC_DIR}/${PETSC_ARCH_REAL}/lib/slepc/conf/slepcvariables not found.  Check your PETSC_ARCH_REAL variable and your SLEPc installation.])
	]
	)
AC_CHECK_FILE( ${SLEPC_DIR}/${PETSC_ARCH_COMPLEX}/lib/slepc/conf/slepcvariables,
	[
		AC_MSG_NOTICE([SLEPc complex architecture-specific include file OK])
	],
	[
		AC_MSG_FAILURE([SLEPc complex architecture-specific include file ${SLEPC_DIR}/${PETSC_ARCH_COMPLEX}/lib/slepc/conf/slepcvariables not found.  Check your PETSC_ARCH_COMPLEX variable and your SLEPc installation.])
	]
	)




# Set output variables to propagate into Makefiles
AC_SUBST([INCLUDEFLAGS],"-I. -I../common -I../atmosphere -I/usr/local/include -I/usr/include")
AC_SUBST([LIBS],"-lgsl -lgslcblas -lm -lfftw3")
AC_SUBST([WARNINGFLAGS],"-Wwrite-strings -Wno-strict-aliasing -Wno-unknown-pragmas")
AC_SUBST([CXXFLAGS], "-fpic -c -Wall")
LDFLAGS="${LDFLAGS} -L/usr/lib"
AC_SUBST([STATICLIBS],"../../lib/libatmosphere.a ../../lib/libcommon.a")
AC_SUBST([LDFLAGS],$LDFLAGS)
AC_SUBST([PETSC_DIR],$PETSC_DIR)
AC_SUBST([PETSC_INCLUDE_FILE_GENERIC],"$PETSC_INCLUDE_FILE_GENERIC")
AC_SUBST([SLEPC_INCLUDE_FILE_GENERIC],"$SLEPC_INCLUDE_FILE_GENERIC")

AC_SUBST([PETSC_ARCH_REAL],$PETSC_ARCH_REAL)
AC_SUBST([PETSC_ARCH_COMPLEX],$PETSC_ARCH_COMPLEX)


AC_CONFIG_FILES([
		Makefile
		src/common/Makefile
		src/atmosphere/Makefile
		src/raytrace/Makefile
		src/modess/Makefile
		src/modbb/Makefile
		src/modess_rd_1wcm/Makefile
		src/pade_pe/Makefile
		src/wmod/Makefile
		src/cmodess/Makefile
		src/cmodbb/Makefile
		src/tdpape/Makefile
		src/wnlrt/Makefile
])

#		test/Makefile

AC_OUTPUT

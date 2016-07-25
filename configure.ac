#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([NCPA Propagation Modeling Suite], [1.3.18], [claus@olemiss.edu], ncpaprop)
#AM_INIT_AUTOMAKE

AC_CONFIG_SRCDIR([src/atmosphere/Atmosphere.h])
#AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])

# Look for bash and perl.  If not found, check with the user, cause
# that's unusual for a Linux-like system
AC_PATH_PROG([BASHPATH],[bash],)
AS_IF([test "x$BASHPATH" = "x"],[
	AC_MSG_FAILURE(["No bash interpreter found on your path!  Please check this.])
])
AC_PATH_PROG([PERLPATH],[perl],)
AS_IF([test "x$PERLPATH" = "x"],[
	AC_MSG_FAILURE(["No perl interpreter found on your path!  Please check this.])
])

# Check to see if auto-prereqs are enabled
AC_ARG_ENABLE([auto-dependencies],
	AS_HELP_STRING([--enable-auto-dependencies],
		[Attempt to satisfy missing dependencies automatically]
	)
)
AC_ARG_ENABLE([library-guess],
	AS_HELP_STRING([--disable-library-guess],
		[Disable guesswork for Fink-installed library path]
	)
)

AC_PROG_RANLIB

AC_CANONICAL_HOST
# check host OS for support
case $host_os in
	darwin* )
		# for MacOS, need to know if port or fink is installed
		AS_IF([test "x$PACKMAN" = "x"],[
			AC_CHECK_PROGS([PACKMAN],[port fink],)
		])
		AC_SUBST([PETSC_USE_SHARED_LIBRARIES],[--with-shared-libraries=0])
		;;
	linux* )
	
		# For Linux, need to know if apt-get or yum is installed
		AS_IF([test "x$PACKMAN" = "x"],[
			AC_CHECK_PROGS([PACKMAN],[apt-get yum],)
		])
		
		AC_SUBST([PETSC_USE_SHARED_LIBRARIES],[--with-shared-libraries])
		;;
	*)
		AC_MSG_ERROR([Your platform is not currently supported])
		;;
esac
AC_MSG_NOTICE([enable-auto-dependencies = ${enable_auto_dependencies}])
AC_MSG_NOTICE([enable-library-guess = ${enable_library_guess}])

# check to see if auto-dependencies are enabled and we have a package manager
AS_IF([test "x$enable_auto_dependencies" = "xyes"],[
#	AS_CASE(["x$PACKMAN"],
#		["x"],[
#			AX_NCPA_MESSAGES_PACKMAN_NOT_FOUND
#			AC_MSG_FAILURE([No package manager found])
#		]
#	)
	AC_MSG_NOTICE([Using ${PACKMAN} as package manager])
])
AS_IF([test "x$PACKMAN" = "xfink"],[
	fullfink=$(which fink)
	AC_MSG_NOTICE([Full path of fink is $fullfink])
	finkbase=$(echo $fullfink | perl -e 'while(<>) { print $1 if m#(.+)/bin/fink$#; }')
	AC_MSG_NOTICE([Using fink at $finkbase])
])

# Checks for programs.  First, the g++ C++ compiler
AC_CHECK_PROGS([CXX_TEST],[g++-4 g++],)
AS_IF([test "x${CXX_TEST}" = "x"], [
	# Are we going to try to auto-correct?
	AS_IF([test "x$enable_auto_dependencies" = "xyes"],[
		AS_CASE([$PACKMAN],
			[*port],[
				sudo $PACKMAN install gcc48
				sudo $PACKMAN select --set gcc mp-gcc48
			],
			[*fink],[
				sudo $PACKMAN install gcc48
			],
			[*apt-get],[
				sudo $PACKMAN -y install g++
			],
			[*yum],[
				sudo $PACKMAN install gcc-c++
			],[
				AX_NCPA_MESSAGES_PACKMAN_NOT_FOUND
				AC_MSG_FAILURE([No package manager found])
			]
		)
	],[
		AX_NCPA_MESSAGES_GXX_NOT_FOUND
	])
])

# Now fortran compiler
AC_CHECK_PROGS([FC_TEST],[gfortran],)
AS_IF([test "x${FC_TEST}" = "x"], [
	# Are we going to try to auto-correct?
	AS_IF([test "x$enable_auto_dependencies" = "xyes"],[
		AS_CASE([$PACKMAN],
			[*port],[
				sudo $PACKMAN install gcc48
				sudo $PACKMAN select --set gcc mp-gcc48
			],
			[*fink],[
				sudo $PACKMAN install gcc48
			],
			[*apt-get],[
				sudo $PACKMAN -y install gfortran
			],
			[*yum],[
				sudo $PACKMAN install gcc-gfortran
			],[
				AX_NCPA_MESSAGES_PACKMAN_NOT_FOUND
				AC_MSG_FAILURE([No package manager found])
			]
		)
	],[
		AX_NCPA_MESSAGES_GFORTRAN_NOT_FOUND
	])
])




AC_PROG_CXX([g++-4 g++])
AC_PROG_CC([gcc-4 gcc])
AC_PROG_FC([gfortran])
AC_PROG_CPP
dnl AC_PROG_INSTALL

# Check for versions of compilers to be in agreement
# This is mostly a concern for g++ and gfortran to agree
# get the g++ version
AS_VAR_SET(GXX_VERSION,$(${CXX} -dumpversion | head -1 | perl -e '$v=<>; if ($v=~/(\d+\.\d+(?:\.\d+)?)/) { print $1; }'))
if test "x$GXX_VERSION" = "x" ; then
	AC_MSG_FAILURE([g++ version cannot be determined])
else
	AC_MSG_NOTICE([g++ version ${GXX_VERSION}])
fi

# do the same check for gcc
AS_VAR_SET(GCC_VERSION,$(${CC} -dumpversion | head -1 | perl -e '$v=<>; if ($v=~/(\d+\.\d+(?:\.\d+)?)/) { print $1; }'))
if test "x$GCC_VERSION" = "x" ; then
        AC_MSG_FAILURE([gcc version cannot be determined])
else
        AC_MSG_NOTICE([gcc version ${GCC_VERSION}])
fi

# do the same check for gfortran
AS_VAR_SET(GFORTRAN_VERSION,$(${FC} -dumpversion | head -1 | perl -e '$v=<>; if ($v=~/(\d+\.\d+(?:\.\d+)?)/) { print $1; }'))
if test "x$GFORTRAN_VERSION" = "x" ; then
        AC_MSG_FAILURE([gfortran version cannot be determined])
else
        AC_MSG_NOTICE([gfortran version ${GFORTRAN_VERSION}])
fi

# check for version equality
AX_COMPARE_VERSION([$GCC_VERSION],[eq],[$GXX_VERSION])
if test "x$ax_compare_version" = "xtrue" ; then
	AC_MSG_NOTICE([gcc and g++ versions are equal])
else
	AC_MSG_FAILURE([gcc and g++ versions not equal!])
fi
AX_COMPARE_VERSION([$GCC_VERSION],[eq],[$GFORTRAN_VERSION])
if test "x$ax_compare_version" = "xtrue" ; then
	AC_MSG_NOTICE([gcc and gfortran versions are equal])
else
	AC_MSG_FAILURE([gcc and gfortran versions not equal!])
fi


# Checks for libraries.
AC_CHECK_LIB([fftw3], [fftw_plan_dft_1d],,[
	# Are we going to try to auto-correct?
	AS_IF([test "x$enable_auto_dependencies" = "xyes"],[
		AS_CASE([$host_os],
			[darwin*],[
				mkdir prereq
				cd prereq
				curl -O http://www.fftw.org/fftw-3.3.4.tar.gz
				tar xvzf fftw-3.3.4.tar.gz
				cd fftw-3.3.4
				./configure
				make
				sudo make install
				cd ../..
				rm -r prereq
			],
			[linux*],[
				AS_CASE([$PACKMAN],
					[*apt-get],[
						sudo $PACKMAN -y install libfftw3-dev
					],
					[*yum],[
						sudo $PACKMAN install fftw-devel
					],[
						AX_NCPA_MESSAGES_PACKMAN_NOT_FOUND
						AC_MSG_FAILURE([No package manager found])
				])
		])
		unset ac_cv_lib_fftw3_fftw_plan_dft_1d
		AC_CHECK_LIB([fftw3], [fftw_plan_dft_1d],,[AX_NCPA_MESSAGES_LFFTW3_NOT_FOUND])
	],[
		AX_NCPA_MESSAGES_LFFTW3_NOT_FOUND
	])
])

AC_CHECK_LIB([m],[cos],,[AC_MSG_ERROR([Missing M library])])
AC_CHECK_LIB([gslcblas],[cblas_dgemm],,[
	# Are we going to try to auto-correct?
	AS_IF([test "x$enable_auto_dependencies" = "xyes"],[
		AS_CASE([$host_os],
			[darwin*],[
				mkdir prereq
				cd prereq
				curl -O ftp://ftp.gnu.org/gnu/gsl/gsl-1.16.tar.gz
				tar xvzf gsl-1.16.tar.gz
				cd gsl-1.16
				./configure --disable-shared --disable-dependency-tracking
				make
				sudo make install
				cd ../..
				rm -r prereq
			],
			[linux*],[
				AS_CASE([$PACKMAN],
					[*apt-get],[
						sudo $PACKMAN -y install libgsl0-dev
					],
					[*yum],[
						sudo $PACKMAN install gsl-devel
					],[
                                                AX_NCPA_MESSAGES_PACKMAN_NOT_FOUND
                                                AC_MSG_FAILURE([No package manager found])
                                ])
			],[
				AC_MSG_ERROR([OS not supported])
		])
		unset ac_cv_lib_gslcblas_cblas_dgemm
		AC_CHECK_LIB([gslcblas],[cblas_dgemm],,[AX_NCPA_MESSAGES_LGSL_NOT_FOUND])
	],[
		AX_NCPA_MESSAGES_LGSL_NOT_FOUND
	])
])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm],,[
	# Are we going to try to auto-correct?
	AS_IF([test "x$enable_auto_dependencies" = "xyes"],[
		AS_CASE([$host_os],
			[darwin*],[
				mkdir prereq
				cd prereq
				curl -O ftp://ftp.gnu.org/gnu/gsl/gsl-1.16.tar.gz
				tar xvzf gsl-1.16.tar.gz
				cd gsl-1.16
				./configure --disable-shared --disable-dependency-tracking
				make
				sudo make install
				cd ..
				rm -r prereq
			],
			[linux*],[
				AS_CASE([$PACKMAN],
					[*apt-get],[
						sudo $PACKMAN -y install libgsl0-dev
					],
					[*yum],[
						sudo $PACKMAN install gsl-devel
					],[
                                                AX_NCPA_MESSAGES_PACKMAN_NOT_FOUND
                                                AC_MSG_FAILURE([No package manager found])
                                ])
			],[
				AC_MSG_ERROR([OS not supported])
		])
		unset ac_cv_lib_gsl_gsl_blas_dgemm
		AC_CHECK_LIB([gsl],[gsl_blas_dgemm],,[AX_NCPA_MESSAGES_LGSL_NOT_FOUND])
	],[
		AX_NCPA_MESSAGES_LGSL_NOT_FOUND
	])
])

# If MacPorts is used and the selected g++ version is 4.8, try to guess the
# library path.  For some reason it doesn't get set properly when linking
# things to the PETSc libraries
AS_IF([test "x$PACKMAN" = "xport"],[
	AS_VAR_SET(GXX_VERSION_SHORT,$(echo ${GXX_VERSION} | perl -e '$v=<>; if ($v=~/^(\d+\.\d+)/) { print $1; }'))
	AS_IF([test "x$GXX_VERSION_SHORT" = "x4.8"],[
		LDFLAGS="${LDFLAGS} -L/opt/local/lib/gcc48"
	])
])


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

# Set output variables to propagate into Makefiles
AC_SUBST([PETSC_OS], [$host_os])
AC_SUBST([PETSC_VERSION],[3.2-p5])
AC_SUBST([SLEPC_VERSION],[3.2-p3])
AC_SUBST([PETSC_ARCH_REAL],[arch-${host_os}-c-debug])
AC_SUBST([PETSC_ARCH_COMPLEX],[arch-${host_os}-c-complex])
AC_SUBST([PETSC_DIR],[${abs_top_builddir}/src/extern/petsc-${PETSC_VERSION}])
AC_SUBST([SLEPC_DIR],[${abs_top_builddir}/src/extern/slepc-${SLEPC_VERSION}])
AC_SUBST([CXXFLAGS], "-fpic -c -ggdb -Wall")
AC_SUBST([DORUFLAGS],"-Wwrite-strings -Wno-strict-aliasing -Wno-unknown-pragmas")
#AC_SUBST([INCLUDEFLAGS],"-I. -I../common -I../atmosphere -I/usr/local/include")
#AC_SUBST([LDFLAGS],$LDFLAGS)

AS_IF([test "x$PACKMAN" = "xfink"],[
	#AC_SUBST([LIBS], "${LIBS} -L${finkbase}/lib/gcc4.8/lib -L${finkbase}/lib")
	#AC_SUBST([LDFLAGS], "${LDFLAGS} -L${finkbase}/lib/gcc4.8/lib -L${finkbase}/lib")
	AS_IF([test "x$CXX" = "xg++-4"],[
		AS_IF([test "x$enable_library_guess" != "xno"],[
			fink_gcc_ver=$(echo $GXX_VERSION | perl -e '$v=<>; if ($v =~ /(\d+\.\d+)(?:\.\d+)?/) { print $1; }')
			AC_MSG_NOTICE([Guessing Fink library path as ${finkbase}/lib/gcc${fink_gcc_ver}/lib])
			LDFLAGS="${LDFLAGS} -L${finkbase}/lib/gcc${fink_gcc_ver}/lib"
		])
	])
	LDFLAGS="${LDFLAGS} -L${finkbase}/lib"
	AC_SUBST([INCLUDEFLAGS],"-I. -I../common -I../atmosphere -I/usr/local/include -I/sw/include")
],[
	#AC_SUBST([LDFLAGS], "${LDFLAGS} -L/usr/lib")
	LDFLAGS="${LDFLAGS} -L/usr/lib"
	AC_SUBST([INCLUDEFLAGS],"-I. -I../common -I../atmosphere -I/usr/local/include")
])
AC_SUBST([STATICLIBS],"../../lib/libatmosphere.a ../../lib/libcommon.a")
AC_SUBST([LDFLAGS],$LDFLAGS)

# todo: Make this user-configurable
AC_SUBST([PETSC_USE_MPI],[--with-mpi=0])
AC_SUBST([SLEPC_USE_MPI],[--with-mpi=0])
AC_SUBST([MPICXX],$CXX)

AC_CONFIG_FILES([
		Makefile
		src/extern/Makefile
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
		test/Makefile
])

AC_OUTPUT

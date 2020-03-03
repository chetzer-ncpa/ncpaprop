#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([NCPA Propagation Modeling Suite], [1.4.0], [claus@olemiss.edu], ncpaprop)
#AM_INIT_AUTOMAKE

AC_CONFIG_SRCDIR([src/atmosphere/Atmosphere.h])
#AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])

# set default variables
default_petsc_build_dir="${abs_builddir}/extern/petsc"
default_slepc_build_dir="${abs_builddir}/extern/slepc"
default_slepc_version="3.12.2"

# Need OS name for automatic default architecture name
os_type=$(uname -s | awk '{print tolower($0)}')
default_petsc_arch_name_real="arch-${os_type}-real"
default_petsc_arch_name_complex="arch-${os_type}-complex"




# Set up options for automatic prerequisite building
AC_ARG_WITH([local-petsc],
	AS_HELP_STRING([--with-local-petsc],[Attempt to build a PETSc instance local to the $PACKAGE_NAME installation])
)
AC_ARG_WITH([local-slepc],
	AS_HELP_STRING([--with-local-slepc],[Attempt to build a SLEPc instance local to the $PACKAGE_NAME installation.  Automatically set to 'yes' if --with-local-petsc is specified])
)
AS_IF([test "x$with_local-petsc" = "xyes"],
	[with_local-slepc="yes"] )

AC_ARG_WITH([local-petsc-debug],
	AS_HELP_STRING([--with-local-petsc-debug],[Build the debug version of PETSc/SLEPc instead of the default optimized version])
)
AC_ARG_ENABLE([local-petsc-mpi],
	AS_HELP_STRING([--enable-local-petsc-mpi],[Build PETSc/SLEPc without MPI support])
)
AC_ARG_WITH([local-petsc-mpi-download],
	AS_HELP_STRING([--with-local-petsc-mpi-download],[Download and build OpenMPI libraries as part of the PETSc installation process, rather than looking for system-installed libraries])
)
AC_ARG_WITH([local-petsc-fblas-download],
	AS_HELP_STRING([--with-local-petsc-fblas-download],[Download and build FBLAS-LAPACK libraries as part of the PETSc installation process, rather than looking for system-installed libraries])
)

# These --with options take arguments that must be processed
AC_ARG_VAR([local-petsc-dir],[Override default local PETSc installation directory (${default_petsc_build_dir})])
AC_ARG_VAR([local-slepc-dir],[Override default local SLEPc installation directory (${default_slepc_build_dir})])
AC_ARG_VAR([local-petsc-arch-real],[Override automatic PETSc real architecture name (${default_petsc_arch_name_real} or ${default_petsc_arch_name_real}-debug)])
AC_ARG_VAR([local-petsc-arch-complex],[Override automatic PETSc complex architecture name (${default_petsc_arch_name_complex} or ${default_petsc_arch_name_complex}-debug)])
AC_ARG_VAR([local-slepc-version],[Override automatic (or outdated) SLEPc version number (${default_slepc_version})])


AC_MSG_CHECKING([operating system])
AC_CANONICAL_HOST
AC_MSG_RESULT([$host_os])
# check host OS for support
case $host_os in
	darwin* )
		# for MacOS, need to know if port or fink is installed
		#AS_IF([test "x$PACKMAN" = "x"],[
		#	AC_CHECK_PROGS([PACKMAN],[brew port fink],)
		#])
		#AC_MSG_NOTICE([Detected package manager = ${PACKMAN}])
		#AC_SUBST([PETSC_USE_SHARED_LIBRARIES],[--with-shared-libraries=0])
		;;
	linux* )
		# For Linux, need to know if apt-get or yum is installed
		AC_CHECK_PROGS([default-package-manager],[apt-get yum],[])
		AC_MSG_NOTICE([Detected package manager = ${default-package-manager}])
		
		AS_IF([test "x${default-package-manager}" != "x"],[
			AC_ARG_WITH([auto-dependencies],
				AS_HELP_STRING([--with-auto-dependencies],[Attempt to install dependencies using ${default-package-manager}]))
		])
		

		#AC_SUBST([PETSC_USE_SHARED_LIBRARIES],[--with-shared-libraries])
		;;
	*)
		AC_MSG_ERROR([Your platform is not currently supported])
		;;
esac


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

AC_PROG_RANLIB
AC_PROG_CXX
AC_PROG_CC
AC_PROG_CPP
dnl AC_PROG_INSTALL

# Checks for libraries.
AC_CHECK_LIB([m],[cos],,[AC_MSG_ERROR([Missing M library])])

AC_CHECK_LIB([fftw3], [fftw_plan_dft_1d],,[
	# Are we going to try to auto-correct?
	AS_IF([test "x$enable_auto_dependencies" = "xyes"],[
		AS_CASE([$host_os],
			[darwin*],[
				AC_MSG_ERROR([Auto dependencies not available for MacOS.  Install FFTW using your package manager of choice.])
			],
			[linux*],[
				AS_CASE([$default-package-manager],
					[*apt-get],[
						sudo $default-package-manager -y install libfftw3-dev
					],
					[*yum],[
						sudo $default-package-manager install fftw-devel
					],[
						AC_MSG_ERROR([No supported package manager found (apt and yum currently supported)])
					]
				)
		])
		unset ac_cv_lib_fftw3_fftw_plan_dft_1d
		AC_CHECK_LIB([fftw3], [fftw_plan_dft_1d],,[AC_MSG_ERROR([Could not install FFTW using ${default-package-manager}, please install manually])])
	],[
		AC_MSG_ERROR([Failed to detect FFTW.  Please install the appropriate fftw3 dev package or compile from source.])
	])
])

exit

AC_CHECK_LIB([gslcblas], [cblas_dgemm],,[
	# Are we going to try to auto-correct?
	AS_IF([test "x$enable_auto_dependencies" = "xyes"],[
		AS_CASE([$host_os],
			[darwin*],[
				AC_MSG_ERROR([Auto dependencies not available for MacOS.  Install GSL using your package manager of choice.])
			],
			[linux*],[
				AS_CASE([$default-package-manager],
					[*apt-get],[
						sudo $default-package-manager -y install libgsl-dev
					],
					[*yum],[
						sudo $default-package-manager install gsl-devel
					],[
						AC_MSG_ERROR([No supported package manager found (apt and yum currently supported)])
					]
				)
		])
		unset ac_cv_lib_gslcblas_cblas_dgemm
		AC_CHECK_LIB([gslcblas], [cblas_dgemm],,[AC_MSG_ERROR([Could not install GSL using ${default-package-manager}, please install manually])])
	],[
		AC_MSG_ERROR([Failed to detect GSL.  Please install the appropriate gsl dev package or compile from source.])
	])
])

AC_CHECK_LIB([gsl], [gsl_blas_dgemm],,[
	# Are we going to try to auto-correct?
	AS_IF([test "x$enable_auto_dependencies" = "xyes"],[
		AS_CASE([$host_os],
			[darwin*],[
				AC_MSG_ERROR([Auto dependencies not available for MacOS.  Install GSL using your package manager of choice.])
			],
			[linux*],[
				AS_CASE([$default-package-manager],
					[*apt-get],[
						sudo $default-package-manager -y install libgsl-dev
					],
					[*yum],[
						sudo $default-package-manager install gsl-devel
					],[
						AC_MSG_ERROR([No supported package manager found (apt and yum currently supported)])
					]
				)
		])
		unset ac_cv_lib_gsl_gsl_blas_dgemm
		AC_CHECK_LIB([gsl], [gsl_blas_dgemm],,[AC_MSG_ERROR([Could not install GSL using ${default-package-manager}, please install manually])])
	],[
		AC_MSG_ERROR([Failed to detect GSL.  Please install the appropriate gsl dev package or compile from source.])
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

# Try to install PETSc and/or SLEPc if requested
AS_IF([test "x${with_local-petsc" = "xyes"],[
	
	# See if default PETSc directory is overridden
	AS_IF([test "x${local-petsc-dir" = "x"],
		[petsc_build_dir=${local-petsc-dir}],
		[petsc_build_dir="${abs_builddir}/extern/petsc"])

	# check for extra configuration options.  First, disable Fortran
	petsc_extra_config="--with-fc=0"
	
	# See if default PETSc real arch name is overridden
	AS_IF([test "x${local-petsc-arch-real}" = "x"],
		[
			AS_IF([test "x${with_local-petsc-debug}" = "xyes"],
				[
					petsc_local_arch_name_real="${default_petsc_arch_name_real}-debug"
				],
				[
					petsc_local_arch_name_real="${default_petsc_arch_name_real}"
				]
			)
		],[petsc_local_arch_name_real=$local-petsc-arch-real])

	# See if default PETSc complex arch name is overridden
	AS_IF([test "x${local-petsc-arch-complex}" = "x"],
		[
			AS_IF([test "x${with_local-petsc-debug}" = "xyes"],
				[
					petsc_local_arch_name_complex="${default_petsc_arch_name_complex}-debug"
				],
				[
					petsc_local_arch_name_complex="${default_petsc_arch_name_complex}"
				]
			)
		],[petsc_local_arch_name_real=$local-petsc-arch-complex])

	# debug version?
	AS_IF([test "x${with_local-petsc-debug}" = "xyes"],
		[petsc_extra_config="${petsc_extra_config} --with-debugging=1"],
		[petsc_extra_config="${petsc_extra_config} --with-debugging=0"]
	)

	# Use MPI?
	AS_IF([test "x${with_local-petsc-mpi-disabled" = "xyes"],
		[$petsc_extra_config="${petsc_extra_config} --with-mpi=0"],
		[
			AS_IF([test "x${with_local-petsc-mpi-download" = "xyes"],
				[$petsc_extra_config="${petsc_extra_config} --download-openmpi"])
		])
	
	# Go get FBLAS-LAPACK?
	AS_IF([test "x${with_local-petsc-fblas-download}" = "xyes"],
		[$petsc_extra_config="${petsc_extra_config} --download-fblaslapack"])

	# check downloader
	downloader=""
	AC_MSG_CHECKING("command-line downloader")
	AC_CHECK_PROGS([downloader],[wget curl])
	AS_IF([test "x${downloader}" = "x"],[AC_MSG_ERROR("Neither wget nor curl detected")],[AC_MSG_FOUND([$downloader])])

	# make the install dir
	mkdir -p ${petsc_build_dir}
	AS_IF([test -d "$petsc_build_dir"],,AC_MSG_ERROR(["Error creating PETSc build directory ${petsc_build_dir}"]))

	# see if it's already been downloaded and/or built
	cd ${petsc_build_dir}
	AS_IF([test -f "configure"],
		[AC_MSG_NOTICE([Configure script exists, using existing download])],
		[
			AC_MSG_NOTICE([Downloading PETSc])
			git clone -b maint https://gitlab.com/petsc/petsc.git .
		]
	)

	AC_MSG_CHECKING([for existing ${petsc_local_arch_name_real} architecture])
	AC_CHECK_FILE( ${petsc_build_dir}/${petsc_local_arch_name_real}/lib/petsc/conf/petscvariables,
		[AC_MSG_RESULT([yes])],
		[
			AC_MSG_RESULT([no, building])

			# build the real version with indicated options
			./configure PETSC_ARCH=${petsc_local_arch_name_real} PETSC_DIR=${petsc_build_dir} ${petsc_extra_config}
			AS_IF([test $? != 0],
				[AC_MSG_ERROR([Error running configure command: ./configure PETSC_ARCH=${petsc_local_arch_name_real} PETSC_DIR=${petsc_build_dir} ${petsc_extra_config}])]
			)
			make all
			AS_IF([test $? != 0],
				[AC_MSG_ERROR([Error running make command for real PETSc: make all])]
			)

			# see if it made the include file
			AC_CHECK_FILE( ${petsc_build_dir}/${petsc_local_arch_name_real}/lib/petsc/conf/petscvariables,
				[
					AC_MSG_NOTICE([PETSc real version build OK])
				],
				[
					AC_MSG_ERROR([Error building PETSc real version])
				]
			)
		]
	)

	AC_MSG_CHECKING([for existing ${petsc_local_arch_name_complex} architecture])
	AC_CHECK_FILE( ${petsc_build_dir}/${petsc_local_arch_name_complex}/lib/petsc/conf/petscvariables
		[AC_MSG_RESULT([yes])],
		[
			AC_MSG_RESULT([no, building])

			# build the complex version with indicated options
			./configure PETSC_ARCH=${petsc_local_arch_name_complex} PETSC_DIR=${petsc_build_dir} ${petsc_extra_config} --with-scalar-type=complex
			AS_IF([test $? != 0],
				[AC_MSG_ERROR([Error running configure command: ./configure PETSC_ARCH=${petsc_local_arch_name_complex} PETSC_DIR=${petsc_build_dir} ${petsc_extra_config} --with-scalar-type=complex])]
			)
			make all
			AS_IF([test $? != 0],
				[AC_MSG_ERROR([Error running make command for complex PETSc: make all])]
			)

			# see if it made the include file
			AC_CHECK_FILE( ${petsc_build_dir}/${petsc_local_arch_name_complex}/lib/petsc/conf/petscvariables,
				[
					AC_MSG_NOTICE([PETSc complex version build OK])
				],
				[
					AC_MSG_ERROR([Error building PETSc complex version])
				]
			)
		]
	)

	PETSC_DIR=${petsc_build_dir}
	PETSC_ARCH_REAL=${petsc_local_arch_name_real}
	PETSC_ARCH_COMPLEX=${petsc_local_arch_name_complex}
])









# Checks for environmental variables
AS_IF([test "x$PETSC_DIR" = "x"],[
	AC_MSG_FAILURE([No value provided for PETSC_DIR!])
])
AS_IF([test "x$SLEPC_DIR" = "x"],[
	AC_MSG_FAILURE([No value provided for SLEPC_DIR!])
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

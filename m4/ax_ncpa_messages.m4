# AX_NCPA_MESSAGES_PACKMAN_NOT_FOUND
# Outputs a customized error message with instructions for
# dealing with no package manager
#
AC_DEFUN([AX_NCPA_MESSAGES_PACKMAN_NOT_FOUND],[
	AC_REQUIRE([AC_CANONICAL_HOST])
	AS_CASE([$host_os],
		[linux*],[
			cat <<EOF

Configure failed to detect either apt-get or yum, the two package managers supported by the automatic installation option.  If you are using a package manager other than these two, you will need to manually install the proper packages (or install from source) in order to have:

1.  Matching versions of gcc, g++, and gfortran
2.  FFTW version 3 (including developer libraries)
3.  GSL (including developer libraries)

EOF
		],
		[darwin*],[
			cat <<EOF

Configure failed to detect either port or fink, the two package managers supported by the automatic installation option.  If you are using a package manager other than these two, you will need to manually install the proper packages (or install from source) in order to have:

1.  Matching versions of gcc, g++, and gfortran

EOF
		],[
			AC_MSG_ERROR([No supported OS detected!])
	])
])


# AX_NCPA_MESSAGES_GXX_NOT_FOUND
# Outputs a customized error messages with instructions for 
# installing g++ depending on the detected OS.
#
AC_DEFUN([AX_NCPA_MESSAGES_GXX_NOT_FOUND],[
	AC_REQUIRE([AC_CANONICAL_HOST])
	case $host_os in
		linux* )
			cat <<EOF

Configure failed to detect g++, the GNU C++ compiler.  This is required to compile $PACKAGE_NAME.  Please install the 'g++' package or compile from source.  Depending on your package manager, this may be accomplished with 
   sudo apt-get install g++
or
   sudo yum install gcc-c++

EOF
			AC_MSG_ERROR([missing g++ compiler])
			;;

		darwin* )
			cat <<EOF

Configure failed to detect g++, the GNU C++ compiler.  This is required to compile $PACKAGE_NAME.  To obtain Apple's base version, install the XCode package from the App Store, then from within XCode install the Command Line Developer Tools.  This will provide both gcc and g++, most likely version 4.2.1.  However, it will not install the gfortran Fortran compiler.  To install this, and simultaneously a more up-to-date (and supported) version of gcc and g++, we recommend using MacPorts:

1.  Make sure you've accepted the XCode license agreement using
    xcodebuild -license

2.  Download and install the latest version of MacPorts from macports.org

3.  Install and activate the Gnu compilers version 4.8 using
    sudo port install gcc48
    sudo port select --set gcc mp-gcc48


EOF
			AC_MSG_ERROR([missing g++ compiler])
			;;

		*)
			AC_MSG_ERROR([operating system not recognized.])
			;;
	esac
])


# AX_NCPA_MESSAGES_GFORTRAN_NOT_FOUND
# Outputs a customized error message with instructions for installing
# gfortran, depending on the detected OS.
#
AC_DEFUN([AX_NCPA_MESSAGES_GFORTRAN_NOT_FOUND],[
	AC_REQUIRE([AC_CANONICAL_HOST])
	AS_VAR_SET(GCC_VERSION,$(${CC} -dumpversion))
	case $host_os in
		linux* )
			cat <<EOF

Configure failed to detect gfortran, the GNU Fortran compiler.  This is required to compile $PACKAGE_NAME.  Please install the gfortran package or compile from source.  Depending on your package manager, this may be accomplished with 
   sudo apt-get install gfortran
or
   sudo yum install gcc-gfortran
or you may install from source.  Please be careful to match the version to that of gcc, which is version ${GCC_VERSION}.

EOF
			AC_MSG_ERROR([missing gfortran compiler])
			;;

		darwin* )
			cat <<EOF

Configure failed to detect gfortran, the GNU Fortran compiler.  This is required to compile $PACKAGE_NAME.  Please install the gfortran package or compile from source.  Please be sure to match the version to that of gcc, which is version ${GCC_VERSION}.

To simultaneously install matching versions of gcc and g++, we recommend using MacPorts:

1.  Make sure you've accepted the XCode license agreement using
    xcodebuild -license

2.  Download and install the latest version of MacPorts from macports.org

3.  Install and activate the Gnu compilers version 4.8 using
    sudo port install gcc48
    sudo port select --set gcc mp-gcc48


EOF
			AC_MSG_ERROR([missing gfortran compiler])
			;;

		*)
			AC_MSG_ERROR([operating system not supported.])
			;;
	esac
])

# AX_NCPA_MESSAGES_LFFTW3_NOT_FOUND
# Outputs a customized error message with instructions for installing
# fftw3, depending on the detected OS.
#
AC_DEFUN([AX_NCPA_MESSAGES_LFFTW3_NOT_FOUND],[
        AC_REQUIRE([AC_CANONICAL_HOST])
        case $host_os in
                linux* )
                        cat <<EOF

Configure failed to detect libfftw3, the FFTW version 3 library.  This is required to compile $PACKAGE_NAME.  Please install the fftw3 dev package or compile from source.  Depending on your package manager, this may be accomplished with
   sudo apt-get install libfftw3-dev
or
   sudo yum install fftw-devel

EOF
                        AC_MSG_ERROR([missing libfftw3])
                        ;;

                darwin* )
                        cat <<EOF

Configure failed to detect libfftw3, the FFTW version 3 library.  This is required to compile $PACKAGE_NAME.  At time of release the current version is 3.3.4.  Please download and compile with:

curl -O http://www.fftw.org/fftw-3.3.4.tar.gz
tar xvzf fftw-3.3.4.tar.gz
cd fftw-3.3.4
./configure
make
sudo make install


EOF
                        AC_MSG_ERROR([missing libfftw3])
                        ;;

                *)
                        AC_MSG_ERROR([operating system not supported.])
                        ;;
        esac
])

# AX_NCPA_MESSAGES_LGSL_NOT_FOUND
# Outputs a customized error message with instructions for installing
# GSL, depending on the detected OS.
#
AC_DEFUN([AX_NCPA_MESSAGES_LGSL_NOT_FOUND],[
        AC_REQUIRE([AC_CANONICAL_HOST])
        case $host_os in
                linux* )
                        cat <<EOF

Configure failed to detect libgsl, the GSL library.  This is required to compile $PACKAGE_NAME.  Please install the GSL dev package or compile from source.  Depending on your package manager, this may be accomplished with
   sudo apt-get install libgsl0-dev
or
   sudo yum install gsl-devel

EOF
                        AC_MSG_ERROR([missing libgsl])
                        ;;

                darwin* )
                        cat <<EOF

Configure failed to detect libgsl, the GSL library.  This is required to compile $PACKAGE_NAME.  At release time, the current version was 1.16.  Please compile from source:

curl -O ftp://ftp.gnu.org/gnu/gsl/gsl-1.16.tar.gz
tar xvzf gsl-1.16.tar.gz
cd gsl-1.16
./configure
make
sudo make install


EOF
                        AC_MSG_ERROR([missing libgsl])
                        ;;

                *)
                        AC_MSG_ERROR([operating system not supported.])
                        ;;
        esac
])


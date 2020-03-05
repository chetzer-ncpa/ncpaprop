AC_DEFUN([AC_CXX_MPI], [
AC_PREREQ(2.50) dnl for AC_LANG_CASE

AC_REQUIRE([AC_PROG_CXX])
AC_ARG_VAR(MPICXX,[MPI C++ compiler command])

AC_CHECK_PROGS(MPICXX, mpic++ mpicxx mpiCC hcp mpxlC_r mpxlC mpCC cmpic++, $CXX)
if test $MPICXX = $CXX; then
 AC_ERROR([cannot find MPI])
fi

        acx_mpi_save_CXX="$CXX"
        CXX="$MPICXX"
        AC_SUBST(MPICXX)

if test x = x"$MPILIBS"; then
#check for mpi header
dnl We have to use AC_TRY_COMPILE and not AC_CHECK_HEADER because the
dnl latter uses $CPP, not $CC (which may be mpicc).
    AC_MSG_CHECKING([for mpi++.h])
    AC_TRY_COMPILE([#include <mpi++.h>],[],[AC_MSG_RESULT(yes)], [MPILIBS=""
                AC_MSG_RESULT(no) 
                AC_ERROR([cannot find mpi++.h])])
#check for mpi libs
  AC_CHECK_FUNC(MPI_Init, [MPILIBS=" "])
  AC_CHECK_LIB(mpi, MPI_Init, [MPILIBS="-lmpi"])
  AC_CHECK_LIB(mpich, MPI_Init, [MPILIBS="-lmpich"])
fi

dnl CXX="$acx_mpi_save_CXX"

AC_SUBST(MPILIBS)

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x = x"$MPILIBS"; then
        $2
        :
else
        ifelse([$1],,[AC_DEFINE(HAVE_MPI,1,[Define if you have the MPI library.])],[$1])
        :
fi
])dnl ACX_MPI

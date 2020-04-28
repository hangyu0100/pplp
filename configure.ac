AC_INIT([pplp], [1.0])
AC_CONFIG_MACRO_DIRS([m4])
AM_INIT_AUTOMAKE([-Wall -Wno-extra-portability foreign subdir-objects])

AC_PROG_CPP
AC_PROG_CXX
AC_PROG_LIBTOOL

AC_ARG_ENABLE([parallel],
  AC_HELP_STRING([--enable-parallel],[tbb, openmp or disable [default=tbb]]),
  [case $enableval in 
    tbb | openmp | disable) ;;
    *) AC_MSG_ERROR([Invalid value for --parallel]) ;;
  esac],
  [enable_parallel=tbb]
)

AC_MSG_CHECKING([library for parallelism])
if test "x$enable_parallel" = "xopenmp"; then
  CPPFLAGS="$CPPFLAGS -fopenmp"
  LIBS="-fopenmp $LIBS"
  AC_MSG_RESULT([openmp])
elif test "x$enable_parallel" = "xtbb"; then
  CPPFLAGS="$CPPFLAGS -D_TBB"
  LIBS="-ltbb $LIBS"
  AC_MSG_RESULT([tbb])
else
  AC_MSG_RESULT([disable])
fi 

AC_SEARCH_LIBS([exp], [m], [], [AC_MSG_ERROR([Cannot find m library.])])
AC_SEARCH_LIBS([glp_simplex], [glpk], [], [AC_MSG_ERROR([Cannot find glpk library.])])
AC_SEARCH_LIBS([fmpq_mat_init], [flint], [], [AC_MSG_ERROR([Cannot find flint library.])])
LIBS="-lgmp -lgmpxx -lmpfr $LIBS"

AC_CONFIG_FILES([
 src/Makefile
 test/Makefile
 Makefile
])
AC_OUTPUT

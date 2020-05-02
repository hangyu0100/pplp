sudo apt install uswsusp  AC_INIT([pplp], [1.0])
AC_CONFIG_MACRO_DIRS([m4])
AM_INIT_AUTOMAKE([-Wall -Wno-extra-portability foreign subdir-objects])

AC_PROG_CPP
AC_PROG_CXX
AC_PROG_LIBTOOL

AC_ARG_ENABLE([parallel],
  AC_HELP_STRING([--enable-parallel], 
                 [tbb, openmp or disable [default=tbb]]),
  [case $enableval in 
    tbb | openmp | disable) ;;
    *) AC_MSG_ERROR([Invalid value for --enable-parallel]) ;;
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

AC_ARG_ENABLE([rational-checker],
  AC_HELP_STRING([--enable-rational-checker], 
                 [turn on the rational checker [default=yes]]),
  [case $enableval in
    yes | no) ;;
    *) AC_MSG_ERROR([Invalid value for --enable-rational-checker])
  esac],
  [enable_rational_checker=yes]
)
AC_MSG_CHECKING([whether to use rational checker])
if test "x$enable_rational_checker" = "xyes"; then
  CPPFLAGS="$CPPFLAGS -DPLP_CHECKER_ON"
  AC_MSG_RESULT([yes])
else
  AC_MSG_RESULT([no])
fi

AC_ARG_ENABLE([debug-log],
  AC_HELP_STRING([--enable-debug-log],
                 [choose level of debug logs [default=0, turn off]]),
  [case $enableval in
    0 | 1 | 2 | 3 | 4 | 5) ;;
    *) AC_MSG_ERROR([Invalid value for --enable-debug-log])
  esac],
  [enable_debug_log=0]
)
AC_MSG_CHECKING([the level of debug logs])
if test "x$enable_debug_log" = "x1"; then
  CPPFLAGS="$CPPFLAGS -DPRINT_WARNING -DPRINT_WARNING"
  AC_MSG_RESULT([1, print warnings])
elif test "x$enable_debug_log" = "x2"; then
  CPPFLAGS="$CPPFLAGS -DPRINT_WARNING -DDEBUGINFO_RAYTRACING"
  AC_MSG_RESULT([2, print debug logs for raytracing minimization])
elif test "x$enable_debug_log" = "x3"; then
  CPPFLAGS="$CPPFLAGS -DPRINT_WARNING -DDEBUGINFO_SIMPLEX"
  AC_MSG_RESULT([3, print debug logs for ratinoal Simplex algorithm])
elif test "x$enable_debug_log" = "x4"; then
  CPPFLAGS="$CPPFLAGS -DPRINT_WARNING -DDEBUGINFO_PLP"
  AC_MSG_RESULT([4, print debug logs for PLP])
elif test "x$enable_debug_log" = "x5"; then
  CPPFLAGS="$CPPFLAGS -DPRINT_WARNING -DDEBUGINFO_RAYTRACING -DDEBUGINFO_SIMPLEX -DDEBUGINFO_PLP"
  AC_MSG_RESULT([5, print all debug logs])
else
  AC_MSG_RESULT([0, turn off])
fi

AC_ARG_ENABLE([rational-plp],
  AC_HELP_STRING([--enable-rational-plp], 
                 [use ratinoal PLP algorithm instead of floating-point PLP
                  [default=no]]),
  [case $enableval in
    yes | no) ;;
    *) AC_MSG_ERROR([Invalid value for --enable-ratinoal-plp])
  esac],
  [enable_rational_plp=no]
)
AC_MSG_CHECKING([whether to use ratinoal PLP algorithm])
if test "x$enable_rational_plp" = "xyes"; then
  CPPFLAGS="$CPPFLAGS -DALWAYS_RATIONAL"
  AC_MSG_RESULT([yes])
else
  AC_MSG_RESULT([no])
fi

AC_ARG_WITH(apron,
  AC_HELP_STRING([--with-apron], 
               [use apron library for comparing the results [default=no]]),
  [case $withval in
    yes | no) ;;
    *) AC_MSG_ERROR(Invalid value for --with-apron) 
  esac],
  [with_apron=no]
)
AC_MSG_CHECKING([whether to use apron library])
if test "x$with_apron" = "xyes"; then
  LIBS="$LIBS -lpolkaMPQ -lapron"
  AC_MSG_RESULT([yes])
else
  AC_MSG_RESULT([no])
fi
AM_CONDITIONAL([WITH_APRON],[test "x$with_apron" = "xyes"])

AC_SEARCH_LIBS([exp], [m], [], [AC_MSG_ERROR([Cannot find m library.])])
AC_SEARCH_LIBS([glp_simplex], [glpk], [], [AC_MSG_ERROR([Cannot find glpk library.])])
AC_SEARCH_LIBS([fmpq_mat_init], [flint], [], [AC_MSG_ERROR([Cannot find flint library.])])
LIBS="-lgmp -lgmpxx -lmpfr $LIBS"

WARNING="-Wall -Wextra -pedantic -Wno-deprecated-copy"
CXXFLAGS="$CXXFLAGS $WARNING -std=c++11 -g -fPIC"

AC_CONFIG_FILES([
 src/Makefile
 test/Makefile
 Makefile
])
AC_OUTPUT

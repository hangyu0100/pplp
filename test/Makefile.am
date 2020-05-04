AM_CPPFLAGS = -I ../include -I /usr/include/flint
LDADD = ../src/libpplp.la 

bin_PROGRAMS = minimize projection convexhull
minimize_SOURCES = test_minimize.cc
projection_SOURCES = test_projection.cc
convexhull_SOURCES = test_convexhull.cc

if WITH_APRON
bin_PROGRAMS += minimize_apron convexhull_apron projection_apron
minimize_apron_SOURCES = test_minimize_apron.cc
convexhull_apron_SOURCES = test_convexhull_apron.cc
projection_apron_SOURCES = test_projection_apron.cc
endif

check_PROGRAMS = unit_test
TESTS = ${check_PROGRAMS}
unit_test_SOURCES = unit_test.cc


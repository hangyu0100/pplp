
AM_CPPFLAGS = -I ../include -I /usr/include/flint
LDADD = ../src/libpplp.la 

bin_PROGRAMS = apron convexhull_apron equality exact inclusion minimize \
  minimize_simple plp projection_apron projection \
  witness unit_test
apron_SOURCES = test_apron.cc

convexhull_apron_SOURCES = test_convexhull_apron.cc

equality_SOURCES = test_equality.cc

exact_SOURCES = test_exact.cc

inclusion_SOURCES = test_inclusion.cc

minimize_SOURCES = test_minimize.cc

minimize_simple_SOURCES = test_minimize_simple.cc

plp_SOURCES = test_plp.cc

projection_apron_SOURCES = test_projection_apron.cc

projection_SOURCES = test_projection.cc

witness_SOURCES = test_witness.cc

unit_test_SOURCES = unit_test.cc
check: unit_test$(EXEEXT)

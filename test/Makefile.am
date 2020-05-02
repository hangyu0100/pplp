AM_CPPFLAGS = -I ../include -I /usr/include/flint
LDADD = ../src/libpplp.la 

bin_PROGRAMS = equality exact inclusion minimize \
  minimize_simple projection witness unit_test

if WITH_APRON
bin_PROGRAMS += apron convexhull_apron projection_apron plp 
apron_SOURCES = test_apron.cc
convexhull_apron_SOURCES = test_convexhull_apron.cc
projection_apron_SOURCES = test_projection_apron.cc
plp_SOURCES = test_plp.cc
endif

equality_SOURCES = test_equality.cc
exact_SOURCES = test_exact.cc
inclusion_SOURCES = test_inclusion.cc
minimize_SOURCES = test_minimize.cc
minimize_simple_SOURCES = test_minimize_simple.cc
projection_SOURCES = test_projection.cc
witness_SOURCES = test_witness.cc
unit_test_SOURCES = unit_test.cc
check: unit_test$(EXEEXT)

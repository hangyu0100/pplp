lib_LTLIBRARIES = libpplp.la

libpplp_la_CPPFLAGS = -I../include -I/usr/include/flint
libpplp_la_LDFLAGS = -L/usr/lib -L/usr/local/lib

libpplp_la_SOURCES = \
core/pplp_double.cc \
core/pplp_ioInterface.cc \
core/pplp_point.cc \
core/pplp_polyhedron.cc \
core/pplp_tool.cc \
plp/pplp_plp.cc \
plp/pplp_simplex.cc \
plp/pplp_tbbParallel.cc \
raytracing/pplp_ray.cc \
raytracing/pplp_raytracing.cc \
thirdpart/pplp_apronInterface.cc \
thirdpart/pplp_glpkInterface.cc

include_HEADERS = \
../include/pplp_apronInterface.h \
../include/pplp_double.h \
../include/pplp_glpkInterface.h \
../include/pplp_ioInterface.h \
../include/pplp_plp.h \
../include/pplp_point.h \
../include/pplp_polyhedron.h \
../include/pplp_ray.h \
../include/pplp_raytracing.h \
../include/pplp_simplex.h \
../include/pplp_simple_concurrent_vector.h \
../include/pplp_tbbParallel.h \
../include/pplp_tool.h

libpplp_la_LIBADD = -lglpk -lm  -lpolkaMPQ -lapron -lgmpxx -lflint -lmpfr -lgmp

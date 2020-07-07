Install: 
  ./configure
  make
  make install
  
Required libraries: GLPK, TBB, GMP (with c++), MPFR, FLINT, APRON (optional)

To use this library:
Compile with 
  LDLIBS+=-lpplp -ltbb -lglpk -lm  -lpolkaMPQ -lapron -lgmpxx -lflint -lmpfr -lgmp
  CPPFLAGS+=-D _TBB (if use TBB)

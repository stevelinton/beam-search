bin_PROGRAMS = addchain ascode aascode addchain2

AM_CFLAGS = -g -O3 -march=native $(OPENMP_CFLAGS)

BEAM = beam.c beam.h
addchain_SOURCES = addchain.c $(BEAM)
addchain2_SOURCES = addchain2.c $(BEAM)
ascode_SOURCES = ascode.c $(BEAM)
aascode_SOURCES = aascode.c $(BEAM)

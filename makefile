# Linux/OSX makefile for spharm-interp MEX compilation.
# cut down from cryo work 2015-2016; Alex Barnett.

# ===== system-dependent settings : please adjust =====
FC = gfortran
FFLAGS = -O3 -march=native -std=legacy

# No openmp (single-core):
#OMPFLAGS =
#OMBLIBS =
# Or use following for openmp: (OSX may instead need link to libiomp5.so)
OMPFLAGS = --openmp
OMPLIBS = -lgomp

# Linux users, set this to your mex executable (eg in MATLAB source tree):
MEX = mex
# Mac users should use something like this:
#MEX = /Applications/MATLAB_R2014a.app/bin/mex

# Linux users:
MEXFLAGS = -largeArrayDims -lgfortran -lm $(OMPLIBS)
# Mac users use this:
#MEXFLAGS = -largeArrayDims -L/usr/local/gfortran/lib -lgfortran -lm
# =====================================

# Fortran sources
FSRCS = spharmbasics.f yrecursion.f
OBJS    = $(FSRCS:.f=.o)

default: mex

.f.o:
	$(FC) -fPIC $(FFLAGS) $(OMPFLAGS) -c $<

# Compile the interface, using whatever libs needed:
mex: gateway.c $(OBJS)
	$(MEX) gateway.c $(OBJS) $(MEXFLAGS)

# rebuild gateway, experts only (requires mwrap to be installed):
gateway.c: spharm.mw
	mwrap -c99complex -mex gateway -mb -list spharm.mw
	mwrap -c99complex -mex gateway -c gateway.c spharm.mw

# Remove all generated files (except gateway.c which is part of the repo):
clean:
	rm -f *.o *.mex*

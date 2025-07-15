# Build a LBFGS-B library:

FC = gfortran
FILES = blas.f \
	lbfgsb.f \
	linpack.f \
	timer.f

LIB   = liblbfgs.a

RANLIB = ranlib

ARCHIVE = ../libext/$(LIB)($(addsuffix .o, $(basename $(FILES))))

all: lib

lib: $(LIB)
$(LIB): $(ARCHIVE)
	$(RANLIB) ../libext/$@

clean:
	rm -f liblbfgs.a *.o

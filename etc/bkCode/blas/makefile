# Build an adjoint library:

AR=ar r
FC=gfortran
LINK=gfortran
FLAG=-O3
INC=.

BLAS_SOURCE=\
	dasum.f \
	daxpy.f \
	dcopy.f \
	ddot.f \
	dnrm2.f \
	drot.f \
	drotg.f \
	dscal.f \
	dswap.f \
	idamax.f \
	isamax.f \
	sasum.f \
	saxpy.f \
	sdot.f \
	sscal.f \
	sswap.f

lib_date: $(BLAS_SOURCE)
	$(FC) $(FLAG) -c -I$(INC) $?
	$(AR) ../lib/liblas.a $(?:.f=.o)
	/bin/rm $(?:.f=.o)
	touch $@

clean:
	rm -f lib_date

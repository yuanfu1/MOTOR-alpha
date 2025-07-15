# Build an adjoint library:

AR=ar r
FC=gfortran
LINK=gfortran
FLAG=-O3
INC=.

CH_SOURCE=\
	dqrdc.f \
	dqrsl.f \

lib_date: $(CH_SOURCE)
	$(FC) $(FLAG) -c -I$(INC) $?
	$(AR) ../lib/liblinpack.a $(?:.f=.o)
	/bin/rm $(?:.f=.o)
	touch $@

clean:
	rm -f lib_date *.o


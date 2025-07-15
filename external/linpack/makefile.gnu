# Build an adjoint library:

AR=ar r
FC=gfortran
LINK=gfortran
FLAG=-O3
INC=.

CH_SOURCE=\
	dgeco.f \
	dgefa.f \
	dgesl.f \
	dqrdc.f \
	dqrsl.f \
	schdc.f \
	sgbco.f \
	sgbfa.f \
	sgbsl.f \
	sgeco.f \
	sgedi.f \
	sgefa.f \
	sgesl.f \
	spoco.f \
	spofa.f \
	sqrdc.f \
	sqrsl.f

lib_date: $(CH_SOURCE)
	$(FC) $(FLAG) -c -I$(INC) $?
	$(AR) ./liblinpack.a $(?:.f=.o)
	/bin/rm $(?:.f=.o)
	touch $@

clean:
	rm -f lib_date *.o *.a


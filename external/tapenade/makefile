AR=ar -r
FC=gfortran-9
LibName=libtapenadeWRF

tapenade_c=\
	adStack.c \
	adBuffer.F \

lib_data: $(tapenade_c)
	$(FC) -c $?
	$(AR) ./$(LibName).a *.o
	/bin/rm *.o
	touch $@

clean:
	rm -f lib_data *.a

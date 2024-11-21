CPP=g++
LINKER=g++
FFLAGS=-O
LFLAGS=
OBJECTFILES=integrate.o

default: integrate Sq_Nwt Sq_Frc

integrate: integrate.o
	echo "Compiling integrate"; \
	$(LINKER) $(LFLAGS) integrate.o -o integrate

integrate.o: integrate.cc
	$(CPP) $(FFLAGS) integrate.cc -c

clean:
	rm *.o integrate
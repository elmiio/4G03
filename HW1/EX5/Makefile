CPP=g++
LINKER=g++
FFLAGS=-O
LFLAGS=
OBJECTFILES=Sq_Tlr.o Sq_Nwt.o Sq_Frc.o timing.o

default: Sq_Tlr Sq_Nwt Sq_Frc

Sq_Tlr: Sq_Tlr.o
	echo "Compiling Sq_Tlr"; \
	$(LINKER) $(LFLAGS) Sq_Tlr.o -o Sq_Tlr
 
Sq_Nwt: Sq_Nwt.o
	echo "Compiling Sq_Nwt"; \
	$(LINKER) $(LFLAGS) Sq_Nwt.o -o Sq_Nwt
 
Sq_Frc: Sq_Frc.o
	echo "Compiling Sq_Frc"; \
	$(LINKER) $(LFLAGS) Sq_Frc.o -o Sq_Frc

Sq_Tlr.o: Sq_Tlr.cc
	$(CPP) $(FFLAGS) Sq_Tlr.cc -c

Sq_Nwt.o: Sq_Nwt.cc
	$(CPP) $(FFLAGS) Sq_Nwt.cc -c

Sq_Frc.o: Sq_Frc.cc
	$(CPP) $(FFLAGS) Sq_Frc.cc -c

clean:
	rm *.o Sq_Tlr Sq_Nwt Sq_Frc
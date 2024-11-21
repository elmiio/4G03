CPP=g++
LINKER=g++
FFLAGS=-O
LFLAGS=
OBJECTFILES_3_1=3_1.o
OBJECTFILES_3_3=3_3.o

default: 3_1 3_3

3_1: $(OBJECTFILES_3_1)
	echo "Compiling"; \
	$(LINKER) $(LFLAGS) $(OBJECTFILES_3_1) -o 3_1
 
3_3: $(OBJECTFILES_3_3)
	echo "Compiling"; \
	$(LINKER) $(LFLAGS) $(OBJECTFILES_3_3) -o 3_3

$(OBJECTFILES): 3_1.cc
	$(CPP) $(FFLAGS) 3_1.cc -c
 
 $(OBJECTFILES): 3_3.cc
	$(CPP) $(FFLAGS) 3_3.cc -c

clean:
	rm *.o 3_1 3_3
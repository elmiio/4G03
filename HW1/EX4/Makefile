CPP=g++
LINKER=g++
FFLAGS=-O
LFLAGS=
OBJECTFILES=zeros.o

default: zeros

zeros: $(OBJECTFILES)
	echo "Compiling"; \
	$(LINKER) $(LFLAGS) $(OBJECTFILES) -o zeros

$(OBJECTFILES): zeros.cc
	$(CPP) $(FFLAGS) zeros.cc -c

clean:
	rm *.o zeros
CPP=g++
LINKER=g++
FFLAGS=-O -std=c++11
LFLAGS=-lm
OBJECTFILES=mcarlo.o

default: mcarlo

mcarlo: $(OBJECTFILES)
	echo "Compiling"; \
	$(LINKER) $(LFLAGS) $(OBJECTFILES) -o mcarlo

$(OBJECTFILES): mcarlo.cc
	$(CPP) $(FFLAGS) mcarlo.cc -c

clean:
	rm *.o mcarlo
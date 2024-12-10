CXX = g++
CXXFLAGS = -O3 -std=c++11 -Wall -Wextra 

SOURCES = main.cc

OBJECTS = $(SOURCES:.cc=.o)

EXEC = main

INCLUDES = -I.

all: $(EXEC)

$(EXEC): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $(EXEC)

%.o: %.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJECTS)
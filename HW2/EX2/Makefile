CXX = g++
CXXFLAGS = -O2 -std=c++11 -Wall -Wextra

SOURCES = SQIsing.cc MCSweeps.cc
OBJECTS = $(SOURCES:.cc=.o)
EXEC = SQIsing
INCLUDES = -I.

all: $(EXEC)

$(EXEC): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(OBJECTS) -o $(EXEC)

%.o: %.cc
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f $(OBJECTS) $(EXEC)
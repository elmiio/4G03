CXX = g++
CXXFLAGS = -O2 -std=c++11 -Wall -Wextra 

SOURCES = numerov.cc func.cc
		  
OBJECTS = $(SOURCES:.cc=.o)
	
EXEC = numerov
	    
INCLUDES = -I.
		   
all: $(EXEC)
	
$(EXEC): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $(EXEC)
								  
%.o: %.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@
													  
clean:
	rm -f $(OBJECTS)
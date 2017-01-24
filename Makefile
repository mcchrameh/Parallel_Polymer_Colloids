CC= mpic++   
CFLAGS=-c -Wall  -std=c++11   -g -Wextra  #-Werror  ### -Wno-write-strings 
LDFLAGS= -lm -pg
SOURCES=CDS_BASE.cpp   simulation.cpp #vtk_export.cpp 
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=Excutable

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(OBJECTS) -o $@  $(LDFLAGS)

.cpp.o:
	$(CC)  $(CFLAGS) $< -o $@  $(LDFLAGS)


clean:	 
	rm -f $(OBJECTS)
	rm -f $(EXECUTABLE)

run:
	mpirun -np 4 ./$(EXECUTABLE) 

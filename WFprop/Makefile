CC              = g++
CC_FLAGS        = -g3 -O3 -Wall -I  /usr/local/include/
LD_FLAGS	= -L/usr/local/lib/ -lgsl -lgslcblas
REC_OBJECTS	= simulation_recombine.o utilities_recombine.o

rec: $(REC_OBJECTS)
	$(CC) $(CC_FLAGS) $(REC_OBJECTS) -o comp_sim_recombine $(LD_FLAGS)
simulation_recombine.o: simulation_recombine.cpp
	$(CC) $(CC_FLAGS) -c simulation_recombine.cpp
utilities_recombine.o: utilities_recombine.cpp
	$(CC) $(CC_FLAGS) -c utilities_recombine.cpp


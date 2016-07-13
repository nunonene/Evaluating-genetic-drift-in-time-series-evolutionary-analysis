CC		= g++
CC_FLAGS	= -g3 -Wall -O3 -I /usr/local/include/
LD_FLAGS	= -L/usr/local/lib/ -lgsl -lgslcblas -lcblas 
filterHDobj	= emission.o common-functions.o drift.o minimization.o MatOp.o


.PHONY: clean all
all: DMS
	rm -f ./*.o

DMS: $(filterHDobj) DMS.o
	$(CC) $(CC_FLAGS) $^ -o DMS $(LD_FLAGS)



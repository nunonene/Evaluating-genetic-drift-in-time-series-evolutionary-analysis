CC              = g++
CC_FLAGS        = -g3 -O3 -Wall -I/usr/local/include -I/usr/local/opt/openblas/include
LD_FLAGS	= -L/usr/local/lib -lgsl -L/usr/local/opt/openblas/lib -lopenblas
REC_OBJECTS	= main.o TMatBinomial.o MatOp.o
rec: $(REC_OBJECTS)
	$(CC) $(CC_FLAGS) $(REC_OBJECTS) -o TransitionMatrixPow_OpenBLAS $(LD_FLAGS)
main.o: main.cpp
	$(CC) $(CC_FLAGS) -c main.cpp
TMatBinomial.o: TMatBinomial.cpp
	$(CC) $(CC_FLAGS) -c TMatBinomial.cpp
MatOp.o: MatOp.cpp
	$(CC) $(CC_FLAGS) -c MatOp.cpp


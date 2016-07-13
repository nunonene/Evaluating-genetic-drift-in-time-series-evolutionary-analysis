CC              = icpc
CC_FLAGS        = -g -O3 -xHOST -Wall -I.
LD_FLAGS	= -lgsl -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm
REC_OBJECTS	= main.o TMatBinomial.o MatOp.o
rec: $(REC_OBJECTS)
	$(CC) $(CC_FLAGS) $(REC_OBJECTS) -o TransitionMatrixPow_MKL $(LD_FLAGS)
main.o: main.cpp
	$(CC) $(CC_FLAGS) -c main.cpp
TMatBinomial.o: TMatBinomial.cpp
	$(CC) $(CC_FLAGS) -c TMatBinomial.cpp
MatOp.o: MatOp.cpp
	$(CC) $(CC_FLAGS) -c MatOp.cpp


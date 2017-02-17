Evaluating-genetic-drift-in-time-series-evolutionary-analysis
=============================================================

Code used in Nene, N.R., Mustonen,V., Illingworth, C.R. Evaluating genetic drift in time-series evolutionary analysis (https://arxiv.org/abs/1611.06152).


InferenceCode folder
-------------------
Forward-backward/predict-update HMM algorithm for drift model parameter estimation (N or sigma) through likelihood maximization and posterior calculation. Executable: DMS .Required external C libraries: gsl (see makefile) 

RunMatPower folder
------------------
Routines for generating matrix powers either using Open_BLAS or MKL (see respective makefiles). Pre-computed matrices for Wright-Fisher propagation are available in InferenceCode/MatPow.


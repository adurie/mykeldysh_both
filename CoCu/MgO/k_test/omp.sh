#!/bin/bash
file_run=$1
icpc $file_run -std=c++14 -qopenmp -O3 -DMKL_ILP64 -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -lpthread -lm -ldl /home/alex/NAG_omp/csl6i26ddl/lib/libnagc_mkl.so
#-xAvx -simd -lifcoremt 
#-g -qopt-report=5 -vec -simd -I~/INTEL/advisor/include

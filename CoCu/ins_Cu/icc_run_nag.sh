#!/bin/bash
file_run=$1
icpc $file_run -std=c++14 -O3 -DMKL_LP64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl -xAvx -simd -lifcoremt /home/alex/NAG/cll6i262dl/lib/libnagc_mkl.so
#-g -qopt-report=5 -vec -simd -I~/INTEL/advisor/include

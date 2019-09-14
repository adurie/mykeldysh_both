#!/bin/bash
file_run=$1
mpiicc $file_run -std=c++14 -O3 -DMKL_LP64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl /home/alex/NAG/cll6i262dl/lib/libnagc_mkl.so

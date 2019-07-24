#!/bin/bash
file_run=$1
icc $file_run -o sc.exe -std=c++14 -O3 -DMKL_LP64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl

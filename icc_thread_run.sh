#!/bin/bash
file_run=$1
icc $file_run -std=c++14 -O3 -DMKL_LP64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread -lm -ldl -fopenmp

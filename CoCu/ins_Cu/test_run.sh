#!/bin/bash
file_run=$1
icpc $file_run -std=c++14 -fast -DMKL_LP64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl -lgsl -o test.exe

#!/bin/bash
file_run=$1
g++ $file_run -std=c++14 -O3 /home/alex/NAG/cll6i262dl/lib/libnagc_nag.so

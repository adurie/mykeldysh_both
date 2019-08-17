#!/bin/bash
file=keldysh_MgO.cpp
for k in {7..1}
do
	let x=k-1
	sed -i "s/double V = 0.0$k/double V = 0.0$x/" $file 
	source icc_run_nag.sh $file 
	mpirun -np 30 ./a.out
done
sed -i "s/double V = 0.00/double V = -0.01/" $file 
source icc_run_nag.sh $file 
mpirun -np 30 ./a.out
for k in {1..4}
do
	let x=k+1
	sed -i "s/double V = -0.0$k/double V = -0.0$x/" $file
	source icc_run_nag.sh $file 
	mpirun -np 30 ./a.out
done

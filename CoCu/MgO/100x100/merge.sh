#! /bin/bash
./convert.exe *'=0.06'* *'=0.05'* *'=0.04'* *'=0.03'* *'=0.02'* *'=0.01'* *'=0.00'* *'=-0.01'* *'=-0.02'* *'=-0.03'* *'=-0.04'* *'=-0.05'*

for i in {2..11}
do
	let k=i+2
	awk "{print \$1, \$$i}" grand_out.txt > $k.txt
done
xmgrace -p appliedv.par 4.txt 5.txt 6.txt 7.txt 8.txt 9.txt 10.txt 11.txt 12.txt 13.txt &

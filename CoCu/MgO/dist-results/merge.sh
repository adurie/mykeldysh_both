#! /bin/bash
for i in {2..11}
do
	let k=i+2
	awk "{print \$1, \$$i}" grand_out.txt > $k.txt
done

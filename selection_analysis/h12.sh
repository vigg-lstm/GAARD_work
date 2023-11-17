#!/bin/bash

pops=(
	Avrankou
	Aboisso
	Baguida
	Korle-Bu
	Madina
	Obuasi
)

for chrom in 2L 2R 3L 3R X
do
	for pop in ${pops[@]}
	do
		python H12.py $pop $chrom 2000 &
	done
	wait
done

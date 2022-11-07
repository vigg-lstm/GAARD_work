#!/bin/bash

pops=(
	Avrankou.coluzzii.Delta
	Baguida.gambiae.Delta
	Baguida.gambiae.PM
	Korle-Bu.coluzzii.Delta
	Korle-Bu.coluzzii.PM
	Madina.gambiae.Delta
	Madina.gambiae.PM
	Obuasi.gambiae.Delta
	Obuasi.gambiae.PM
)

for chrom in 2L #2R 3L 3R X
do
	for pop in ${pops[@]}
	do
		python H12.py $pop $chrom 200 &
	done
	wait
done

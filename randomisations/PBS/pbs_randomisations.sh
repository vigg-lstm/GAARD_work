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

for chrom in 2L 2R 3L 3R X
do
	for pop in ${pops[@]}
	do
		python PBS.py $pop $chrom 200 &
		# We wait a few minutes between runs so that they're not all downloading at once
		sleep 5m
	done
	wait
done

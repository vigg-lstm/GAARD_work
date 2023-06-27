#!/bin/bash
# The difference with download_and_bracken_GAARD_sample.sh is that this version will check 
# whether a fastq file has successfully downloaded and, if not, will try to download it
# again. Also, it will run kraken with memory mapping. 

# define sample id (WA-XXXX)
s=$1
dict_samrun=$2
kraken_db=$3
bracken_db=$4
nt=$5

# directories
dir="."
out_dir=${dir}/$s
mkdir -p $out_dir

# list of runs per s sample  (_1 and _2)
sl=$( grep -w $s $dict_samrun | cut -f14 )

# function: download fastq files from ENA if they don't exist already
# (checks for complete files, incomplete downloads will be restarted)
read_download_ena () {
	echo "### Download and filter fastq.gz files for sample " $s " ###"
	for r in $sl ; do 
		for pair in 1 2 ; do
			if [ ! -f "${out_dir}/${r}_${pair}.fastq.gz" ] ; then
				path1=${r:0:6}
				path2=${r:0:10}
				path3=${r:9:10}
				URL="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${path1}/00${path3}/${path2}/${path2}_${pair}.fastq.gz"
				echo $URL
				# We then put the same test into a while loop. That's because if the file was present before
				# we even tried to download it, then we want to report this (in the else clause below). The 
				# reason for the while loop is that sometimes the download will fail but the script will 
				# happily carry on. If even one of the downloads pass, bracken will still work and we won't
				# realise that something went wrong. 
				download_attempt=1
				while [ ! -f "${out_dir}/${r}_${pair}.fastq.gz" ] ; do
					if [ $download_attempt -gt 1 ] ; then
						echo -e "\t Download attempt $download_attempt"
					fi
					curl $URL -s -S --retry 100 --retry-delay 10 \
					> ${out_dir}/${r}_${pair}.fastq.gz.tmp \
					&& mv ${out_dir}/${r}_${pair}.fastq.gz.tmp ${out_dir}/${r}_${pair}.fastq.gz
					download_attempt=$[$download_attempt+1]
				done
			else echo "${out_dir}/${r}_${pair}.fastq.gz found"
			fi
		done
	done
	echo "### The following fastq files were downloaded ###"
	ls -l ${out_dir}/*.fastq.gz
	echo "### Concatenating all fastq files ###"
	zcat ${out_dir}/*.fastq.gz > ${out_dir}/${s}_joined.fastq
	echo "### Running bracken ###"
	~/software/kraken2-2.0.8-beta/kraken2 --memory-mapping --db ${kraken_db} --threads ${nt} --report ${out_dir}/${s}.kreport ${out_dir}/${s}_joined.fastq > ${out_dir}/${s}.kraken 2> ${out_dir}/${s}_kreport.log
	# We run bracken at both the species and genus levels
	~/software/bracken/Bracken-2.5/bracken -d ${bracken_db} -i ${out_dir}/${s}.kreport -o ${out_dir}/${s}.bracken -r 150 > ${s}_bracken.log
	~/software/bracken/Bracken-2.5/bracken -d ${bracken_db} -i ${out_dir}/${s}.kreport -o ${out_dir}/${s}_genus.bracken -r 150 -l G > ${s}_genus_bracken.log
	echo "### Remove intermediate files ###"
	rm ${out_dir}/*.fastq*
}

read_download_ena

echo "Done"

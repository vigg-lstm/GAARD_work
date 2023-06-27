# This script identifies, downloads and runs bracken on a list of fastq files. 

import time # required to output the time at which the script was run
from sys import stdout # this import is need to flush python output to the stdout (instead of leaving it
# in the buffer
from sys import argv # this import is needed in order for the script to handle command line arguments
import socket # this import allows us to access the name of the computer that the script is being run on
import pandas as pd
import numpy as np
import os
import re
from subprocess import call

if (len(argv) == 5):
	phenotypes_file = argv[1]
	conditions = argv[2]
	path_to_db = argv[3]
	threads = argv[4]
else:
	raise Exception("Fail. There should be four command line arguments (phenotypes_file, conditions, path_to_db, threads).")


print('Running ' + argv[0]  + ' at ' + time.strftime("%H:%M") + ' on ' + time.strftime("%d/%m/%Y") + ' using machine ' + socket.gethostname() + '\n\n')

print('Script was run using the following arguments:\n')
print('\tphenotypes_file = ' + phenotypes_file)
print('\tconditions = ' + conditions)
print('\tpath_to_db = ' + path_to_db)
print('\tthreads = ' + threads + '\n')
stdout.flush()

# Write a function that will execute a string in bash, but return an error in Python if the bash process does
# not return 0
def bash_process(command_string, process_name, shell = False):
	print('\tUsing command ' + command_string)
	stdout.flush()
	if shell:
		result = call(command_string, shell = shell, executable = '/bin/bash')
	else:
		result = os.system(command_string)
	if result != 0:
		raise Exception('bash process ' + process_name + ' failed with error code ' + str(result))

# Use the following accession number table
accession_table_fn = '~/data/ML/GAARD_bracken/data/sample_sequencing_runs_1258specimens.csv'
accession_table = pd.read_csv(accession_table_fn, sep = '\t', names = ['GAARD_SampleID', 'Country', 'Location', 'Plate', 'Well', 'MalGEN_ID', 'Date_of_collection', 'StudyID', 'Type', 'Info1', 'Info2', 'Info3', 'Sample_accession', 'Run_accession'], index_col = 'GAARD_SampleID')
sequenced_samples = accession_table.index

# Load the phenotypes data
phenotypes = pd.read_csv(phenotypes_file, sep = '\t')
sample_names = phenotypes.query(conditions).specimen
print(str(len(sample_names)), 'samples fit the stipulated conditions.')
# Find out how many of these were sequenced
sequenced_sample_names = np.intersect1d(accession_table.index, sample_names)
print(str(len(sequenced_sample_names)), 'of these were sequenced.\n')

# Copy the kraken database to RAM
if not os.path.exists('/dev/shm/kraken_database'):
	bash_process('mkdir /dev/shm/kraken_database', 'Create folder /dev/shm/kraken_database.')
	bash_process('cp ' +  path_to_db + '/*.k2d /dev/shm/kraken_database', 'Copy kraken database in ' + path_to_db + ' to /dev/shm/kraken_database.')
else :
	print('Kraken database already present in /dev/shm')

# Now download the files and run bracken
for s in sequenced_sample_names:
	print('\n')
	print(s)
	if os.path.exists(s):
		print('\tDirectory already exists. Skipping this sample')
		continue
	# Download the fastq files for this sample and run bracken
	download_command = '~/data/ML/GAARD_bracken/download_and_bracken_GAARD_sample3.sh ' + s + ' ' + accession_table_fn + ' /dev/shm/kraken_database ' + path_to_db + ' ' + threads
	bash_process(download_command, 'fastq download and align for sample ' + s)

print('\n\nScript finished running at ' + time.strftime("%H:%M") + ' on ' + time.strftime("%d/%m/%Y") + '\n')


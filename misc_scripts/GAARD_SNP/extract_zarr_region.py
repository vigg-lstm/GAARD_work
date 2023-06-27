# This is a script to extract a given region from a zarr file and output the genotypes, REF and ALT alleles in a csv file
import allel
import zarr
import os
import numpy as np
import pandas as pd
import socket # this import allows us to access the name of the computer that the script is being run on
import time
from sys import stdout, argv
from re import sub

if len(argv) == 6:
	zarr_folder = argv[1]
	study_id = argv[2]
	chrom = argv[3]
	region = argv[4]
	output_fn = argv[5]
	sample_names = 'all'
elif len(argv) == 7:
	zarr_folder = argv[1]
	study_id = argv[2]
	chrom = argv[3]
	region = argv[4]
	output_fn = argv[5]
	sample_names = argv[6]
else:
	raise Exception('Fail. There should be five or six command line arguments (zarr_folder, study_id, chrom, region, output_fn (, sample_names)).')

print('Running ' + argv[0]  + ' at ' + time.strftime("%H:%M") + ' on ' + time.strftime("%d/%m/%Y") + ' using machine ' + socket.gethostname() + '\n\n')

print('Script was run using the following arguments:\n')
print('\tzarr_folder = ' + zarr_folder)
print('\tstudy_id = ' + study_id)
print('\tchrom = ' + chrom)
print('\tregion = ' + region)
print('\toutput_fn = ' + output_fn)
print('\tsample_names = ' + sample_names + '\n')
stdout.flush()

region_start, region_end = [int(x) for x in region.split('-')]

# Load the variant positions
variants_fn = zarr_folder + '/sites/' + chrom + '/variants'
variants = zarr.open(os.path.expanduser(variants_fn), mode = 'r')
pos = allel.SortedIndex(variants['POS'])
focal_indices = pos.locate_range(region_start, region_end)
focal_pos = pos[focal_indices]
alt = variants['ALT'][focal_indices].astype('str')
ref = variants['REF'][focal_indices].astype('str')


# Load the genotype calls 
gt_path = zarr_folder + '/' + study_id + '/' + chrom + '/calldata/GT'
genotypes = allel.GenotypeArray(zarr.open(os.path.expanduser(gt_path), mode = 'r')[focal_indices])
if sample_names != 'all':
	meta_fn = zarr_folder + '/GAARD_preQC_metadata/sample_metadata.csv'
	meta = pd.read_csv(meta_fn, sep = '\t', index_col = 'External_ID')
	sample_names = sample_names.split(',')
	samples_zarr = zarr_folder + '/' + study_id + '/samples'
	# Get the sample names
	malgen_IDs = list(meta.loc[sample_names, 'MalGEN_ID'])                                                      
	# Convert the sample names to Malgen IDs
	malgen_sample_names = [x.decode('utf-8') for x in zarr.open(os.path.expanduser(samples_zarr), mode = 'r')]
	# Get the indices of the vobs sample names that match the malgen_IDs
	sample_boolean = np.isin(np.array([sub('-.*', '', x) for x in malgen_sample_names]), malgen_IDs)
	genotypes = genotypes[:, sample_boolean, :]


allele_counts = genotypes.count_alleles(3)

# Join the output together in a data frame
output_df = pd.concat([pd.DataFrame(focal_pos, columns = ['Pos']), 
                       pd.DataFrame(ref, columns = ['allele0']), 
                       pd.DataFrame(alt, columns = ['allele1','allele2','allele3']), 
                       pd.DataFrame(allele_counts, columns = ['counts0','counts1','counts2','counts3'])], 
                      axis = 1)

output_df.to_csv(output_fn, sep = '\t', index = False)





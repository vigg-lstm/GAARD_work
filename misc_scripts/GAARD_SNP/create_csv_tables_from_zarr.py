import allel
import numpy as np
import zarr
import os
import pandas as pd
import socket # this import allows us to access the name of the computer that the script is being run on
import time
from sys import stdout, argv
from re import sub

if (len(argv) == 4):
	root_folder = argv[1]
	study_id = argv[2]
	conditions = argv[3]
elif (len(argv) == 6):
	root_folder = argv[1]
	study_id = argv[2]
	location = argv[3]
	species = argv[4]
	insecticide = argv[5]
	conditions = 'location == "' + location + '" & species == "' + species + '" & insecticide == "' + insecticide + '"'
else:
	raise Exception("""Fail. There should be three or five command line arguments.
                       With three arguments: (root_folder, study_id, conditions).
                       With five arguments:  (root_folder, study_id, location, species, insecticide).""")

print('Running ' + argv[0]  + ' at ' + time.strftime("%H:%M") + ' on ' + time.strftime("%d/%m/%Y") + ' using machine ' + socket.gethostname() + '\n\n')

print('Script was run using the following arguments:\n')
print('\troot_folder = ' + root_folder)
print('\tstudy_id = ' + study_id)
print('\tconditions = ' + conditions + '\n')
stdout.flush()


# Load the metadata and sample name conversion table
meta_fn = root_folder + '/GAARD_preQC_metadata/sample_metadata.csv'
meta = pd.read_csv(meta_fn, sep = '\t', index_col = 0)

chroms = ('2L', '2R', '3L', '3R', 'X')
zarr_callset = root_folder + '/' + study_id
# Get the sample names
qc_sample_names = [x.decode('utf-8') for x in zarr.open(os.path.expanduser(zarr_callset + '/samples'), mode = 'r')]
# Convert the sample names to GAARD IDs
qc_gaard_sample_names = list(meta.loc[[sub('-.*', '', x) for x in qc_sample_names], 'External_ID'])

# Load our metadata which containes species, location and insecticide information, which we 
# can use to subset the samples. 
phenotypes_fn = root_folder + '/GAARD_preQC_metadata/sample_phenotypes.csv'
phenotypes = pd.read_csv(phenotypes_fn, sep = '\t', index_col = 0)

# Get the samples of interest to us
this_study_gaard_sample_names = phenotypes.query(conditions).index
# Now subset the ones that passed VObs QC
gaard_sample_names = np.array(this_study_gaard_sample_names.intersection(qc_gaard_sample_names))
which_samples = np.isin(qc_gaard_sample_names, gaard_sample_names)

# Load the accessibility map
accessibility_fn = root_folder + '/site_filters/dt_20200416/gamb_colu'
acc = zarr.open(os.path.expanduser(accessibility_fn), mode = 'r')

# Load the variant positions
variants_fn = root_folder + '/sites'
pos = zarr.open(os.path.expanduser(variants_fn), mode = 'r')

# Create functions to filter sites 
def has_mut(x):
	return(np.sum(x <= 0, (1,2)) < np.prod(x.shape[1:3]))

def has_ref(x):
	return(np.sum(x != 0, (1,2)) < np.prod(x.shape[1:3]))

def no_missing(x):
	return(np.sum(x == -1, (1,2)) == 0)

# Create a function to load and filter the SNPs for a given chromosome
def load_and_filter(zarr_path, chrom, verbose = False):
	
	if verbose:
		print('Loading Chromosome ' + chrom)
		stdout.flush()
	
	# Load the variants and positions from the zarr file
	variants_zarr = zarr.open(os.path.expanduser(zarr_path + '/' + chrom + '/calldata/GT'), mode = 'r')
	these_variants = allel.GenotypeArray(variants_zarr.get_orthogonal_selection((np.array(range(variants_zarr.shape[0])), np.where(which_samples)[0], np.array(range(2)))))
	
	# We don't use the allel.is_segregating method because this will also include samples where two 
	# different variants (but no ref) are segregating. Since we will turn the matrix into a presence/
	# absence of the ref, we don't want those sites. 
	hasmut = has_mut(these_variants)
	hasref = has_ref(these_variants)
	nomissing = no_missing(these_variants)
	
	# We try one filter with and one without accessibility filtering
	filter_vector = hasmut & hasref & nomissing 
	filter_vector_acc = filter_vector & np.array(acc[chrom]['variants']['filter_pass'])
	 
	# Get the variant positions
	filtered_variant_pos = pos[chrom]['variants']['POS'][:][filter_vector]
	filtered_variant_pos_acc = pos[chrom]['variants']['POS'][:][filter_vector_acc]
	
	variants_for_output = these_variants[filter_vector, :, :]
	variants_for_output_acc = these_variants[filter_vector_acc, :, :]
	
	# Turn the data into number of mutant alleles at each locus
	mutsum_variants_for_output = variants_for_output.to_n_alt()
	mutsum_variants_for_output_acc = variants_for_output_acc.to_n_alt()
	# Turn the data into presence / absence of a variant in each sample at each locus
	presabs_variants_for_output = np.invert(variants_for_output.is_hom_ref()).astype('int')
	presabs_variants_for_output_acc = np.invert(variants_for_output_acc.is_hom_ref()).astype('int')
	
	# Create an LD pruning filter. 
	LD_prune = allel.locate_unlinked(mutsum_variants_for_output)
	LD_prune_acc = allel.locate_unlinked(mutsum_variants_for_output_acc)
	
	# Unpruned mutsum tables
	mutsum_df = pd.concat([pd.DataFrame([chrom]*len(filtered_variant_pos), columns = ['Chrom']), pd.DataFrame(filtered_variant_pos, columns = ['Pos']), pd.DataFrame(mutsum_variants_for_output, columns = gaard_sample_names)], 1)
	mutsum_df_acc = pd.concat([pd.DataFrame([chrom]*len(filtered_variant_pos_acc), columns = ['Chrom']), pd.DataFrame(filtered_variant_pos_acc, columns = ['Pos']), pd.DataFrame(mutsum_variants_for_output_acc, columns = gaard_sample_names)], 1)
	# Pruned mutsum tables
	mutsum_df_pruned = pd.concat([pd.DataFrame([chrom]*np.sum(LD_prune), columns = ['Chrom']), pd.DataFrame(filtered_variant_pos.compress(LD_prune), columns = ['Pos']), pd.DataFrame(mutsum_variants_for_output.compress(LD_prune, axis = 0), columns = gaard_sample_names)], 1)
	mutsum_df_acc_pruned = pd.concat([pd.DataFrame([chrom]*np.sum(LD_prune_acc), columns = ['Chrom']), pd.DataFrame(filtered_variant_pos_acc.compress(LD_prune_acc), columns = ['Pos']), pd.DataFrame(mutsum_variants_for_output_acc.compress(LD_prune_acc, axis = 0), columns = gaard_sample_names)], 1)
	# Unpruned presence / absence tables
	presabs_df = pd.concat([pd.DataFrame([chrom]*len(filtered_variant_pos), columns = ['Chrom']), pd.DataFrame(filtered_variant_pos, columns = ['Pos']), pd.DataFrame(presabs_variants_for_output, columns = gaard_sample_names)], 1)
	presabs_df_acc = pd.concat([pd.DataFrame([chrom]*len(filtered_variant_pos_acc), columns = ['Chrom']), pd.DataFrame(filtered_variant_pos_acc, columns = ['Pos']), pd.DataFrame(presabs_variants_for_output_acc, columns = gaard_sample_names)], 1)
	# Pruned presence / absence tables
	presabs_df_pruned = pd.concat([pd.DataFrame([chrom]*np.sum(LD_prune), columns = ['Chrom']), pd.DataFrame(filtered_variant_pos.compress(LD_prune), columns = ['Pos']), pd.DataFrame(presabs_variants_for_output.compress(LD_prune, axis = 0), columns = gaard_sample_names)], 1)
	presabs_df_acc_pruned = pd.concat([pd.DataFrame([chrom]*np.sum(LD_prune_acc), columns = ['Chrom']), pd.DataFrame(filtered_variant_pos_acc.compress(LD_prune_acc), columns = ['Pos']), pd.DataFrame(presabs_variants_for_output_acc.compress(LD_prune_acc, axis = 0), columns = gaard_sample_names)], 1)
	
	return([mutsum_df, mutsum_df_acc, mutsum_df_pruned, mutsum_df_acc_pruned, presabs_df, presabs_df_acc, presabs_df_pruned, presabs_df_acc_pruned])


all_variants_list = [load_and_filter(zarr_callset, chrom, True) for chrom in chroms]

mutsum_df = pd.concat([x[0] for x in all_variants_list], 0)
mutsum_df_acc = pd.concat([x[1] for x in all_variants_list], 0)
mutsum_df_pruned = pd.concat([x[2] for x in all_variants_list], 0)
mutsum_df_acc_pruned = pd.concat([x[3] for x in all_variants_list], 0)
presabs_df = pd.concat([x[4] for x in all_variants_list], 0)
presabs_df_acc = pd.concat([x[5] for x in all_variants_list], 0)
presabs_df_pruned = pd.concat([x[6] for x in all_variants_list], 0)
presabs_df_acc_pruned = pd.concat([x[7] for x in all_variants_list], 0)

mutsum_df.to_csv('filtered_SNPs.csv', sep = '\t', index = False)
mutsum_df_acc.to_csv('filtered_SNPs_acc.csv', sep = '\t', index = False)
mutsum_df_pruned.to_csv('filtered_SNPs_pruned.csv', sep = '\t', index = False)
mutsum_df_acc_pruned.to_csv('filtered_SNPs_acc_pruned.csv', sep = '\t', index = False)
presabs_df.to_csv('presabs_filtered_SNPs.csv', sep = '\t', index = False)
presabs_df_acc.to_csv('presabs_filtered_SNPs_acc.csv', sep = '\t', index = False)
presabs_df_pruned.to_csv('presabs_filtered_SNPs_pruned.csv', sep = '\t', index = False)
presabs_df_acc_pruned.to_csv('presabs_filtered_SNPs_acc_pruned.csv', sep = '\t', index = False)



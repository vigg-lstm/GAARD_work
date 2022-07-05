import time # required to output the time at which the script was run
from sys import stdout # this import is need to flush python output to the stdout (instead of leaving it
# in the buffer
from sys import argv # this import is needed in order for the script to handle command line arguments
import socket # this import allows us to access the name of the computer that the script is being run on
import re
from os import system
import malariagen_data
import allel
import scipy.cluster
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import seaborn as sns

if (len(argv) == 4):
	study_pop = argv[1]
	region = argv[2]
	repo_folder = argv[3]
	output_filename = f'{study_pop}_{region}.csv'
elif (len(argv) == 5):
	study_pop = argv[1]
	region = argv[2]
	repo_folder = argv[3]
	output_filename = argv[4]
else:
	raise Exception("Fail. There should be three or four command line arguments (study_pop, region, repo_folder (, output_filename)).")


print('Running ' + argv[0]  + ' at ' + time.strftime("%H:%M") + ' on ' + time.strftime("%d/%m/%Y") + ' using machine ' + socket.gethostname() + '\n\n')

print('Script was run using the following arguments:\n')
print('\tstudy_pop = ' + study_pop)
print('\tregion = ' + region)
print('\trepo_folder = ' + repo_folder)
print('\toutput_filename = ' + output_filename + '\n')
stdout.flush()

ag3 = malariagen_data.Ag3('gs://vo_agam_release', pre = True)
ag3.sample_sets(release = '3.2')
ag3

# Get the metadata and phenotypes, and drop the males
meta = pd.merge(
    pd.read_csv(
        repo_folder + '/data/combined/all_samples.samples.meta.csv', 
        sep = ',', 
        index_col = 'partner_sample_id'
    ).drop(['location', 'country'], axis = 1),
    pd.read_csv(
        repo_folder + '/data/combined/sample_phenotypes.csv', 
        sep = '\t', 
        index_col = 'specimen'
    ),
    left_index = True, 
    right_index = True
).query('sex_call == "F"')

# Get sib groups and kick out sibs 
sib_groups = pd.read_csv(repo_folder + '/NGSrelate/full_relatedness/sib_group_table.csv', sep = '\t', index_col = 'sample.name')
kick_out = sib_groups.query('keep == False').index
meta.drop(kick_out, inplace = True)

# Subset to only the phenotypes that we are interested in
meta['population'] = meta[['location', 'species', 'insecticide']].apply('_'.join, 1)
meta.query(f'population == "{study_pop}"', inplace = True)

# Get the haplotypes from the region of interest
haps = ag3.haplotypes(
    region = region,
    analysis = 'gamb_colu', 
    sample_sets = '3.2',
    sample_query = f'sample_id in {list(meta.sample_id)}',
)

haparray = allel.GenotypeArray(haps['call_genotype']).to_haplotypes()
hapsamples = np.array(np.repeat(haps['sample_id'], 2))
dist = allel.pairwise_distance(haparray, metric = 'hamming')

site_filter = ag3.snp_calls(region=region, sample_sets="3.2")['variant_filter_pass_gamb_colu']
n_bases = np.sum(site_filter.values)
dist_dxy = dist * haparray.n_variants / n_bases

def find_clusters(dist, n, threshold=0.001, method='complete'):
        # build hierarchy
        clust = scipy.cluster.hierarchy.linkage(dist, method=method)
        # find clusters
        f = scipy.cluster.hierarchy.fcluster(clust, threshold,
                                             criterion='distance')
        # compute cluster sizes
        fsz = np.bincount(f)
        # sort largest first
        fsort = np.argsort(fsz)[::-1]
        # take largest n
        fsort = fsort[:n]
        # get haplotype indices for each cluster
        clusters = [set(np.nonzero(f == i)[0]) for i in fsort]
        return clusters

focal_clusters = find_clusters(dist_dxy, n=50, threshold=0.001)

min_cluster_size = 20
large_focal_clusters = [cluster for cluster in focal_clusters if len(cluster) > min_cluster_size]
print(f'There were {len(large_focal_clusters)} clusters larger than {min_cluster_size} haplotypes.')
stdout.flush()

def assign_sample_haplotypes(cluster_ids):
    cluster_samples = hapsamples[list(cluster_ids)]
    sample_hap_genotype = meta.sample_id.apply(lambda x: np.sum(np.isin(cluster_samples, x)))
    return(sample_hap_genotype)

for i in range(len(large_focal_clusters)):
    meta[f'cluster_{i+1}'] = assign_sample_haplotypes(large_focal_clusters[i])

output_table = meta[['phenotype'] + [f'cluster_{i+1}' for i in range(len(large_focal_clusters))]]

output_table.to_csv(output_filename, sep = '\t', index_label = 'sample_name')

# Now we want to work out the sequence of each haplotype
window_alleles = pd.concat((
    pd.DataFrame({'Chrom': [haps.contigs[x] for x in haps['variant_contig'].values]}),
    haps['variant_position'].to_pandas().rename('Pos'),
    haps['variant_allele'].astype('str').to_pandas().set_axis(['Ref', 'Alt'], axis = 1)),
    axis = 1
)

remove_rows = np.full(window_alleles.shape[0], True)
for i in range(len(large_focal_clusters)):
	cluster_haparray = haparray[:, list(large_focal_clusters[i])]
	# We want a table that shows, for each position, the ref, the alt, and the number of
	# each in the haplotypes
	n_alt = np.sum(cluster_haparray, axis = 1)
	window_alleles[f'n_Ref_cluster_{i+1}'] = cluster_haparray.shape[1] - n_alt
	window_alleles[f'n_Alt_cluster_{i+1}'] = n_alt
	remove_rows = remove_rows & (n_alt == 0)

non_cluster = set(range(haparray.shape[1])).difference(*large_focal_clusters)
non_cluster_haparray = haparray[:, list(non_cluster)]
n_alt_non_cluster = np.sum(non_cluster_haparray, axis = 1)
window_alleles[f'n_Ref_non_cluster'] = non_cluster_haparray.shape[1] - n_alt_non_cluster
window_alleles[f'n_Alt_non_cluster'] = n_alt_non_cluster

# We want to remove the rows where all clusters are fully Ref. 
window_alleles = window_alleles.loc[~remove_rows, :]

haparray_output_filename = re.sub('\.[^.]+$', '_hap_sequence.csv', output_filename)
window_alleles.to_csv(haparray_output_filename, sep = '\t', index = False)

window_alleles['ID'] = '.'
window_alleles['QUAL'] = '.'
window_alleles['FILTER'] = '.'
window_alleles['INFO'] = '.'

vcf_table = window_alleles.rename(
	columns = {
		'Chrom': '#CHROM',
		'Pos': 'POS',
		'Ref': 'REF',
		'Alt': 'ALT'
	}
)[['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']]

vcf_filename = re.sub('\.[^.]+$', '_temp.vcf', output_filename)
snpeff_filename = re.sub('\.[^.]+$', '.vcf', output_filename)
vcf_table.to_csv(vcf_filename, sep = '\t', index = False)

# Run SNPEff
snp_eff_command = ('java -jar ~/software/snpEff/snpEff.jar eff Anopheles_gambiae -noStats ' +
	               vcf_filename + 
                   ' > ' +
                   snpeff_filename)

print('\nRunning snpEff on haplotype SNPs, using command ' + snp_eff_command + '\n')
stdout.flush()
system(snp_eff_command)
system('rm ' + vcf_filename)

# Now make a dendrogram plot
def truspan(cluster, r):
	# get the index of the cluster haps in the dendrogram list of all haps
	cluster_leaves = sorted([r['leaves'].index(i) for i in cluster])
	# are these indices monotonic - they should be!
	x = np.asarray(cluster_leaves)
	dx = np.diff(x)
	mon = np.all(dx == 1)
	assert mon
	return min(cluster_leaves), max(cluster_leaves)


def plot_dendrogram(zhier, ax, h, method='complete', color_threshold=0, above_threshold_color='k'):
	
	# plot dendrogram
	sns.despine(ax=ax, offset=5, bottom=True, top=False)
	r = scipy.cluster.hierarchy.dendrogram(zhier, no_labels=True, count_sort=True, 
										   color_threshold=color_threshold, 
										   above_threshold_color=above_threshold_color,
										   ax=ax)
	xmin, xmax = ax.xaxis.get_data_interval()
	xticklabels = np.array(list(range(0, h.shape[1], 200)) + [h.shape[1]])
	xticks = xticklabels / h.shape[1]
	xticks = (xticks * (xmax - xmin)) + xmin
	ax.set_xticks(xticks)
	ax.set_xticklabels(xticklabels)
	#ax.set_xlabel('Haplotypes')
	ax.xaxis.set_label_position('top')
	ax.set_ylim(bottom=-.0001)
	
	ax.set_ylabel(r'$d_{xy}$')
	ax.autoscale(axis='x', tight=True)
	return z, r

phen_colours = {'dead': 'orangered', 'alive': 'deepskyblue'}

def draw_hap_cluster_plot(z, r, h, cluster_labels, vspans, add_legend = True):
	
	gs = GridSpec(3, 1, height_ratios=[5.0, 0.5, 0.5])
	fig = plt.figure(figsize=(15, 4))
	
	ax1 = plt.subplot(gs[0])
	sns.despine(ax=ax1, offset=5, bottom=True, top=False)
	_ = plot_dendrogram(z, ax1, h)
	ax1.spines['top'].set_visible(False)
	ax1.set_xticks([])
	
	ax_pops = fig.add_subplot(gs[1])
	
	hap_phen = meta.set_index('sample_id').loc[hapsamples, 'phenotype'].values
	
	x = hap_phen.take(r['leaves'])
	hap_clrs = [phen_colours[p] for p in x]
	ax_pops.broken_barh(xranges=[(i, 1) for i in range(h.shape[1])], 
						yrange=(0, 1),
						color=hap_clrs);
	sns.despine(ax=ax_pops, offset=5, left=True, bottom=True)
	
	ax_pops.set_xticks([])
	ax_pops.set_yticks([])
	ax_pops.set_xlim(0, h.shape[1])
	ax_pops.yaxis.set_label_position('left')
	ax_pops.set_ylabel('Population', rotation=0, ha='right', va='center')
	
	if add_legend:
		unique_phenotypes = np.unique(x)
		plot_x_range = ax1.get_xlim()[1] - ax1.get_xlim()[0]
		plot_y_range = ax1.get_ylim()[1] - ax1.get_ylim()[0]
		legend_x = ax1.get_xlim()[0] + plot_x_range * 0.95
		for i, k in enumerate(unique_phenotypes[::-1]):
			legend_y = ax1.get_ylim()[0] + plot_y_range * (0.6 + 0.1 * i)
			ax1.add_patch(mpl.patches.Rectangle((legend_x,legend_y), plot_x_range/50,plot_y_range/15, color = phen_colours[k]))
			ax1.text(legend_x*1.03,legend_y + plot_y_range / 50, k)
	
	# cluster brackets
	ax_clu = fig.add_subplot(gs[2])
	sns.despine(ax=ax_clu, bottom=True, left=True)
	ax_clu.set_xlim(0, h.shape[1])
	ax_clu.set_ylim(0, 1)
	for lbl, (xmin, xmax) in zip(cluster_labels, vspans):
		if lbl:
			# hack to get the "fraction" right, which controls length of bracket arms
			fraction = -2 / (xmax - xmin)
			ax_clu.annotate("", ha='left', va='center',
							xy=(xmin, 1), xycoords='data',
							xytext=(xmax, 1), textcoords='data',
							arrowprops=dict(arrowstyle="-",
											connectionstyle="bar,fraction=%.4f" % fraction,
											),
							)
			ax_clu.text((xmax + xmin)/2, 0.1, lbl, va='top', ha='center')
	
	ax_clu.set_xticks([])
	ax_clu.set_yticks([])
	ax1.set_title(f'{study_pop}_{region}')



z = scipy.cluster.hierarchy.linkage(dist_dxy, method="complete")
r = scipy.cluster.hierarchy.dendrogram(
        z, no_labels=True, count_sort=True, 
        color_threshold=0, no_plot=True,
        above_threshold_color='k')
v_span = [truspan(cluster, r) for cluster in large_focal_clusters]
cluster_lab = ['Cluster' + str(i+1) for i in range(len(large_focal_clusters))]

dendrogram_filename = re.sub('\.[^.]+$', '_dendrogram.png', output_filename)
draw_hap_cluster_plot(z, r, haparray, cluster_labels=cluster_lab, vspans=v_span)
plt.savefig(dendrogram_filename, format = 'png', dpi = 100)



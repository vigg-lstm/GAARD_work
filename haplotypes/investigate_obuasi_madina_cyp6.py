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

study_pops = ['Madina_gambiae_Delta', 'Obuasi_gambiae_Delta']
region1 = '2R:28454838-28474602'
region2 = '2R:28475225-28496115'
min_cluster_size = 20

print('Running ' + argv[0]  + ' at ' + time.strftime("%H:%M") + ' on ' + time.strftime("%d/%m/%Y") + ' using machine ' + socket.gethostname() + '\n\n')

ag3 = malariagen_data.Ag3('gs://vo_agam_release', pre = True)
ag3.sample_sets(release = '3.2')
ag3

# Get the metadata and phenotypes, and drop the males
meta = pd.merge(
    pd.read_csv(
        '../data/combined/all_samples.samples.meta.csv',
        sep = ',',
        index_col = 'partner_sample_id'
    ).drop(['location', 'country'], axis = 1),
    pd.read_csv(
        '../data/combined/sample_phenotypes.csv',
        sep = '\t',
        index_col = 'specimen'
    ),
    left_index = True,
    right_index = True
).query('sex_call == "F"')

# Get sib groups and kick out sibs
sib_groups = pd.read_csv('../NGSrelate/full_relatedness/sib_group_table.csv', sep = '\t', index_col = 'sample.name')
kick_out = sib_groups.query('keep == False').index
meta.drop(kick_out, inplace = True)

# Subset to only the phenotypes that we are interested in
meta['population'] = meta[['location', 'species', 'insecticide']].apply('_'.join, 1)
meta.query(f'population in {study_pops}', inplace = True)

# Write a function to identify haplotyde clusters from a distance matrix.
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

# Some plotting functions
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
	return r

loc_colours = {'Madina': 'purple', 'Obuasi': 'green'}

def draw_hap_cluster_plot(z, r, h, cluster_labels, vspans, colour_scheme, labels_for_colours, add_legend = True, title = 'dendrogram'):

	gs = GridSpec(3, 1, height_ratios=[5.0, 0.5, 0.5])
	fig = plt.figure(figsize=(15, 4))

	ax1 = plt.subplot(gs[0])
	sns.despine(ax=ax1, offset=5, bottom=True, top=False)
	_ = plot_dendrogram(z, ax1, h)
	ax1.spines['top'].set_visible(False)
	ax1.set_xticks([])

	ax_pops = fig.add_subplot(gs[1])

	x = labels_for_colours.take(r['leaves'])
	hap_clrs = [colour_scheme[p] for p in x]
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
		unique_label = np.unique(x)
		plot_x_range = ax1.get_xlim()[1] - ax1.get_xlim()[0]
		plot_y_range = ax1.get_ylim()[1] - ax1.get_ylim()[0]
		legend_x = ax1.get_xlim()[0] + plot_x_range * 0.95
		for i, k in enumerate(unique_label[::-1]):
			legend_y = ax1.get_ylim()[0] + plot_y_range * (0.6 + 0.1 * i)
			ax1.add_patch(mpl.patches.Rectangle((legend_x,legend_y), plot_x_range/50,plot_y_range/15, color = colour_scheme[k]))
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
	ax1.set_title(f'{title}')

# Get the haplotypes from the regions of interest
haps1 = ag3.haplotypes(
    region = region1,
    analysis = 'gamb_colu',
    sample_sets = '3.2',
    sample_query = f'sample_id in {list(meta.sample_id)}',
)

haparray1 = allel.GenotypeArray(haps1['call_genotype']).to_haplotypes()
hapsamples1 = np.array(np.repeat(haps1['sample_id'], 2))
hapnames1 = np.concatenate([[x + 'a', x + 'b'] for x in haps1['sample_id'].values])

dist1 = allel.pairwise_distance(haparray1, metric = 'hamming')

site_filter1 = ag3.snp_calls(region=region1, sample_sets="3.2")['variant_filter_pass_gamb_colu']
n_bases1 = np.sum(site_filter1.values)
dist_dxy1 = dist1 * haparray1.n_variants / n_bases1

focal_clusters1 = find_clusters(dist_dxy1, n=50, threshold=0.001)

large_focal_clusters1 = [cluster for cluster in focal_clusters1 if len(cluster) > min_cluster_size]


z1 = scipy.cluster.hierarchy.linkage(dist_dxy1, method="complete")
r1 = scipy.cluster.hierarchy.dendrogram(
        z1, no_labels=True, count_sort=True,
        color_threshold=0, no_plot=True,
        above_threshold_color='k')
v_span1 = [truspan(cluster, r1) for cluster in large_focal_clusters1]
cluster_lab1 = ['Cluster' + str(i+1) for i in range(len(large_focal_clusters1))]

hap_labels = meta.set_index('sample_id').loc[hapsamples1, 'location'].values
draw_hap_cluster_plot(z1, r1, haparray1, cluster_labels=cluster_lab1,
                      vspans=v_span1, colour_scheme = loc_colours,
                      labels_for_colours = hap_labels, title = 'Window1'
)

# So the biggest cluster is shared between the two populations. Is it the same as the
# largest clusters identified in each pop respecitvely?
meta_madina = meta.query(f'population == "Madina_gambiae_Delta"')
haps1_madina = ag3.haplotypes(
    region = region1,
    analysis = 'gamb_colu',
    sample_sets = '3.2',
    sample_query = f'sample_id in {list(meta_madina.sample_id)}',
)
haparray1_madina = allel.GenotypeArray(haps1_madina['call_genotype']).to_haplotypes()
hapnames1_madina = np.concatenate([[x + 'a', x + 'b'] for x in haps1_madina['sample_id'].values])
dist1_madina = allel.pairwise_distance(haparray1_madina, metric = 'hamming')
dist_dxy1_madina = dist1_madina * haparray1_madina.n_variants / n_bases1
focal_cluster1_madina = find_clusters(dist_dxy1_madina, n=1, threshold=0.001)[0]
focal_cluster1_haps_madina = hapnames1_madina[list(focal_cluster1_madina)]

meta_obuasi = meta.query(f'population == "Obuasi_gambiae_Delta"')
haps1_obuasi = ag3.haplotypes(
    region = region1,
    analysis = 'gamb_colu',
    sample_sets = '3.2',
    sample_query = f'sample_id in {list(meta_obuasi.sample_id)}',
)
haparray1_obuasi = allel.GenotypeArray(haps1_obuasi['call_genotype']).to_haplotypes()
hapnames1_obuasi = np.concatenate([[x + 'a', x + 'b'] for x in haps1_obuasi['sample_id'].values])
dist1_obuasi = allel.pairwise_distance(haparray1_obuasi, metric = 'hamming')
dist_dxy1_obuasi = dist1_obuasi * haparray1_obuasi.n_variants / n_bases1
focal_cluster1_obuasi = find_clusters(dist_dxy1_obuasi, n=1, threshold=0.001)[0]
focal_cluster1_haps_obuasi = hapnames1_obuasi[list(focal_cluster1_obuasi)]

hap_assignment1 = np.full(len(hapnames1), '', dtype = object)
hap_assignment1[np.isin(hapnames1, hapnames1_madina)] = 'Madina_noncluster'
hap_assignment1[np.isin(hapnames1, hapnames1_obuasi)] = 'Obuasi_noncluster'
hap_assignment1[np.isin(hapnames1, focal_cluster1_haps_madina)] = 'Madina_cluster'
hap_assignment1[np.isin(hapnames1, focal_cluster1_haps_obuasi)] = 'Obuasi_cluster'
cluster_colours = {'Madina_cluster': 'purple',
                   'Obuasi_cluster': 'green',
                   'Madina_noncluster': 'orchid',
                   'Obuasi_noncluster': 'springgreen'}
draw_hap_cluster_plot(z1, r1, haparray1, cluster_labels=cluster_lab1,
                      vspans=v_span1, colour_scheme = cluster_colours,
                      labels_for_colours = hap_assignment1, title = 'Window1'
)

# Indeed, the two clusters are identical
hap_clusters_split = set(np.concatenate([focal_cluster1_haps_obuasi, focal_cluster1_haps_madina]))
hap_clusters_together = set(hapnames1[list(large_focal_clusters1[0])])
if (hap_clusters_split == hap_clusters_together):
	print('In window1, the main clusters was identical whether Obuasi and Madina were analysed together or separately.')

# Now lets do the same with the clusters that's over Cyp6aa1 in Obuasi
haps2 = ag3.haplotypes(
    region = region2,
    analysis = 'gamb_colu',
    sample_sets = '3.2',
    sample_query = f'sample_id in {list(meta.sample_id)}',
)

haparray2 = allel.GenotypeArray(haps2['call_genotype']).to_haplotypes()
hapsamples2 = np.array(np.repeat(haps2['sample_id'], 2))
hapnames2 = np.concatenate([[x + 'a', x + 'b'] for x in haps2['sample_id'].values])

dist2 = allel.pairwise_distance(haparray2, metric = 'hamming')

site_filter2 = ag3.snp_calls(region=region2, sample_sets="3.2")['variant_filter_pass_gamb_colu']
n_bases2 = np.sum(site_filter2.values)
dist_dxy2 = dist2 * haparray2.n_variants / n_bases2

focal_clusters2 = find_clusters(dist_dxy2, n=50, threshold=0.001)

large_focal_clusters2 = [cluster for cluster in focal_clusters2 if len(cluster) > min_cluster_size]
focal_cluster2_haps = hapnames2[list(large_focal_clusters2[0])]

z2 = scipy.cluster.hierarchy.linkage(dist_dxy2, method="complete")
r2 = scipy.cluster.hierarchy.dendrogram(
        z2, no_labels=True, count_sort=True,
        color_threshold=0, no_plot=True,
        above_threshold_color='k')
v_span2 = [truspan(cluster, r2) for cluster in large_focal_clusters2]
cluster_lab2 = ['Cluster' + str(i+1) for i in range(len(large_focal_clusters2))]

hap_labels = meta.set_index('sample_id').loc[hapsamples2, 'location'].values
draw_hap_cluster_plot(z2, r2, haparray2, cluster_labels=cluster_lab2,
                      vspans=v_span2, colour_scheme = loc_colours,
                      labels_for_colours = hap_labels, title = 'Window2'
)
plt.show()

# So the biggest cluster is shared between the two populations. Is it the same as the
# largest clusters identified in each pop respecitvely?
haps2_madina = ag3.haplotypes(
    region = region2,
    analysis = 'gamb_colu',
    sample_sets = '3.2',
    sample_query = f'sample_id in {list(meta_madina.sample_id)}',
)
haparray2_madina = allel.GenotypeArray(haps2_madina['call_genotype']).to_haplotypes()
hapnames2_madina = np.concatenate([[x + 'a', x + 'b'] for x in haps2_madina['sample_id'].values])
dist2_madina = allel.pairwise_distance(haparray2_madina, metric = 'hamming')
dist_dxy2_madina = dist2_madina * haparray2_madina.n_variants / n_bases2
focal_cluster2_madina = find_clusters(dist_dxy2_madina, n=1, threshold=0.001)[0]
focal_cluster2_haps_madina = hapnames2_madina[list(focal_cluster2_madina)]

haps2_obuasi = ag3.haplotypes(
    region = region2,
    analysis = 'gamb_colu',
    sample_sets = '3.2',
    sample_query = f'sample_id in {list(meta_obuasi.sample_id)}',
)
haparray2_obuasi = allel.GenotypeArray(haps2_obuasi['call_genotype']).to_haplotypes()
hapnames2_obuasi = np.concatenate([[x + 'a', x + 'b'] for x in haps2_obuasi['sample_id'].values])
dist2_obuasi = allel.pairwise_distance(haparray2_obuasi, metric = 'hamming')
dist_dxy2_obuasi = dist2_obuasi * haparray2_obuasi.n_variants / n_bases2
focal_cluster2_obuasi = find_clusters(dist_dxy2_obuasi, n=1, threshold=0.001)[0]
focal_cluster2_haps_obuasi = hapnames2_obuasi[list(focal_cluster2_obuasi)]


hap_assignment2 = np.full(len(hapnames2), '', dtype = object)
hap_assignment2[np.isin(hapnames2, hapnames2_madina)] = 'Madina_noncluster'
hap_assignment2[np.isin(hapnames2, hapnames2_obuasi)] = 'Obuasi_noncluster'
hap_assignment2[np.isin(hapnames2, focal_cluster2_haps_madina)] = 'Madina_cluster'
hap_assignment2[np.isin(hapnames2, focal_cluster2_haps_obuasi)] = 'Obuasi_cluster'
cluster_colours = {'Madina_cluster': 'purple',
                   'Obuasi_cluster': 'green',
                   'Madina_noncluster': 'orchid',
                   'Obuasi_noncluster': 'springgreen'}
draw_hap_cluster_plot(z2, r2, haparray2, cluster_labels=cluster_lab2,
                      vspans=v_span2, colour_scheme = cluster_colours,
                      labels_for_colours = hap_assignment2, title = ''
)
dendrogram_filename = 'Madina_Obuasi_delta_shared_Cyp6_cluster_dendrogram.png'
plt.savefig(dendrogram_filename, format = 'png', dpi = 100)

# Now colour-code by alive/dead instead.
hap_phen2 = meta.set_index('sample_id').loc[hapsamples2, 'phenotype'].values
hap_assignment2_deadalive = np.full(len(hapnames2), '', dtype = object)
hap_assignment2_deadalive[np.isin(hapnames2, hapnames2_madina) & (hap_phen2 == 'alive')] = 'Madina_alive'
hap_assignment2_deadalive[np.isin(hapnames2, hapnames2_obuasi) & (hap_phen2 == 'alive')] = 'Obuasi_alive'
hap_assignment2_deadalive[np.isin(hapnames2, hapnames2_madina) & (hap_phen2 == 'dead')] = 'Madina_dead'
hap_assignment2_deadalive[np.isin(hapnames2, hapnames2_obuasi) & (hap_phen2 == 'dead')] = 'Obuasi_dead'
phen_loc_colours = {'Madina_alive': 'deepskyblue',
                    'Obuasi_alive': 'lightskyblue',
                    'Madina_dead': 'orangered',
                    'Obuasi_dead': 'coral'}
draw_hap_cluster_plot(z2, r2, haparray2, cluster_labels=cluster_lab2,
                      vspans=v_span2, colour_scheme = phen_loc_colours,
                      labels_for_colours = hap_assignment2_deadalive,
                      title = ''
)
dendrogram_filename_phen = 'Madina_Obuasi_delta_shared_Cyp6_cluster_dendrogram_phenlabel.png'
plt.savefig(dendrogram_filename, format = 'png', dpi = 100)

# Indeed, the two clusters are identical
hap_clusters_split = set(np.concatenate([focal_cluster2_haps_obuasi, focal_cluster2_haps_madina]))
hap_clusters_together = set(focal_cluster2_haps)
if (hap_clusters_split == hap_clusters_together):
	print('In window2, the main clusters was identical whether Obuasi and Madina were analysed together or separately.')

# Since we have the same sweep in Obuasi and Madina, it makes sense to combine the association
# analysis for that haplotype.

def assign_sample_haplotypes(hapsamples, cluster_ids):
    cluster_samples = hapsamples[list(cluster_ids)]
    sample_hap_genotype = meta.sample_id.apply(lambda x: np.sum(np.isin(cluster_samples, x)))
    return(sample_hap_genotype)

meta['Cyp6_main_cluster'] = assign_sample_haplotypes(hapsamples2, large_focal_clusters2[0])

output_table = meta[['location', 'insecticide', 'phenotype', 'Cyp6_main_cluster']]

output_filename = 'Madina_Obuasi_delta_shared_Cyp6_cluster.csv'
output_table.to_csv(output_filename, sep = '\t', index_label = 'sample_name')

# Now we want to work out the sequence of each haplotype
window_alleles = pd.concat((
    pd.DataFrame({'Chrom': [haps2.contigs[x] for x in haps2['variant_contig'].values]}),
    haps2['variant_position'].to_pandas().rename('Pos'),
    haps2['variant_allele'].astype('str').to_pandas().set_axis(['Ref', 'Alt'], axis = 1)),
    axis = 1
)

cluster_haparray = haparray2[:, list(large_focal_clusters2[0])]
# We want a table that shows, for each position, the ref, the alt, and the number of each
# in the haplotypes
n_alt = np.sum(cluster_haparray, axis = 1)
window_alleles['n_Ref'] = cluster_haparray.shape[1] - n_alt
window_alleles['n_Alt'] = n_alt

non_cluster = set(range(haparray2.shape[1])).difference(large_focal_clusters2[0])
non_cluster_haparray = haparray2[:, list(non_cluster)]
n_alt_non_cluster = np.sum(non_cluster_haparray, axis = 1)
window_alleles[f'n_Ref_non_cluster'] = non_cluster_haparray.shape[1] - n_alt_non_cluster
window_alleles[f'n_Alt_non_cluster'] = n_alt_non_cluster

window_alleles

# We want to remove the rows where all clusters are fully Ref.
window_alleles = window_alleles.loc[n_alt > 0, :]
window_alleles

haparray_output_filename = 'Madina_Obuasi_delta_shared_Cyp6_cluster_hap_sequence.csv'
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

vcf_filename = 'Madina_Obuasi_delta_shared_Cyp6_cluster_temp.vcf'
snpeff_filename = 'Madina_Obuasi_delta_shared_Cyp6_cluster.vcf'
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

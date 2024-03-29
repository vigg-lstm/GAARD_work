import os
import re
import numpy as np
import pandas as pd

populations = ['Avrankou_coluzzii_Delta',
               'Baguida_gambiae_Delta',
               'Korle-Bu_coluzzii_Delta',
               'Madina_gambiae_Delta',
               'Obuasi_gambiae_Delta',
               'Baguida_gambiae_PM',
               'Korle-Bu_coluzzii_PM',
               'Madina_gambiae_PM',
               'Obuasi_gambiae_PM']

populations_italics = {p: re.sub('_(.+)_', r'\_*\1*\_', p) for p in populations}

all_files = os.listdir('../../haplotypes')

all_windows = np.unique([re.findall('^[^.]+(?=.vcf)', x)[0] for x in all_files if '.vcf' in x])

output_file = 'window_Fst_haplotype_summary.md'

f = open(output_file, 'w')

f.write('# F<sub>ST</sub> windows of interest\n\n')

f.write(
'For each sample set, we first provide a summary plot of F<sub>ST</sub> bewteen '\
'resistant and susceptible samples across the genome, with F<sub>ST</sub> shown '\
'in red and the results of the 200 randomisations shown behind in grey. Windows '\
'identified as peaks are highlighted by points, colour-coded by whether they are '\
'significantly higher than expected based on the simulations (green) or not (purple). '\
'For each significant window (green points), we then provide four additional '\
'figures. The first shows the significant (*P* < 0.01) SNPs found in the region '\
'of that window, and their -log10(Pvalue) of association with phenotype. Red points '\
'indicate non-synonymous SNPs, blue points indicate all other SNPs. Point shape '\
'indicates whether the mutant allele at that SNP is associated with increased '\
'(circle) or decreased (triangle) resistance. Dark points in the centre of the plot '\
'show SNPs within the significant window, light points on the sides show SNPs in '\
'region 10,000 bp either side of the window. The second plot shows the hierarchical '\
'cluster of haplotypes in the window, with brackets indicating large haplotype '\
'clusters. The third plot shows, each haplotype cluster, indicating the number of '\
'haplotypes comprising the cluster (n), the *P* value of association of the cluster '\
'with resistance phenotype, and the SNPs that differentiate that cluster from other '\
'haplotypes in the region, colour-coded by type (red = non-synonymous, yellow = '\
'synonymous, blue = non-coding). The final figure shows the contingency tables of '\
'phenotype association for each haplotype cluster. \n\n'
)

f.write('[Legend](#Plot_legend)  \n')
for pop in populations:
	f.write(f'[{populations_italics[pop]}](#{pop.lower()})  \n')

f.write('\n___\n\n')

f.write('## Plot legend\n\n')
f.write('<a id="Plot_legend">\n\n')
f.write('![legend](./Fst_example.svg)\n\n')

f.write('\n___\n\n')


for pop in populations:
	f.write(f'## {populations_italics[pop]}\n\n')
	f.write(f'![{pop}_peak_filter](../../randomisations/Fst/{pop}_peak_filter_plot.png)\n\n&nbsp;\n\n')

	pop_windows = [w for w in all_windows if re.search(pop, w)]
	pop_windows_link = [re.sub(':', '%3A', w) for w in pop_windows]
	pop_windows_section = [re.sub(':', '_', w) for w in pop_windows]
	window_code = [re.sub('.*_([23]?[LRX]):', r'\1_', w) for w in pop_windows]
	# Getting the right overview file is a bit complicated since they were named by the middle SNP
	# rather than the edges. We need to set up a few objects to deal with this.
	overview_files = np.array([x for x in os.listdir(f'../../focal_gwas/Fst/{pop}') if re.search('.png$', x)])
	overview_file_chrom = [re.findall('^[23]?[LRX]', x)[0] for x in overview_files]
	overview_file_window_num = [int(re.findall(r'(?<=_)\d+', x)[0]) for x in overview_files]

	for i in range(len(pop_windows)):
		f.write(f'[{pop_windows[i]}](#{pop_windows_section[i].lower()})  \n')

	f.write('\n')

	for i in range(len(pop_windows)):
		f.write(f'### {pop_windows_section[i]}\n\n&nbsp;\n\n')
		# Getting the right overview file is a bit complicated since they were named by the middle SNP
		# rather than the edges.
		start, stop = [int(x) for x in re.findall(r'(?<=[:-])\d+', pop_windows[i])]
		chrom = re.findall('(?<=_)[23]?[LRX](?=:)', pop_windows[i])[0]
		overview_file = overview_files[np.isin(np.array(overview_file_chrom), chrom) &
		                               (np.array(overview_file_window_num) > start) &
		                               (np.array(overview_file_window_num) < stop)
		][0]
		f.write(f'![{pop_windows[i]}_overview](../../focal_gwas/Fst/{pop}/{overview_file})\n\n&nbsp;\n\n')

		f.write(f'![{pop_windows[i]}_dendrogram](../../haplotypes/{pop_windows_link[i]}_dendrogram.png)\n\n&nbsp;\n\n')
		f.write(f'![{pop_windows[i]}_haplotypes](../../haplotypes/{pop_windows_link[i]}.png)\n\n&nbsp;\n\n')
		haps = pd.read_csv(f'../../haplotypes/{pop_windows[i]}.csv', sep = '\t')
		all_clusters = [re.findall('cluster_\d+', x)[0] for x in haps.columns if re.search('cluster', x)]
		for c in all_clusters:
			contable = pd.crosstab(haps.phenotype, haps[c])
			if 0 not in contable.columns:
				contable.loc[:, 0] = [0, 0]
			if 1 not in contable.columns:
				contable.loc[:, 1] = [0, 0]
			if 2 not in contable.columns:
				contable.loc[:, 2] = [0, 0]
			f.write(f'#### Contingency table for {c}:\n\n')
			f.write(f'| &nbsp;&nbsp;&nbsp;  phenotype  &nbsp;&nbsp;&nbsp; |  &nbsp;&nbsp;&nbsp;  wt  &nbsp;&nbsp;&nbsp;  |  &nbsp;&nbsp;&nbsp;  het  &nbsp;&nbsp;&nbsp;  |  &nbsp;&nbsp;&nbsp;  hom  &nbsp;&nbsp;&nbsp;  |\n')
			f.write(f'| :---:     | :---: | :---: | :---: |\n')
			f.write(f'|  alive    | {contable.loc["alive", 0]} | {contable.loc["alive", 1]} | {contable.loc["alive", 2]} |\n')
			f.write(f'|  dead     | {contable.loc["dead", 0]} | {contable.loc["dead", 1]} | {contable.loc["dead", 2]} |\n')
			f.write('\n\n&nbsp;\n\n')

	f.write('___\n\n')

f.close()

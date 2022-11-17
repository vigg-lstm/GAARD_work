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

all_files = {p: np.sort(os.listdir(f'../../focal_gwas/H12/{p}')) for p in populations}

output_file = 'window_H12_summary.md'

f = open(output_file, 'w')

f.write('# H12 windows of interest\n\n')

f.write(
'For each sample set, we first provide a summary plot of the difference in '\
'H<sub>12</sub> (H<sub>12</sub> in susceptible samples subtracted from that in '\
'resistant samples) across the genome, with H<sub>12</sub> shown in green and the '\
'results of the 200 randomisations shown behind in grey. Windows identified as peaks '\
'are highlighted by points, colour-coded by whether they are significantly higher '\
'than expected based on the simulations (green) or not (purple). For each significant '\
' window (green points), we then provide its own plot showing the significant '\
'(*P* < 0.01) SNPs found in the region of that window, and their -log10(Pvalue) '\
'of association with phenotype. Red points indicate non-synonymous SNPs, blue points '\
'indicate all other SNPs. Point shape indicates whether the mutant allele at that SNP '\
'is associated with increased (circle) or decreased (triangle) resistance. Dark points '\
'in the centre of the plot show SNPs within the significant window, light points on '\
'the sides show SNPs in the region 10,000 bp either side of the window. \n\n'
)

for pop in populations:
	f.write(f'[{populations_italics[pop]}](#{pop.lower()})  \n')

f.write('\n___\n\n')

for pop in populations:
	f.write(f'## {populations_italics[pop]}\n\n')
	f.write(f'![{pop}_peak_filter](../../randomisations/H12/{re.sub("_", ".", pop)}_peak_filter_plot.png)\n\n&nbsp;\n\n')

	overview_files = [w for w in all_files[pop] if re.search('.png', w)]
	window_code = [re.sub(".png", "", w) for w in overview_files]
	pop_windows = [f'{pop}_{re.sub("_", ":", re.sub(".png", "", w))}' for w in window_code]
	pop_windows_link = [re.sub(':', '%3A', w) for w in pop_windows]
	pop_windows_section = [re.sub(':', '_', w) for w in pop_windows]

	for i in range(len(pop_windows)):
		f.write(f'[{pop_windows[i]}](#{pop_windows_section[i].lower()})  \n')

	f.write('\n')

	for i in range(len(pop_windows)):
		f.write(f'### {pop_windows_section[i]}\n\n&nbsp;\n\n')
		f.write(f'![{pop_windows[i]}_overview](../../focal_gwas/H12/{pop}/{overview_files[i]})\n\n&nbsp;\n\n')

	f.write('___\n\n')

f.close()

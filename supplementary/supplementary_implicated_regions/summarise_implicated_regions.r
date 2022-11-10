# This script will create a summary table of all the genomic regions implicated in resistance based
# on each of the genome-wide anlysis methods. 
library(data.table)
library(stringr)
library(ape)

study.pops <- setNames(nm = c('Avrankou_coluzzii_Delta',
                              'Baguida_gambiae_Delta',
                              'Baguida_gambiae_PM',
                              'Korle-Bu_coluzzii_Delta',
                              'Korle-Bu_coluzzii_PM',
                              'Madina_gambiae_Delta',
                              'Madina_gambiae_PM',
                              'Obuasi_gambiae_Delta',
                              'Obuasi_gambiae_PM'))

# Function to remove the 2La region from a table 
filter.2La <- function(regions.table){
	region.2La <- c(20524000, 42165600)
	output.table <- regions.table[chrom != '2L' |
	                              pos < region.2La[1] |
	                              pos > region.2La[2]
	]
	output.table
}

# Functions to identify all the genes found in a table of genomic regions
find.genes.from.region <- function(chrom, start.pos, end.pos, gff){
	gff[seqid == chrom & 
	    start < end.pos & 
	    end > start.pos, 
		full.gene.name] 
}

find.genes.from.regions.table <- function(regions.table, gff){
	gene.counts <- regions.table %>%
	               with(., mapply(find.genes.from.region, 
	                              chrom, 
	                              start, 
	                              end, 
	                              MoreArgs = list(gff), 
	                              SIMPLIFY = F)
	               ) %>%
	               unlist() %>%
	               table()
	if (dim(gene.counts) == 0)
		return(data.table(character(), numeric()))
	else
		return(data.table(gene.counts))
}

# Load the gff file
gff.path <- '../../data/Anopheles-gambiae-PEST_BASEFEATURES_AgamP4.9.gff3'
genes.gff <- as.data.table(read.gff(gff.path, GFF3 = T))[type == 'gene']
genes.gff[, gene.id := unique(str_extract(attributes, 'AGAP\\d{6}'))]
genes.gff[, gene.name := str_extract(attributes, '(?<=Name=)[^;]+') %>%
	                     str_to_title]
genes.gff[, full.gene.name := ifelse(is.na(gene.name), gene.id, paste(gene.name, ' (', gene.id, ')', sep = ''))]

# First, the global GWAS:
load.gwas.snp.clumps <- function(pop){
	fn <- paste('~/data/ML/GAARD_SNP', 
	            pop,
	            'classical_analysis_sibgroups/top_snp_clumps_annotated.csv',
	            sep = '/'
	)
	if (file.exists(fn))
		clumps.table <- fread(fn)
	else 
		clumps.table <- NULL
	clumps.table
}

extract.gwas.genes <- function(clumps.table){
	if (!is.null(clumps.table)){
		gwas.snp.genes <- clumps.table$genes %>%
                          strsplit('\\|') %>%
                          unlist() %>%
                          table() %>%
                          data.table()
	}
	else{
		gwas.snp.genes <- data.table(character(), numeric())
	}
	gwas.snp.genes
}

gwas.snp.clumps <- lapply(study.pops, load.gwas.snp.clumps)

# How do we report our findings?
# One option is we have a genes x population table, where each cell records whether and how a gene
# has been implicated in resisance in that population. Where there is more than one source of 
# evidence, we get a ;-separated string. 
implicated.genes.gwas <- lapply(study.pops, 
                                function(pop) {
                                    extract.gwas.genes(gwas.snp.clumps[[pop]]) %>%
		                            setnames(c('gene', pop))
                                }
                         ) %>%
                         Reduce(function(x, y) {merge(data.table(x), data.table(y), by = 'gene', all = T)}, .) %>%
                        .[, c(list(gene = gene), lapply(.SD, function(x) ifelse(is.na(x), '', 'GWAS'))), .SDcols = study.pops]

# Alternatively, we have a separate table for each population, and each table is a genes x source-
# -of-evidence table with each cell just being true/false
implicated.genes.gwas.2 <- lapply(study.pops, 
                                  function(pop) {
                                    extract.gwas.genes(gwas.snp.clumps[[pop]]) %>%
		                              setnames(c('gene', 'GWAS')) %>%
                                  .[, GWAS := GWAS > 0]
                             }
                      ) 
# I think we will end up using the second, but let's roll with both for now. 

# Next, the Fst. The results in the haplpotypes folder are derived from the significan Fst regions
fst.regions <- fread('../../haplotypes/haplotype_significance_tests.csv') %>%
                     .[, population := factor(population, levels = ..study.pops)] %>%
                     .[, c('chrom', 'start', 'end') := tstrsplit(window, ':|-')] %>%
                     .[(direction == 'alive') & (logregP < 0.05)] %>%
                     .[, .(population, chrom, start = as.numeric(start), end = as.numeric(end))] %>%
                     unique() %>%
					 split(by = 'population')

fst.regions.extra <- fread('../../haplotypes/Ace1_haplotypes/haplotype_significance_tests_ace1.csv') %>%
                     .[, population := factor(population, levels = ..study.pops)] %>%
                     .[, c('chrom', 'start', 'end') := tstrsplit(window, ':|-')] %>%
                     .[(direction == 'alive') & (logregP < 0.05)] %>%
                     .[, .(population, chrom, start = as.numeric(start), end = as.numeric(end))] %>%
                     unique() %>%
					 split(by = 'population')

fst.regions[['Korle-Bu_coluzzii_PM']] <- rbind(fst.regions[['Korle-Bu_coluzzii_PM']], fst.regions.extra[['Korle-Bu_coluzzii_PM']])
fst.regions[['Madina_gambiae_PM']] <- rbind(fst.regions[['Madina_gambiae_PM']], fst.regions.extra[['Madina_gambiae_PM']])

rm(fst.regions.extra)

implicated.genes.fst <- lapply(study.pops, 
                               function(pop) {
                                   find.genes.from.regions.table(fst.regions[[pop]], genes.gff) %>%
		                           setnames(c('gene', pop))
                               }
                        ) %>%
                        Reduce(function(x, y) {merge(data.table(x), data.table(y), by = 'gene', all = T)}, .) %>%
                        .[, c(list(gene = gene), lapply(.SD, function(x) ifelse(is.na(x), '', 'Fst'))), .SDcols = study.pops]

implicated.genes.fst.2 <- lapply(study.pops, 
                                 function(pop) {
                                     find.genes.from.regions.table(fst.regions[[pop]], genes.gff) %>%
                                     setnames(c('gene', 'Fst')) %>%
                                     .[, Fst := Fst > 0]
                                 }
                          ) 

# Get the H12 implicated genes. 
h12.p.thresh <- 0.01
h12.regions.fn <- '../../randomisations/H12/h12_filtered_windows.RDS'
h12.regions <- readRDS(h12.regions.fn) %>%
               lapply(function(D) D$diff[is.peak == T & pval < h12.p.thresh, 
                                         .(chrom = chromosome, pos = midpoint, start = startpoint, end = endpoint)]
               ) %>%
               lapply(filter.2La)

implicated.genes.h12 <- lapply(study.pops, 
                               function(pop) {
                                   dotpop <- gsub('_', '.', pop)
                                   find.genes.from.regions.table(h12.regions[[dotpop]], genes.gff) %>%
		                           setnames(c('gene', pop))
                               }
                        ) %>%
                        Reduce(function(x, y) {merge(data.table(x), data.table(y), by = 'gene', all = T)}, .) %>%
                        .[, c(list(gene = gene), lapply(.SD, function(x) ifelse(is.na(x), '', 'H12'))), .SDcols = study.pops]

implicated.genes.h12.2 <- lapply(study.pops, 
                                 function(pop) {
                                     dotpop <- gsub('_', '.', pop)
                                     find.genes.from.regions.table(h12.regions[[dotpop]], genes.gff) %>%
                                     setnames(c('gene', 'H12')) %>%
                                     .[, H12 := H12 > 0]
                                 }
                          ) 

# Get the PBS implicated genes. 
pbs.p.thresh <- 0.01
pbs.regions.fn <- '../../randomisations/PBS/pbs_filtered_windows.RDS'
pbs.regions <- readRDS(pbs.regions.fn) %>%
               lapply(function(D) D[is.peak == T & pval < pbs.p.thresh, 
                                    .(chrom = chromosome, pos = midpoint, start = startpoint, end = endpoint)]
               ) %>%
               lapply(filter.2La)

implicated.genes.pbs <- lapply(study.pops, 
                               function(pop) {
                                   dotpop <- gsub('_', '.', pop)
                                   find.genes.from.regions.table(pbs.regions[[dotpop]], genes.gff) %>%
		                           setnames(c('gene', pop))
                               }
                        ) %>%
                        Reduce(function(x, y) {merge(data.table(x), data.table(y), by = 'gene', all = T)}, .) %>%
                        .[, c(list(gene = gene), lapply(.SD, function(x) ifelse(is.na(x), '', 'PBS'))), .SDcols = study.pops]

implicated.genes.pbs.2 <- lapply(study.pops, 
                                 function(pop) {
                                     dotpop <- gsub('_', '.', pop)
                                     find.genes.from.regions.table(pbs.regions[[dotpop]], genes.gff) %>%
                                     setnames(c('gene', 'PBS')) %>%
                                     .[, PBS := PBS > 0]
                                 }
                          ) 

# Now combine the sources of evidence
all.genes <- unique(c(implicated.genes.gwas$gene, 
                      implicated.genes.fst$gene,
                      implicated.genes.h12$gene,
                      implicated.genes.pbs$gene))
ig.gwas <- implicated.genes.gwas[all.genes, ..study.pops][, lapply(.SD, function(x) replace(x, is.na(x), ''))]
ig.fst <- implicated.genes.fst[all.genes, ..study.pops][, lapply(.SD, function(x) replace(x, is.na(x), ''))]
ig.h12 <- implicated.genes.h12[all.genes, ..study.pops][, lapply(.SD, function(x) replace(x, is.na(x), ''))]
ig.pbs <- implicated.genes.pbs[all.genes, ..study.pops][, lapply(.SD, function(x) replace(x, is.na(x), ''))]
implicated.genes <- matrix(paste(as.matrix(ig.gwas), 
                                 as.matrix(ig.fst),
                                 as.matrix(ig.h12),
                                 as.matrix(ig.pbs),
                                 sep = ';'),
                           length(all.genes),
                           ncol(ig.gwas),
						   dimnames = list(NULL, study.pops)
                    ) %>%
                    gsub('^;+|;+$', '', .) %>%
                    data.table(.) %>%
					.[, gene := ..all.genes] %>%
                    setcolorder(c('gene', study.pops))

implicated.genes.2 <- lapply(study.pops,
                             function(pop){
								 Reduce(function(x,y) merge(x,y,by = 'gene', all = T),
                                        list(implicated.genes.gwas.2[[pop]],
                                             implicated.genes.fst.2[[pop]],
                                             implicated.genes.h12.2[[pop]],
                                             implicated.genes.pbs.2[[pop]])
                                 ) %>%
                                 .[, lapply(.SD, function(x) replace(x, is.na(x), F))]
                             }
                      ) 

fwrite(implicated.genes, 'full_implicated_genes_table.csv', sep = '\t')
for (pop in study.pops)
	fwrite(implicated.genes.2[[pop]], paste('implicated_genes_table_', pop, '.csv', sep = '') , sep = '\t')


# Now we want to plot a figure that summarises the genomic regions implicated by each method. 
chrom.sizes <- c('2R' = 61545105, '2L' = 49364325, '3R' = 53200684, '3L' = 41963435, 'X' = 24393108)

plot.regions.on.chromosome <- function(regions, 
	                                   chrom.sizes, gaps = 5e6,
	                                   gene.cex = 0.9, gene.col = 'grey20',
	                                   chrom.col = NULL, chrom.cex = 1.4,
	                                   chrom.offset = 0, rect.col = 'blue'){
	ce <- cumsum(chrom.sizes + c(0, 0, gaps, 0, gaps))
	cs <- ce - chrom.sizes
	plot(c(cs[1], ce[5]), c(-6.5,1.3), xlim = c(cs[1] + gaps/2, ce[5] - gaps/2), type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
	
	# Add the regions
	if (nrow(regions) > 0){
		regions[, c('genome.start', 'genome.end') := .SD + cs[.BY$chrom], 
				  .SDcols = c('start', 'end'), 
				  by = 'chrom'
		]
		# Give every region a minimum width so that it's visible on the plot
		min.region.size <- 1e6
		region.size.deficit <- min.region.size - (regions$end - regions$start)
		regions[region.size.deficit > 0, ':='(genome.start = genome.start - region.size.deficit/2,
											  genome.end = genome.end + region.size.deficit/2)]
		
		regions[, rect(genome.start, -1, genome.end, 1, 
					   border = NA, col = rect.col, lwd = 5)]
	}
	
	# Draw the outline of the chromosomes
	chrom.col <- if (!is.null(chrom.col)) chrom.col else c('2R' = 'black', '2L' = 'black', '3R' = 'black', '3L' = 'black', 'X' = 'black')
	chrom.y <- -4.5 - chrom.offset
	lines(c(ce['2R'], ce['2R'] - gaps/2, cs['2R'], cs['2R'], ce['2R'] - gaps/2, ce['2R']), 
		  c(-0.2, -1, -1, 1, 1, 0.2), lwd = 2, col = chrom.col['2R'])
	text((cs['2R'] + ce['2R'])/2, chrom.y, '2R', adj = 0.5, xpd = NA, cex = chrom.cex)
	lines(c(cs['2L'], cs['2L'] + gaps/2, ce['2L'], ce['2L'], cs['2L'] + gaps/2, cs['2L']), 
		  c(0.2, 1, 1, -1, -1, -0.2), lwd = 2, col = chrom.col['2L'])
	text((cs['2L'] + ce['2L'])/2, chrom.y, '2L', adj = 0.5, xpd = NA, cex = chrom.cex)
	lines(c(ce['3R'], ce['3R'] - gaps/2, cs['3R'], cs['3R'], ce['3R'] - gaps/2, ce['3R']), 
		  c(-0.2, -1, -1, 1, 1, 0.2), lwd = 2, col = chrom.col['3R'])
	text((cs['3R'] + ce['3R'])/2, chrom.y, '3R', adj = 0.5, xpd = NA, cex = chrom.cex)
	lines(c(cs['3L'], cs['3L'] + gaps/2, ce['3L'], ce['3L'], cs['3L'] + gaps/2, cs['3L']), 
		  c(0.2, 1, 1, -1, -1, -0.2), lwd = 2, col = chrom.col['3L'])
	text((cs['3L'] + ce['3L'])/2, chrom.y, '3L', adj = 0.5, xpd = NA, cex = chrom.cex)
	lines(c(cs['X'], cs['X'], ce['X'] - gaps/2, ce['X'], ce['X'], ce['X'] - gaps/2, cs['X']), 
		  c(-1, 1, 1, 0.2, -0.2, -1, -1), lwd = 2, col = chrom.col['X'])
	text((cs['X'] + ce['X'])/2, chrom.y, 'X', adj = 0.5, xpd = NA, cex = chrom.cex)
}

plot.implicated.regions <- function(gwas.regions,
                                    fst.regions,
                                    h12.regions,
                                    pbs.regions,
	                                title = NULL){
	par(mfrow = c(4,1), mar = c(0.7,1,2.2,1))
	if (!is.null(title))
		par(oma = c(0,0,3,0))
	plot.regions.on.chromosome(gwas.regions, chrom.sizes, rect.col = 'purple')
	mtext('GWAS', line = 0.7, col = 'purple', font = 2)
	plot.regions.on.chromosome(fst.regions, chrom.sizes, rect.col = 'orange')
	mtext('Fst', line = 0.7, col = 'orange', font = 2)
	plot.regions.on.chromosome(h12.regions, chrom.sizes, rect.col = 'darkgreen')
	mtext('H12', line = 0.7, col = 'darkgreen', font = 2)
	plot.regions.on.chromosome(pbs.regions, chrom.sizes, rect.col = 'blue')
	mtext('PBS', line = 0.7, col = 'blue', font = 2)
	if (!is.null(title))
		mtext(title, outer = T, line = 1, font = 2, cex = 1.8)
}

# For the GWAS, we need to create a table of region extents
get.clump.regions <- function(clump.table){
	if (is.null(clump.table)){
		output.table <- data.table(chrom = character(),
		                           start = numeric(),
		                           end = numeric())
	}
	else {
		output.table <- clump.table[, c('chrom', 'snp.pos') := tstrsplit(ID, ':')] %>%
		                .[, .(chrom = sub(':.*', '', .BY[[1]]), 
			                  start = min(as.numeric(snp.pos)),
			                  end = max(as.numeric(snp.pos))
			                ),
			                by = 'approx.pos'
		                ] %>%
		                .[, approx.pos := NULL]
	}
	output.table
}

gwas.regions <- lapply(gwas.snp.clumps, get.clump.regions)

for (pop in study.pops){
	filename <- paste(pop, 'implicated_regions.svg', sep = '_')
	svg(filename, width = 7, height = 4)
	plot.implicated.regions(gwas.regions[[pop]], 
							fst.regions[[pop]],
							h12.regions[[gsub('_', '.', pop)]],
							pbs.regions[[gsub('_', '.', pop)]],
							sub('_', ' ', pop))
	dev.off()
}


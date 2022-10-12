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

# First, the global GWAS:
# For this we need to be on Neptune

# Next, the Fst. The windows implicated by the Fst analysis can be found in the focal_gwas folder
fst.regions <- fread(paste('../../haplotypes/haplotype_significance_tests.csv', sep = '/')) %>%
                     .[, population := factor(population, levels = ..study.pops)] %>%
                     .[, c('chrom', 'start', 'end') := tstrsplit(window, ':|-')] %>%
                     .[(direction == 'alive') & (logregP < 0.05)] %>%
                     .[, .(population, chrom, start, end)] %>%
                     unique() %>%
					 split(by = 'population')

# Now, for each of those regions, we need to work out which genes they contain. 
gff.path <- '../../data/Anopheles-gambiae-PEST_BASEFEATURES_AgamP4.9.gff3'
genes.gff <- as.data.table(read.gff(gff.path, GFF3 = T))[type == 'gene']
genes.gff[, gene.id := unique(str_extract(attributes, 'AGAP\\d{6}'))]
genes.gff[, gene.name := str_extract(attributes, '(?<=Name=)[^;]+') %>%
	                     str_to_title]
genes.gff[, full.gene.name := ifelse(is.na(gene.name), gene.id, paste(gene.name, ' (', gene.id, ')', sep = ''))]

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

# How do we report our findings
# One option is we have a genes x population table, where each cell records whether and how a gene
# has been implicated in resisance in that population. Where there is more than one source of 
# evidence, we get a ;-separated string. 
implicated.genes <- lapply(study.pops, 
                           function(pop) {
                               find.genes.from.regions.table(fst.regions[[pop]], genes.gff) %>%
		                       setnames(c('gene', pop))
                           }
                    ) %>%
                    Reduce(function(x, y) {merge(data.table(x), data.table(y), by = 'gene', all = T)}, .) %>%
                    {ifelse(is.na(.), 'No', 'Fst')}

# Alternatively, we have a separate table for each population, and each table is a genes x source-
# -of-evidence table with each cell just being true/false
implicated.genes.2 <- lapply(study.pops, 
                             function(pop) {
                                 find.genes.from.regions.table(fst.regions[[pop]], genes.gff) %>%
                                 setnames(c('gene', 'Fst')) %>%
                                 .[, Fst := Fst > 0]
                             }
                      ) 

# I think we will end up going for the second one, but let's run with both for now. 




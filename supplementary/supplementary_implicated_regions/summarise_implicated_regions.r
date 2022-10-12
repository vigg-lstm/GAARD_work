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
	if (file.exists(fn)){
		clumps.table <- fread(fn)
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

# How do we report our findings
# One option is we have a genes x population table, where each cell records whether and how a gene
# has been implicated in resisance in that population. Where there is more than one source of 
# evidence, we get a ;-separated string. 
implicated.genes.gwas <- lapply(study.pops, 
                                function(pop) {
                                    load.gwas.snp.clumps(pop) %>%
		                            setnames(c('gene', pop))
                                }
                         ) %>%
                         Reduce(function(x, y) {merge(data.table(x), data.table(y), by = 'gene', all = T)}, .) %>%
                        .[, c(list(gene = gene), lapply(.SD, function(x) ifelse(is.na(x), '', 'GWAS'))), .SDcols = study.pops]

# Alternatively, we have a separate table for each population, and each table is a genes x source-
# -of-evidence table with each cell just being true/false
implicated.genes.gwas.2 <- lapply(study.pops, 
                                  function(pop) {
                                      load.gwas.snp.clumps(pop) %>%
		                              setnames(c('gene', 'GWAS')) %>%
                                  .[, GWAS := GWAS > 0]
                             }
                      ) 
# I think we will end up using the second, but let's roll with both for now. 

# Next, the Fst. The results in the haplpotypes folder are derived from the significan Fst regions
fst.regions <- fread(paste('../../haplotypes/haplotype_significance_tests.csv', sep = '/')) %>%
                     .[, population := factor(population, levels = ..study.pops)] %>%
                     .[, c('chrom', 'start', 'end') := tstrsplit(window, ':|-')] %>%
                     .[(direction == 'alive') & (logregP < 0.05)] %>%
                     .[, .(population, chrom, start, end)] %>%
                     unique() %>%
					 split(by = 'population')

# Now, for each of those regions, we need to work out which genes they contain. 
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

# Now combine the two sources of evidence
all.genes <- unique(c(implicated.genes.gwas$gene, implicated.genes.fst$gene))
ig1 <- implicated.genes.gwas[all.genes, ..study.pops][, lapply(.SD, function(x) replace(x, is.na(x), ''))]
ig2 <- implicated.genes.fst[all.genes, ..study.pops][, lapply(.SD, function(x) replace(x, is.na(x), ''))]
implicated.genes <- matrix(paste(as.matrix(ig1), 
                                 as.matrix(ig2),
                                 sep = ';'),
                           length(all.genes),
                           ncol(ig1),
						   dimnames = list(NULL, study.pops)
                    ) %>%
                    gsub('^;|;$', '', .) %>%
                    data.table(.) %>%
					.[, gene := ..all.genes] %>%
                    setcolorder(c('gene', study.pops))

implicated.genes.2 <- lapply(study.pops,
                             function(pop){
                                 merge(implicated.genes.gwas.2[[pop]],
                                       implicated.genes.fst.2[[pop]],
                                       by = 'gene', 
                                       all = T 
                                 ) %>%
                                 .[, lapply(.SD, function(x) replace(x, is.na(x), F))]
                             }
                      ) 



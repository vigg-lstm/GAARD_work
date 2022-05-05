library(data.table)
library(magrittr)
library(lme4)
library(stringr)

phenotypes.fn <- '../../data/combined/sample_phenotypes.csv'

sib.groups.fn <- '../../NGSrelate/full_relatedness/sib_group_table.csv'

phenotypes <- fread(phenotypes.fn, sep = '\t', key = 'specimen', colClasses = c(phenotype = 'factor'))
sib.groups <- fread(sib.groups.fn, sep = '\t', key = 'sample.name')

# We will restrict the analysis to individuals in sib groups

sib.groups$phenotype <- phenotypes[sib.groups$sample.name, phenotype]
sib.groups[, population := paste(location, species, insecticide, sep = '_')]

sib.group.test <- function(pop, verbose = T){
	if (verbose)
		cat(pop, '\n')
	sibs <- sib.groups[(population == pop), ]
	if (length(unique(sibs$sib.group.id)) == 1){
		if (verbose)
			cat('Only one sib group present. Returning NA p-values.\n')
		return(NA)
	}
	# Remove samples that are not in sib groups within this dataset. This is super hacky, but I couldn't
	# find a more elegant way of doing this. 
	sibs <- sibs[, .SD[ifelse(rep(.N, .N) > 1, 1:.N, 0), ], by = 'sib.group.id']
	fisher.p <- sibs %>%
	            {table(.$phenotype, .$sib.group.id)} %>%
	            fisher.test()
	if (verbose)
		print(fisher.p)
	fisher.p$p.value
}

populations <- setNames(nm = sort(unique(sib.groups$population)))
sig.testing <- sapply(populations, sib.group.test)

cat('\n\nFisher testing of sib groups P.values:\n\n')
print(sig.testing)


library(data.table)
library(Biostrings)
library(future.apply)
plan(tweak(multisession, workers = 20))

# Load the genome
cat('Loading AgamP4 genome.\n')
genome.path <- '~/data/ML/GAARD/data/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa'
genome <- readDNAStringSet(genome.path)
names(genome) <- sub(' .*', '', names(genome))
main.chromosomes <- genome[c('2R', '2L', '3R', '3L', 'X')]
chrom.size <- width(main.chromosomes)
names(chrom.size) <- names(main.chromosomes)

# Set the populations. 
pops <- setNames(nm = c('Avrankou_coluzzii_Delta', 
                        'Baguida_gambiae_Delta',
                        'Korle-Bu_coluzzii_Delta',
                        'Madina_gambiae_Delta',
                        'Obuasi_gambiae_Delta',
                        'Baguida_gambiae_PM',
                        'Korle-Bu_coluzzii_PM',
                        'Madina_gambiae_PM',
                        'Obuasi_gambiae_PM'))

# A function to calculate the logistic regression P value. 
logreg.test <- function(genotype, phenotype){
	genotype <- as.numeric(genotype)
	null.model <- glm(phenotype ~ 1, family = 'binomial')
	model1 <- glm(phenotype ~ genotype, family = 'binomial')
	anova(model1, null.model, test = 'Chisq')[['Pr(>Chi)']][2]
}

# The threshold MAF that we will use to filter SNPs
MAF.thresh <- 5

# Load phenotype data
phenotype.filename <- '~/data/ML/GAARD_SNP/pre_QC/metadata/sample_phenotypes.csv'
cat('Loading phenotype data from ', phenotype.filename, '.\n', sep = '')
phenotype.table <- fread(phenotype.filename, key = 'specimen')

# Load the VObs meta data (which includes sex calls)
vobs.meta <- fread('~/data/ML/GAARD_SNP/data/all_samples.samples.meta.csv')
males <- vobs.meta[sex_call == 'M', partner_sample_id]

# A function to get the residuals of a value relative to phenotype
get.residuals <- function(x, phenotypes){
	x[phenotypes == 'alive'] <- x[phenotypes == 'alive'] - mean(x[phenotypes == 'alive'])
	x[phenotypes == 'dead'] <- x[phenotypes == 'dead'] - mean(x[phenotypes == 'dead'])
	x
}

# A function to load asaia data
load.bracken <- function(sample.name, bracken.folder){
	fn <- paste(bracken.folder, '/', sample.name, '/', sample.name, '_genus.bracken', sep = '')
	fread(fn, key = 'name')
}

# A function to calculate minor allele frequency
maf <- function(genotype){
	half.hap.num <- length(genotype)
	half.hap.num - abs(sum(genotype) - half.hap.num)
}

load.SNPs <- function(pop){
	# Load SNP data after accessisibility filtering. 
	SNP.filename <- paste('..', pop, 'csvs/filtered_SNPs_acc.csv', sep = '/')
	cat('Loading SNP data from ', SNP.filename, '.\n', sep = '')
	snp.table <- fread(SNP.filename)

	# Remove any males
	these.males <- intersect(males, names(snp.table))
	if (length(these.males) > 0)
		snp.table <- snp.table[, (these.males) := NULL]
	
	# Get the phenotypes for these samples. 
	sample.names <- colnames(snp.table)[3:ncol(snp.table)]
	phenotypes <- setNames(as.factor(phenotype.table[sample.names, phenotype]), sample.names)
	
	# Generate a SNP id for each SNP. 
	snp.table[, ID := paste(Chrom, Pos, sep = ':')]
	
	# Calculate the MAF for each SNP
	cat('Filtering by minor allele frequency.\n')
	snp.table[, MAF := apply(.SD, 1, maf), .SDcols = sample.names]
	# Remove SNPs where maf is lower than 5
	snp.table = snp.table[(MAF >= MAF.thresh), ]
	
	# Load the bracken Asaia data
	bracken.folder = paste('~/data/ML/GAARD_bracken/', pop, sep = '')
	bracken.by.sample <- lapply(sample.names, load.bracken, bracken.folder = bracken.folder)
	names(bracken.by.sample) <- sample.names
	bracken.by.sample.summarised <- list()
	for (i in 1:length(bracken.by.sample))
		bracken.by.sample.summarised[[names(bracken.by.sample)[i]]] <- bracken.by.sample[[i]][, {x <- list(name, taxonomy_id, fraction_total_reads); names(x) <- c('name', 'taxonomy_id', names(bracken.by.sample)[i]); x}]
	# Merge all of the bracken tables together
	bracken <- Reduce(function(X1, X2) merge(X1, X2, by = c('name', 'taxonomy_id'), all = T), bracken.by.sample.summarised)
	bracken[is.na(bracken)] <- 0
	asaia <- setNames(as.numeric(bracken[name == 'Asaia', ..sample.names]), colnames(bracken[, ..sample.names]))
	# Get the asaia residuals
	asaia.residuals <- get.residuals(asaia, phenotypes = phenotypes)
	# Write a function that will get the correlation p.value against asaia, using the residuals. 
	get.asaia.cor <- function(x){
		x.residuals <- get.residuals(x, phenotypes = phenotypes)
		cor.test(x.residuals, asaia.residuals)$p.value
	}
	# Filter out reads correlated with Asaia. Need to do this on the whole genome this time. 
	snp.table[, asaia.cor := apply(.SD, 1, get.asaia.cor), .SDcols = sample.names]
	snp.table <- snp.table[asaia.cor > 0.05, ]
	
	# Calculate the logistic regression P-value
	cat('Calculating logistic regression P-value.\n')
	snp.table[, logregP := future_apply(.SD, 1, logreg.test, phenotype = phenotypes), .SDcols = sample.names]
	snp.table$logregP[is.na(snp.table$logregP)] <- 1
	
	list(snp.table, phenotypes)
}

snp.tables <- lapply(pops, load.SNPs)

for (pop in pops)
	saveRDS(snp.tables[[pop]], paste('filtered_snp_tables_asaia_', pop, '.rds', sep = ''))


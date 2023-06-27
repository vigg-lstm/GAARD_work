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

# The threshold MAF that we will use to filter SNPs
MAF.thresh <- 2

# Load phenotype data
phenotype.filename <- '~/data/ML/GAARD_SNP/pre_QC/metadata/sample_phenotypes.csv'
cat('Loading phenotype data from ', phenotype.filename, '.\n', sep = '')
phenotype.table <- fread(phenotype.filename, key = 'specimen')

# Load the VObs meta data (which includes sex calls)
vobs.meta <- fread('~/data/ML/GAARD_SNP/data/all_samples.samples.meta.csv')
males <- vobs.meta[sex_call == 'M', partner_sample_id]

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
	# Remove SNPs where maf is lower than the threshold
	snp.table = snp.table[(MAF >= MAF.thresh), ]
	
	list(snp.table, phenotypes)
}

snp.tables <- lapply(pops, load.SNPs)

save.image('filtered_snp_tables.Rdata')

for (pop in pops)
	saveRDS(snp.tables[[pop]], paste('filtered_snp_tables_', pop, '.rds', sep = ''))


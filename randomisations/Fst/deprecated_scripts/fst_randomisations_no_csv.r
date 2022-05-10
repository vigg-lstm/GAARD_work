# This is the same as fst_randomisations.r but only save .Rdata environments rather than csv 
# files. This is so as to avoid overly large output files when there are a large number of
# randomisations. 
# We also remove the hmm calculation, as I don't think we need this for all 10000 randomisations.
# We can always come back and do those calculations. later. 
library(data.table)
library(RColorBrewer)
library(magrittr)
library(stringr)
library(reticulate)

# Load the SNP data and functions from the Fst analysis
load('~/data/ML/GAARD_SNP/fst_analyses/fst_1000_window.Rdata')

arg.values <- commandArgs(trailingOnly=T)
cat('Running with populations:', arg.values[1], '\n', sep = '\n')

{if (arg.values[1] == 'all'){
	populations <- names(snp.tables)
}
else {
	populations <- strsplit(arg.values[1], ',')[[1]]
	snp.tables[!(names(snp.tables) %in% populations)] <- NULL
}}

num.randomisations <- as.integer(arg.values[2])

# Load a conda env to run the reticulate code (which is inside the windowed.fst function)
use_condaenv('gaard')
allel <- import('allel')

# Print out the variables from that environment so we know what they are
cat('window.size = ', window.size, '\n')

# Print out the windowed.fst and windowed.fst.wrapper functions to show what they do
cat('\nwindowed.fst :\n')
print(windowed.fst)

# Load the phenotype randomisations
phenotype.filename <- '../phenotype_10000_randomisations.csv'
cat('Loading randomised phenotype data from ', phenotype.filename, '.\n', sep = '')
phenotype.table <- fread(phenotype.filename, key = 'specimen')
phenotype.table$population <- gsub('\\.', '_', phenotype.table$population)
randomisations <- colnames(phenotype.table)[grepl('^r\\d+$', colnames(phenotype.table))]

samples.by.pop <- with(phenotype.table, split(specimen, population))

cat('Calculating window positions\n')
calculate.window.pos <- function(pop){
	# Get the number of windows that will be found on each chromosome
	chrom.snp.counts <- tapply(snp.tables[[pop]][[1]]$Chrom, snp.tables[[pop]][[1]]$Chrom, length)
	# We take the floor rather than the ceiling because the movingFst function will only report a value
	# for full sized windows. 
	window.num <- floor(chrom.snp.counts / window.size)
	window.end <- cumsum(window.num)
	window.start <- window.end - window.num + 1
	snp.tables[[pop]][[1]]$window <- unlist(sapply(names(window.num), function(chrom) rep(window.start[chrom]:window.end[chrom], each = window.size, length.out = chrom.snp.counts[chrom])))
	window.pos <- snp.tables[[pop]][[1]][, .(window.chrom = unique(Chrom), window.pos = floor(median(Pos))), by = window][, window := NULL]
	window.pos
}

window.pos <- lapply(setNames(nm = populations), calculate.window.pos)

# Write a function to calculate allele counts from genotypes
get.allele.counts <- function(genotypes){
	mut.counts <- as.integer(rowSums(genotypes))
	wt.counts <- as.integer(ncol(genotypes) * 2 - mut.counts)
	data.frame('wt'= wt.counts, 'mut'= mut.counts)
}

# Calculate the total allele count in each population, this will speed things up later

cat('Calculating total allele counts\n')
allele.counts.total <- lapply(setNames(nm = populations), function(pop) snp.tables[[pop]][[1]][, get.allele.counts(.SD), .SDcols = samples.by.pop[[pop]]])

# A function to calculate windowed Fst. Same as for the main data, except an extra element of parallelisation.
windowed.fst <- function(pop, phenotypes, window.size){
	
	samples.alive <- names(phenotypes)[phenotypes == 'alive']
	
	cat('\t\tCalculating allele counts.\n')
	allele.counts.alive <- snp.tables[[pop]][[1]][, .(Chrom, get.allele.counts(.SD)), .SDcols = samples.alive]
	# For the dead, we can just subtract the alive from the total, saving us a bit of time. 
	allele.counts.dead <- data.table(Chrom = snp.tables[[pop]][[1]]$Chrom, wt = allele.counts.total[[pop]]$wt - allele.counts.alive$wt, mut = allele.counts.total[[pop]]$mut - allele.counts.alive$mut)
	
	# Calculate moving Fst. For the HMM, we want it as a column matrix. 
	cat('\t\tCalculating Fst.\n')
	moving.Fst <- lapply(setNames(nm = c('2L', '2R', '3L', '3R', 'X')), function(chrom) matrix(allel$moving_patterson_fst(allele.counts.alive[Chrom == chrom, .(wt, mut)], allele.counts.dead[Chrom == chrom, .(wt, mut)], as.integer(window.size))))
	
	windowed.data <- data.table(moving.Fst = unlist(moving.Fst))
	windowed.data
}


# Write a function to pull out a phenotype vector 
get.phenotype.iteration <- function(pop, iteration){
	cat('\t', iteration, '\n', sep = '')
	with(phenotype.table[population == pop, c('specimen', get('iteration'))], setNames(get(iteration), specimen))
}

# Write a new version of the windowed.fst.wrapper function
windowed.fst.wrapper <- function(pop, randomisations, window.size){
	cat('\nRunning', length(randomisations), 'fst calculations for', pop, '\n')
	windowed.data <- lapply(setNames(nm = c('phenotype', randomisations)), function(r) windowed.fst(pop, get.phenotype.iteration(pop, r), window.size))
	do.call(cbind, windowed.data)
}

fst.table <- lapply(setNames(nm = populations), windowed.fst.wrapper, randomisations = randomisations[1:num.randomisations], window.size = window.size)

# Save the image (after deleting the SNP table, which is huge and available in a different workspace)
# If a third argument was passed to the script, use that in the output filename, if not, use arg 1. 
rm(snp.tables)
{if (is.na(arg.values[3]))
	save.image(paste('fst_randomisations_', arg.values[1], '.Rdata', sep = ''))
else
	save.image(paste('fst_randomisations_', arg.values[3], '.Rdata', sep = ''))
}


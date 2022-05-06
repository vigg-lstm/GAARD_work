library(data.table)
library(RColorBrewer)
library(magrittr)
library(stringr)
library(reticulate)

# Load the SNP data and functions from the Fst analysis
load('~/data/ML/GAARD_SNP/fst_analyses/fst_1000_window/fst_1000_window.Rdata')

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

# Print out some variables from that environment so we know what they are
cat('window.size = ', window.size, '\n')

# Load the phenotype randomisations
phenotype.filename <- '../phenotype_10000_randomisations.csv'
cat('Loading randomised phenotype data from ', phenotype.filename, '.\n', sep = '')
phenotype.table <- fread(phenotype.filename, key = 'specimen')
phenotype.table$population <- gsub('\\.', '_', phenotype.table$population)
randomisations <- colnames(phenotype.table)[grepl('^r\\d+$', colnames(phenotype.table))]

samples.by.pop <- with(phenotype.table, split(specimen, population))

calculate.window.pos <- function(pop){
	# Get the number of windows that will be found on each chromosome
	chrom.snp.counts <- tapply(snp.tables[[pop]][[1]]$Chrom, snp.tables[[pop]][[1]]$Chrom, length)
	window.num <- ceiling(chrom.snp.counts / window.size)
	window.end <- cumsum(window.num)
	window.start <- window.end - window.num + 1
	snp.tables[[pop]][[1]]$window <- unlist(sapply(names(window.num), function(chrom) rep(window.start[chrom]:window.end[chrom], each = window.size, length.out = chrom.snp.counts[chrom])))
	window.pos <- snp.tables[[pop]][[1]][, .(window.chrom = unique(Chrom), window.pos = floor(median(Pos))), by = window][, window := NULL]
	window.pos
}

cat('Calculating window positions\n')
window.pos <- lapply(setNames(nm = populations), calculate.window.pos)
window.num.table <- sapply(window.pos, function(x) table(x$window.chrom))

# Calculate the total allele count in each population, this will speed things up later
cat('Calculating total allele counts\n')
allele.counts.total <- lapply(setNames(nm = populations), function(pop) snp.tables[[pop]][[1]][, get.allele.counts(.SD), .SDcols = samples.by.pop[[pop]]])

# Write a function to pull out a phenotype vector 
get.phenotype.iteration <- function(pop, iteration){
	with(phenotype.table[population == pop, c('specimen', get('iteration'))], setNames(get(iteration), specimen))
}

# A function to calculate windowed Fst. Similar as for the main data, but we pass phenotype as a separate argument. 
windowed.fst.with.dropped.samples <- function(pop, iteration, dropped.samples.matrix, window.size){
	
	cat('\tCalculating Fst for population', pop, 'with phenotype randomisation"', iteration, '"\n')
	phenotypes <- get.phenotype.iteration(pop, iteration)
	snp.table <- snp.tables[[pop]][[1]]
	samples.alive <- names(phenotypes)[phenotypes == 'alive']
	samples.dead <- names(phenotypes)[phenotypes == 'dead']
	
	cat('\t\tCalculating allele counts.\n')
	total.allele.counts.alive <- snp.table[, .(Chrom, get.allele.counts(.SD)), .SDcols = samples.alive]
	# For the dead, we can just subtract the alive from the total, saving us a bit of time. 
	total.allele.counts.dead <- data.table(Chrom = snp.table$Chrom, wt = allele.counts.total[[pop]]$wt - total.allele.counts.alive$wt, mut = allele.counts.total[[pop]]$mut - total.allele.counts.alive$mut)
	
	# For each permutation, we get the allele counts of the dropped samples and subtract them from the total
	dropped.samples.fst <- function(dropped.samples){
		cat('\t\tCalculating Fst after dropping samples ', paste(dropped.samples, collapse = ', '), '\n')
		dropped.samples.alive <- intersect(dropped.samples, samples.alive)
		dropped.samples.dead <- intersect(dropped.samples, samples.dead)
		dropped.samples.allele.counts.alive <- snp.table[, get.allele.counts(.SD), .SDcols = dropped.samples.alive]
		dropped.samples.allele.counts.dead <- snp.table[, get.allele.counts(.SD), .SDcols = dropped.samples.dead]
		# If no alive samples were dropped, the dropped.samples.allele.counts.alive table will have no rows. 
		if (nrow(dropped.samples.allele.counts.alive)){
			allele.counts.alive <- (total.allele.counts.alive[, .(wt, mut)] - dropped.samples.allele.counts.alive) %>%
								   .[, Chrom := total.allele.counts.alive$Chrom]
		}
		else {
			allele.counts.alive <- total.allele.counts.alive
		}
		if (nrow(dropped.samples.allele.counts.dead)){
		allele.counts.dead <- (total.allele.counts.dead[, .(wt, mut)] - dropped.samples.allele.counts.dead) %>%
		                       .[, Chrom := total.allele.counts.dead$Chrom]
		}
		else {
			allele.counts.dead <- total.allele.counts.dead
		}
		# Calculate moving.Fst by chromosome
		moving.Fst <- lapply(setNames(nm = unique(window.pos[[1]]$window.chrom)), 
		                     function(chrom) matrix(allel$moving_patterson_fst(allele.counts.alive[Chrom == chrom, .(wt, mut)], allele.counts.dead[Chrom == chrom, .(wt, mut)], as.integer(window.size)))) %>%
		              # Pad with NAs where the final window of a chrom had < 1000 SNPs
					  {lapply(names(.), 
	                          function(chrom) c(.[[chrom]], rep(NA, window.num.table[chrom, pop] - length(.[[chrom]]))))} %>%
		              unlist()
		moving.Fst
	}
	
	cat('\t\tCalculating mean Fst of sib permutations\n')
	all.moving.Fst <- apply(dropped.samples.matrix, 1, dropped.samples.fst)
	mean.moving.Fst <- apply(all.moving.Fst, 1, mean)
	
	windowed.data <- data.table(moving.Fst = mean.moving.Fst)
	windowed.data
}

# Write a new version of the windowed.fst.permuted.sibdrop function. num.permutations is the number of sibgroup permuations. 
windowed.fst.permuted.sibdrop <- function(pop, randomisations, window.size, num.permutations = 100){
	permutations <- create.permutations(sib.groups.list[[pop]])
	# The last row of the permutations matrix is the total size of each sib group. We can use it
	# to get the start index of each group
	sib.group.start.index <- c(1, permutations[nrow(permutations), -ncol(permutations)]) %>%
	                         cumsum
	
	if (nrow(permutations) > num.permutations)
		permutations <- permutations[sample(1:nrow(permutations), num.permutations), ]
	sample.keep.index <- t(t(permutations) + sib.group.start.index) - 1
	samples.to.drop <- t(apply(sample.keep.index, 1, function(x) sib.groups.list[[pop]][-x, sample.name]))
	
	cat('\nRunning', length(randomisations), 'fst calculations for', pop, '\n')
	windowed.data <- lapply(setNames(nm = c('phenotype', randomisations)), 
	                        function(r) windowed.fst.with.dropped.samples(pop, r, samples.to.drop, window.size))
	do.call(cbind, windowed.data)
}

set.seed(42)

fst.tables <- lapply(setNames(nm = populations), 
                     windowed.fst.permuted.sibdrop, 
                     randomisations = randomisations[1:num.randomisations], 
                     window.size = window.size)

# Save the image (after deleting the SNP table, which is huge and available in a different workspace)
# If a third argument was passed to the script, use that in the output filename, if not, use arg 1. 
rm(snp.tables)
{if (is.na(arg.values[3]))
	save.image(paste('fst_randomisations_', arg.values[1], '.Rdata', sep = ''))
else
	save.image(paste('fst_randomisations_', arg.values[3], '.Rdata', sep = ''))
}


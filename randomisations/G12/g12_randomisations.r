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
hmmlearn <- import('hmmlearn.hmm')

# Print out the variables from that environment so we know what they are
cat('window.size = ', window.size, '\n')
cat('num.hmm.states = ', num.hmm.states, '\n')

# Print out the windowed.fst and windowed.fst.wrapper functions to show what they do
cat('\nwindowed.fst :\n')
print(windowed.fst)

# Load the phenotype randomisations
phenotype.filename <- '../phenotype_randomisations.csv'
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
windowed.fst <- function(pop, phenotypes, window.size, num.hmm.states){
	
	samples.alive <- names(phenotypes)[phenotypes == 'alive']
	
	cat('\t\tCalculating allele counts.\n')
	allele.counts.alive <- snp.tables[[pop]][[1]][, .(Chrom, get.allele.counts(.SD)), .SDcols = samples.alive]
	# For the dead, we can just subtract the alive from the total, saving us a bit of time. 
	allele.counts.dead <- data.table(Chrom = snp.tables[[pop]][[1]]$Chrom, wt = allele.counts.total[[pop]]$wt - allele.counts.alive$wt, mut = allele.counts.total[[pop]]$mut - allele.counts.alive$mut)
	
	# Calculate moving Fst. For the HMM, we want it as a column matrix. 
	cat('\t\tCalculating Fst.\n')
	moving.Fst <- lapply(setNames(nm = c('2L', '2R', '3L', '3R', 'X')), function(chrom) matrix(allel$moving_patterson_fst(allele.counts.alive[Chrom == chrom, .(wt, mut)], allele.counts.dead[Chrom == chrom, .(wt, mut)], as.integer(window.size))))
	
	# Now fit the HMM
	cat('\t\tFitting HMM.\n')
	run.hmm <- function(Fst, variance, n_states = as.integer(num.hmm.states)){
		# I tried training the variance, but it ends up increasing it, which is unlikely to be
		# a true reflection of the data since we have used a dataset-wide calculation of the 
		# variance, which should be a maximum estimate of the within-state variance. Also, it is
		# mostly the non-0 states that have increased variance after training, so you end up 
		# with transitions to states when the data are closer to 0 (because the variance of those
		# is more permissive). 
		means <- matrix(seq(-1, 1, 1/(n_states - 1)*2), n_states, 1)
		# The transition probability between states 
		trans_prob_diff <- 1e-2
		# The probability of remaining in a state. 
		trans_prob_nochange <- 1 - trans_prob_diff * (n_states - 1)
		transmat <- matrix(trans_prob_diff, n_states, n_states)
		diag(transmat) <- trans_prob_nochange
		#
		covars <- matrix(variance*n_states, n_states, 1)
		#
		model <- hmmlearn$GaussianHMM(n_states, 
									  covariance_type='diag', 
									  n_iter=as.integer(0), 
									  init_params='',
									  params='')
		 
		model$means_ <- means
		model$covars_ <- covars
		model$transmat_ <- transmat
		
		model$fit(Fst)
		# Need the -1 because python uses 0 indexing, so the state assigned to the windows
		# with 0 mean Fst will be one smaller that the output of "which" below. 
		base_state <- which(means == 0) - 1
		(model$predict(Fst) - base_state) / base_state
	}
	
	# Do a first pass of the hmm where variance is calculated across all sites
	fst.variance <- var(unlist(moving.Fst))
	hmm.output <- lapply(moving.Fst, run.hmm, variance = fst.variance)
	# Now do a second pass where we calculate the varianc after normalising by the hmm from 
	# the first pass.
	# One problem is that the fst function doesn't calculate a value for windows with fewer
	# than the full complement of SNPs, so unless there is a multiple of window.size SNPs on
	# a chromosome then there will be one window fewer in the hmm output than in the window 
	# means. So we need to potentially add an NA at the end of each hmm.output chromosome.
	new.fst.variance <- (unlist(moving.Fst) - unlist(hmm.output)) %>%
	                    var()
	hmm.output <- lapply(moving.Fst, run.hmm, variance = new.fst.variance)
	
	# We take the floor of the median because data.table is a bit silly here, if some results 
	# are integers and others are floats, it complains that all the output isn't of the same
	# type rather than coverting the integer to a float. 
	windowed.data <- data.table(moving.Fst = unlist(moving.Fst),
	                            hmm.output = unlist(hmm.output))
	windowed.data
}


# Write a function to pull out a phenotype vector 
get.phenotype.iteration <- function(pop, iteration){
	cat('\t', iteration, '\n', sep = '')
	with(phenotype.table[population == pop, c('specimen', get('iteration'))], setNames(get(iteration), specimen))
}

# Write a new version of the windowed.fst.wrapper function
windowed.fst.wrapper <- function(pop, randomisations, window.size, num.hmm.states){
	cat('\nRunning', length(randomisations), 'fst calculations for', pop, '\n')
	lapply(setNames(nm = c('phenotype', randomisations)), function(r) windowed.fst(pop, get.phenotype.iteration(pop, r), window.size, num.hmm.states))
}

randomised.windowed.data <- lapply(setNames(nm = populations), windowed.fst.wrapper, randomisations = randomisations[1:num.randomisations], window.size = window.size, num.hmm.states = num.hmm.states)

for (pop in names(randomised.windowed.data)){
	full.table <- do.call(cbind, randomised.windowed.data[[pop]])
	fst.table <- cbind(window.pos[[pop]], round(full.table[, grepl('Fst$', colnames(full.table)), with = F], 3))
	hmm.table <- cbind(window.pos[[pop]], full.table[, grepl('hmm.output$', colnames(full.table)), with = F])
	fwrite(fst.table, paste(pop, '_randomised_Fst.csv', sep = ''), sep = '\t')
	fwrite(hmm.table, paste(pop, '_randomised_hmm.csv', sep = ''), sep = '\t')
}

# Save the image (after deleting the SNP table, which is huge and available in a different workspace)
# If a third argument was passed to the script, use that in the output filename, if not, use arg 1. 
rm(snp.tables)
{if (is.na(arg.values[3]))
	save.image(paste('fst_randomisations_', arg.values[1], '.Rdata', sep = ''))
else
	save.image(paste('fst_randomisations_', arg.values[3], '.Rdata', sep = ''))
}


# We just want to do the haplotype association for the Ace1 region in those sample sets
# where we explicitly removed the Ace1 region for the purposes of looking for significant
# windows. This is so that we can avoid having a blank in the Ace1 region when asking
# which methods support an association in which genomic region. 
library(data.table)
library(lme4)
library(stringr)
library(ape)

study.pops <- c('Korle-Bu_coluzzii_PM',
                'Madina_gambiae_PM')

python.env <- '~/miniconda3/envs/gaard/bin/python3'

snp.table.fn <- '~/data/ML/GAARD_SNP/filtered_snp_tables/filtered_snp_tables.Rdata'
load(snp.table.fn)
snp.tables <- snp.tables[study.pops]

# Load the results of Fst-based window selection
load('../../randomisations/Fst/fst_filtered_windows.Rdata')
filtered.tables <- filtered.tables[study.pops]

filter.ace1 <- function(regions.table){
	region.ace1 <- c(3100000, 3950000)
	output.table <- regions.table[window.chrom == '2R' &
	                              window.pos >= region.ace1[1] |
	                              window.pos <= region.ace1[2]
	]
	output.table
}

# Now keep only the windows inside the Ace1 region
for (pop in study.pops)
	filtered.tables[[pop]] <- filter.ace1(filtered.tables[[pop]])


# Get the number of windows that will be found on each chromosome
window.pos <- list()
for (pop in study.pops){
	chrom.snp.counts <- tapply(snp.tables[[pop]][[1]]$Chrom, snp.tables[[pop]][[1]]$Chrom, length)
	window.num <- ceiling(chrom.snp.counts / window.size)
	window.end <- cumsum(window.num)
	window.start <- window.end - window.num + 1
	snp.tables[[pop]][[1]]$window <- unlist(sapply(names(window.num), function(chrom) rep(window.start[chrom]:window.end[chrom], each = window.size, length.out = chrom.snp.counts[chrom])))
	window.pos[[pop]] <- snp.tables[[pop]][[1]][, .(window.chrom = unique(Chrom), window.pos = floor(median(Pos))), by = window][, window := NULL]
	window.pos[[pop]]$window <- 1:nrow(window.pos[[pop]])
}

# Now find out which SNPs are found in each of the significant windows
get.window.start.end <- function(chrom, pos, window.pos, snp.table){
	window.pos.this.chrom <- window.pos[window.chrom == chrom]
	window.id <- window.pos.this.chrom[window.pos == pos, window]
	previous.window.id <- window.id - 1
	next.window.id <- window.id + 1
	# Get the window's start point. If it's the first window on the chromosome, then the start point
	# is 1. Otherwise, we take the median of the first snp in this window and the last snp of the 
	# previous one. 
	if (window.id == min(window.pos.this.chrom$window)){
		window.start.point <- 1
	}
	else {
		last.snp.of.previous.window <- max(snp.table[window == window.id - 1, Pos])
		first.snp.of.this.window <- min(snp.table[window == window.id, Pos])
		window.start.point <- floor(mean(c(last.snp.of.previous.window, first.snp.of.this.window)))
	}
	# Now do the same for the end point
	if (window.id == max(window.pos.this.chrom$window)){
		window.end.point <- chrom.size[chrom]
	}
	else {
		last.snp.of.this.window <- max(snp.table[window == window.id, Pos])
		first.snp.of.next.window <- min(snp.table[window == window.id + 1, Pos])
		window.end.point <- ceiling(mean(c(last.snp.of.this.window, first.snp.of.next.window)))
	}
	return(data.table(
		chrom = chrom, 
		start = window.start.point, 
		end = window.end.point, 
		name = paste(chrom, pos, sep = ':')
	))
}

window.ranges <- list()
for (pop in study.pops){
	window.ranges[[pop]] <- mapply(get.window.start.end,
								   filtered.tables[[pop]]$window.chrom, 
								   filtered.tables[[pop]]$window.pos, 
								   MoreArgs = list(window.pos = window.pos[[pop]],
	                                               snp.table = snp.tables[[pop]][[1]][, .(window, Pos)]),
								   SIMPLIFY = F
						    ) %>%
						    rbindlist() %>%
						    setkey(name)
}

# Function to perform logistic regression for genotype-phenotype association.
logreg.test <- function(genotype, phenotype){
	genotype <- as.numeric(genotype)
	null.model <- glm(phenotype ~ 1, family = 'binomial')
	model1 <- glm(phenotype ~ genotype, family = 'binomial')
	pval <- anova(model1, null.model, test = 'Chisq')[['Pr(>Chi)']][2]
	direction <- levels(phenotype)[sign(model1$coefficients['genotype'])/2+1.5]
	data.table(logregP = pval, direction = direction)
}

sig.test.list <- list()
for (pop in study.pops){
	sig.tests.this.pop <- list()
	# Loop through the windows. This code means nrow = 0 results in not entering the loop
	for (i in seq(length.out = nrow(window.ranges[[pop]]))){
		window.range <- paste(
			window.ranges[[pop]][i, chrom],
			paste(window.ranges[[pop]][i, .(start, end)], collapse = '-'),
			sep = ':'
		)
		# Assign haplotypes using python script
		assign.hap.command <- paste(python.env, '../assign_haplotypes.py', pop, window.range, '../..')
		system(assign.hap.command)
		# Load the results
		haplotypes <- fread(paste(pop, '_', window.range, '.csv', sep = ''), sep = '\t', key = 'sample_name')
		# Perform association analysis
		sig.test <- grep('cluster_', names(haplotypes), value = T) %>%
		            setNames(nm = .) %>%
		            lapply(function(cluster) haplotypes[, cbind(cluster.size = sum(get(cluster)), logreg.test(get(cluster), as.factor(phenotype)))]) %>%
		            rbindlist(idcol = 'cluster.name')
		# If there were no haplotype clusters, sig.test will be an empty table and will disappear in the rbindlist
		if (nrow(sig.test) == 0)
			sig.test <- data.table(cluster.name = NA, cluster.size = NA, logregP = NA, direction = NA)
		sig.tests.this.pop[[window.range]] <- sig.test
	}
	sig.test.list[[pop]] <- rbindlist(sig.tests.this.pop, idcol = 'window')
}
sig.tests <- rbindlist(sig.test.list, idcol = 'population')

fwrite(sig.tests, file = 'haplotype_significance_tests_ace1.csv', sep = '\t')

# We don't actually need all the intermediary files
for (pop in study.pops){
	system(paste('rm ', pop, '*', sep = ''))
}


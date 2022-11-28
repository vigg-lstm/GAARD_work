library(data.table)
library(fdrtool)
library(stringr)
library(future.apply)
plan(tweak(multisession, workers = 30))

set.seed(42)

# Logistic regression
logreg.test <- function(genotype, phenotype){
	genotype <- as.numeric(genotype)
	null.model <- glm(phenotype ~ 1, family = 'binomial')
	model1 <- glm(phenotype ~ genotype, family = 'binomial')
	pval <- anova(model1, null.model, test = 'Chisq')[['Pr(>Chi)']][2]
	direction <- sign(model1$coefficients['genotype'])
	list(pval = pval, direction = direction)
}

# A function to calcualte genotype frequencies from an allele frequency
gen.freq <- function(p){
	c((1-p)^2, 
	   2*p*(1-p), 
	   p^2)
}

# Load the significant SNPs that we want to simulate
tsg.fn <- '../target_site_genotypes.csv'
example.sig.snps <- fread(tsg.fn, colClasses = c(phenotype = 'factor')) %>%
                    .[location == 'Avrankou', .(phenotype, Vgsc = Vgsc.1527T_C, Cyp4j5 = Cyp4j5.43F_A)]
example.freqs <- example.sig.snps[, lapply(.SD, function(x) mean(x)/2), by = 'phenotype']
example.genotype.freqs <- example.freqs[, lapply(.SD, gen.freq), by = 'phenotype']

# A function to simulate significant SNPs with a given number of samples. 
simulate.sig.snps <- function(sample.size, genotype.freqs){
	# For each locus, generate a random table of genotypes, in two cohorts (dead and alive) based on the 
	# genotype frequencies given to the function. So we are saying, "let's assume these genotype frequencies
	# (separate frequencies given for dead and alive) and recreate a sample set of size 'sample.size' of 
	# 50/50 dead/alive using those genotype frequencies. 
	genotypes <- genotype.freqs[, 
		lapply(
			.SD, 
			function(x) 
				sample(c(0,1,2), ..sample.size/2, replace = T, prob = x)
		), 
		by = 'phenotype'
	]
	# Calculate the P-values for those loci
	sig.snp.pval <- genotypes[, 
		lapply(
			.SD, 
			function(x) logreg.test(x, phenotype)$pval
		), 
		.SDcols = c('Vgsc', 'Cyp4j5')
	]
	sig.snp.pval[is.na(sig.snp.pval)] <- 1
	sig.snp.pval
}

# A function to simulate the SNPs and then run statistical testing and FDR correction on them.
# n is the number of phenotype-assocated snps to model (this will impact on the FDR, since the 
# more truly associated SNPs there are in the genome, the easier it will be to identify them 
# without including too many false positives). 
simulate.statistical.testing <- function(n, sample.size, snps.in.genome, genotype.freqs){
	simulated.snps <- rbindlist(replicate(n, simulate.sig.snps(sample.size, genotype.freqs), simplify = F))
	# Get a uniform P-value distribution to represent the rest of the genome
	other.snps <- runif(snps.in.genome, 0, 1)
	snp.type <- factor(c(rep('associated', n), rep('unassociated', snps.in.genome)))
	snp.fdr <- simulated.snps[, 
		lapply(
			.SD, 
			function(x) {
				c(x, ..other.snps) %>%
				fdrtool(statistic = 'pvalue', plot = F, verbose = F) %>%
				{.$qval < 0.05} %>%
				factor(levels = c(T,F)) %>%
				table(sig.fdr = .,snp.type)
			}
		)
	]
	# Two of the columns are redundant, so reduce the table and simplify the column names here
	snp.fdr[, .(sig.fdr = Vgsc.sig.fdr, snp.type = Vgsc.snp.type, Vgsc = Vgsc.N, Cyp4j5 = Cyp4j5.N)]
}

num.snps.in.genome <- 7000000
sample.sizes <- c(100, 200, 300, 500, 750, 1000)
num.sig.snps.to.model <- c(1, 3, 10, 30, 100)
num.iterations <- 500


simulation.results <- simulation.summary <- list()
for (sample.size in sample.sizes){
	cat('Sample size', sample.size, '\n')
	sname <- paste('S', sample.size, sep = '')
	simulation.results[[sname]] <- simulation.summary[[sname]] <- list()
	
	for (n in num.sig.snps.to.model){
		cat('\tSimulating', n, 'significant SNPs\n')
		nname <- paste('N', n, sep = '')
		
		simulation.results[[sname]][[nname]] <- future_replicate(
			num.iterations, 
			simulate.statistical.testing(n, sample.size, num.snps.in.genome, example.genotype.freqs),
			simplify = F
		) %>%
		rbindlist()
		
		summed.table <- simulation.results[[sname]][[nname]][,
			lapply(.SD, sum), by = c('sig.fdr', 'snp.type')
		]
		
		this.summary <- summed.table[, .(
			# Sensitivity
			sensitivity = .SD[sig.fdr == T & snp.type == 'associated', .(Vgsc, Cyp4j5)]/(num.iterations*n),
			# Specifity should be 95%, since we set the FDR to 5%
			specificity = .SD[sig.fdr == T, .SD[snp.type == 'associated', .(Vgsc, Cyp4j5)]/lapply(.SD[, .(Vgsc, Cyp4j5)], sum)],
			# False positive rate
			fpr = .SD[sig.fdr == T & snp.type == 'unassociated', .(Vgsc, Cyp4j5)]/(num.iterations*num.snps.in.genome)
		)]
				
		simulation.summary[[sname]][[nname]] <- this.summary
	}
}

cat('Summary of statistical power:\n\n')

print(simulation.summary)

save.image('sample_size_calculation.Rdata')



library(data.table)
library(RColorBrewer)
library(apcluster)
library(qqman)
library(fdrtool)
library(lme4)
library(Biostrings)
library(ape)
library(future.apply)
library(HardyWeinberg)
library(stringr)
plan(tweak(multisession, workers = 20))

study.pop = 'Korle-Bu_coluzzii_PM'
study.id.table <- fread('../../data/study_ids.csv', key = 'population')
study.id <- study.id.table[study.pop, study_id]
chrom.order <- c('2R', '2L', '3R', '3L', 'X')

# Logistic regression
logreg.test <- function(genotype, phenotype){
	genotype <- as.numeric(genotype)
	null.model <- glm(phenotype ~ 1, family = 'binomial')
	model1 <- glm(phenotype ~ genotype, family = 'binomial')
	pval <- anova(model1, null.model, test = 'Chisq')[['Pr(>Chi)']][2]
	direction <- sign(model1$coefficients['genotype'])
	list(pval = pval, direction = direction)
}

# Hardy Weinberg
HW.test <- function(genotype){
	wt.hom <- sum(genotype == 0)
	het <- sum(genotype == 1)
	mut.hom <- sum(genotype == 2)
	HWExact(setNames(c(wt.hom, het, mut.hom), c('AA', 'AB', 'BB')), alternative = 'greater', verbose = F)$pval
}

# The threshold MAF that we will use to filter SNPs
MAF.thresh <- 5

# Load SNP data
snp.data <- readRDS(paste('../../filtered_snp_tables_asaia/filtered_snp_tables_asaia_', study.pop, '.rds', sep = ''))

# Load the sib groups information
sib.groups <- fread('~/data/GAARD_repo/GAARD_work/NGSrelate/full_relatedness/sib_group_table.csv', sep = '\t')
samples.to.remove <- sib.groups[keep == F, sample.name]
sample.names <- setdiff(names(snp.data[[2]]), samples.to.remove)
columns.to.keep <- setdiff(names(snp.data[[1]]), samples.to.remove)
snp.table <- snp.data[[1]][, ..columns.to.keep]
phenotypes <- snp.data[[2]][sample.names]
rm(snp.data)

# Load phenotype data
phenotype.filename <- '~/data/ML/GAARD_SNP/pre_QC/metadata/sample_phenotypes.csv'
cat('Loading phenotype data from ', phenotype.filename, '.\n', sep = '')
phenotype.table <- fread(phenotype.filename, key = 'specimen')
phenotypes <- setNames(as.factor(phenotype.table[sample.names, phenotype]),
                       nm = sample.names)

# Write a function that will get the residuals of a value relative to phenotype
get.residuals <- function(x, phen = phenotypes){
	x[phen == 'alive'] <- x[phen == 'alive'] - mean(x[phen == 'alive'])
	x[phen == 'dead'] <- x[phen == 'dead'] - mean(x[phen == 'dead'])
	x
}

# Function to calculate minor allele frequency. We need to recalculate this since we have removed sibs
maf <- function(genotype){
	half.hap.num <- length(genotype)
	half.hap.num - abs(sum(genotype) - half.hap.num)
}
# Calculate the MAF for each SNP
cat('Filtering by minor allele frequency.\n')
snp.table[, MAF := apply(.SD, 1, maf), .SDcols = sample.names]
# Remove SNPs where maf is lower than 5
snp.table <- snp.table[(MAF >= MAF.thresh), ]

# Get the logistic regression P-value
cat('Calculating logistic regression P-value.\n')
significance.test <- future_apply(snp.table[, ..sample.names], 
                                  1, 
                                  logreg.test, 
                                  phenotype = phenotypes)
snp.table$logregP <- sapply(significance.test, '[[', 'pval')
snp.table$logreg.direction <- sapply(significance.test, '[[', 'direction')
# Set the NA P values to 1
snp.table$logregP[is.na(snp.table$logregP)] <- 1
# Get the FDRs 
cat('Calculating false discovery rates.\n')
snp.table[, logregFDR := fdrtool(logregP, statistic = 'pvalue', plot = F)$qval]

# Create a numeric equivalent for the Chrom column
snp.table[, numChrom := setNames(1:5, chrom.order)[as.character(Chrom)]]



# We filter by logistic regression FDR
sig.snps <- setorder(snp.table[logregFDR < 0.01, ], logregP)
# Recode the direction of the effect as a factor
sig.snps$logreg.direction <- as.factor(c('Protects', 'Null', 'Weakens')[sig.snps$logreg.direction + 2])
sig.snps[, hwe := apply(.SD, 1, HW.test), .SDcols = sample.names]
# 469 SNPs were significant. The vast majority of these are in the Ace1 region.

# We want to take the top 1000 SNPs. Instead of taking the top 1000 SNPs, we get the phenotype-association 
# P value for the 1000th SNP and keep all SNPs with a P-value at least that small. This is to avoid 
# chromosome bias in the case when many SNPs have the same P-value (which would lead to alphabetically 
# smaller chromosomes being favoured). 
threshold.p <- snp.table[order(logregP), ][1000, logregP]
top.snps <- snp.table[logregP <= threshold.p, ]
top.snps$logreg.direction <- as.factor(c('Protects', 'Null', 'Weakens')[top.snps$logreg.direction + 2])
top.snps[, hwe := apply(.SD, 1, HW.test), .SDcols = sample.names]

# Define some colour vectors
chrom.colours <- c('2R' = 'dodgerblue4', '2L' = 'cornflowerblue', '3R' = 'orchid4', '3L' = 'orchid', 'X' = 'palegreen4')
palette <- colorRampPalette(rev(brewer.pal(7, name = 'RdBu')))(101)
palbreaks <- seq(-1, 1, length.out = length(palette) + 1)

# Load the genome
cat('Loading AgamP4 genome.\n')
genome.path <- '~/data/ML/GAARD/data/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa'
genome <- readDNAStringSet(genome.path)
names(genome) <- sub(' .*', '', names(genome))
main.chromosomes <- genome[c('2R', '2L', '3R', '3L', 'X')]
chrom.size <- width(main.chromosomes)
names(chrom.size) <- names(main.chromosomes)

# Write a function to draw the chromosomes on an existing plotting device
add.chromosomes <- function(gaps, cs, ce, gene.cex = 0.9, gene.col = 'grey20', point.cex = 1.2, point.col = 'grey30', chrom.col = NULL, chrom.cex = 1.4, chrom.offset = 0){
	plot(c(cs[1], ce[5]), c(-6.5,1.3), xlim = c(cs[1] + gaps/2, ce[5] - gaps/2), type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
	# show the kdr region
	kdr.region.mean <- cs['2L'] + mean(c(2358158, 2431617))
	points(kdr.region.mean, 0, pch = 19, cex = point.cex, col = point.col)
	text(kdr.region.mean, -1.7, 'Vgsc', srt = 45, adj = 1.1, xpd = T, cex = gene.cex, col = gene.col, font = 2)
	# show the Ace1 gene region
	ace1.region.mean <- cs['2R'] + mean(c(3484107, 3495790))
	points(ace1.region.mean, 0, pch = 19, cex = point.cex, col = point.col)
	text(ace1.region.mean, -1.7, 'Ace1', srt = 45, adj = 1.1, xpd = T, cex = gene.cex, col = gene.col, font = 2)
	# show the CYP6 region
	cyp6.region.mean <- cs['2R'] + mean(c(28463000, 28568000))
	points(cyp6.region.mean, 0, pch = 19, cex = point.cex, col = point.col)
	text(cyp6.region.mean, -1.7, 'Cyp6aa1-\nCyp6p2 ', srt = 45, adj = 1.1, xpd = T, cex = gene.cex, col = gene.col, font = 2)
	# show the GST region
	gst.region.mean <- cs['3R'] + mean(c(28580000, 28605000))
	points(gst.region.mean, 0, pch = 19, cex = point.cex, col = point.col)
	text(gst.region.mean, -1.7, 'Gstu4-\nGste3 ', srt = 45, adj = 1.1, xpd = T, cex = gene.cex, col = gene.col, font = 2)
	# show the CYP6M2-Z1 region
	cyp6m2.z1.region.mean <- cs['3R'] + mean(c(6900000, 7030000))
	points(cyp6m2.z1.region.mean, 0, pch = 19, cex = point.cex, col = point.col)
	text(cyp6m2.z1.region.mean, -1.7, 'Cyp6m2-\nCyp6z1 ', srt = 45, adj = 1.1, xpd = T, cex = gene.cex, col = gene.col, font = 2)
	# show the CYP9K1 region
	cyp9k1.region.mean <- cs['X'] + mean(c(15222000, 15257000))
	points(cyp9k1.region.mean, 0, pch = 19, cex = point.cex, col = point.col)
	text(cyp9k1.region.mean, -1.7, 'Cyp9k1', srt = 45, adj = 1, xpd = T, cex = gene.cex, col = gene.col, font = 2)
	# Plot the outline of the chromosomes
	chrom.col <- if (!is.null(chrom.col)) chrom.col else c('2R' = 'black', '2L' = 'black', '3R' = 'black', '3L' = 'black', 'X' = 'black')
	chrom.y <- -8.5 - chrom.offset
	lines(c(ce['2R'], ce['2R'] - gaps/2, cs['2R'], cs['2R'], ce['2R'] - gaps/2, ce['2R']), 
		  c(-0.2, -1, -1, 1, 1, 0.2), lwd = 2, col = chrom.col['2R'])
	text((cs['2R'] + ce['2R'])/2, chrom.y, '2R', adj = 0.5, xpd = NA, cex = chrom.cex)
	lines(c(cs['2L'], cs['2L'] + gaps/2, ce['2L'], ce['2L'], cs['2L'] + gaps/2, cs['2L']), 
		  c(0.2, 1, 1, -1, -1, -0.2), lwd = 2, col = chrom.col['2L'])
	text((cs['2L'] + ce['2L'])/2, chrom.y, '2L', adj = 0.5, xpd = NA, cex = chrom.cex)
	lines(c(ce['3R'], ce['3R'] - gaps/2, cs['3R'], cs['3R'], ce['3R'] - gaps/2, ce['3R']), 
		  c(-0.2, -1, -1, 1, 1, 0.2), lwd = 2, col = chrom.col['3R'])
	text((cs['3R'] + ce['3R'])/2, chrom.y, '3R', adj = 0.5, xpd = NA, cex = chrom.cex)
	lines(c(cs['3L'], cs['3L'] + gaps/2, ce['3L'], ce['3L'], cs['3L'] + gaps/2, cs['3L']), 
		  c(0.2, 1, 1, -1, -1, -0.2), lwd = 2, col = chrom.col['3L'])
	text((cs['3L'] + ce['3L'])/2, chrom.y, '3L', adj = 0.5, xpd = NA, cex = chrom.cex)
	lines(c(cs['X'], cs['X'], ce['X'] - gaps/2, ce['X'], ce['X'], ce['X'] - gaps/2, cs['X']), 
		  c(-1, 1, 1, 0.2, -0.2, -1, -1), lwd = 2, col = chrom.col['X'])
	text((cs['X'] + ce['X'])/2, chrom.y, 'X', adj = 0.5, xpd = NA, cex = chrom.cex)
}

# Write a function to plot the SNPs along the genome
plot.snps.along.genome <- function(snp.table, filename, gaps = 5000000, chrom.sizes = chrom.size, window.size = 1000000, ...){
	# get a vector of end points for each chromosome (ie: cumulative sizes + gaps)
	ce <- cumsum(chrom.sizes + c(0, 0, gaps, 0, gaps))
	# get a vector of start points for each chromosome (ie: cumulative sizes + gaps)
	cs <- ce - chrom.sizes
	# Get a column that shows the position along the whole genome ,rather than just on the chromosome
	if (!('Genome.pos' %in% names(snp.table)))
		snp.table$Genome.pos <- snp.table$Pos + cs[as.character(snp.table$Chrom)]
	# Work out some genomic windows
	genomic.windows <- mapply(function(x, y) c(seq(x, y, window.size), y), cs, ce)
	
	# Now count the number of reads occuring in each window
	ff <- function(x) hist(subset(snp.table, Chrom == x)$Genome.pos, breaks = genomic.windows[[x]], plot = F)
	window.hist <- lapply(names(genomic.windows), ff)
	names(window.hist) <- names(genomic.windows)
	window.kmer.counts <- lapply(window.hist, function(x) x$counts)
	window.mids <- lapply(window.hist, function(x) x$mids)
	
	# Create the plot
	if (!missing(filename)){
		file.width = 6.5
		file.height = 2.5
		if (grepl('\\.eps', filename))
			postscript(filename, width = file.width, height = file.height, horizontal = F, onefile = FALSE, paper = "special")
		else if (grepl('\\.png', filename))
			png(filename, width = file.width, height = file.height, units = 'in', res = 600)
		else if (grepl('\\.tif', filename))
			tiff(filename, width = file.width, height = file.height, units = 'in', res = 600)
	}
	layout(matrix(c(rep(1,2), rep(2,1)), nrow = 3, ncol = 1))
	par(mar = c(0,4,1,0.8), mgp = c(2, 0.7, 0), ...) 
	# First add the density plot for the significant kmers
	max.counts <- max(unlist(window.kmer.counts))
	plot(c(cs[1], ce[5]), c(0, max.counts), xlim = c(cs[1] + gaps/2, ce[5] - gaps/2), type = 'n', bty = 'n', xaxt = 'n', xlab = '', ylab = 'kmer counts')
	mapply(function(x, y) lines(x, y + 2.5, col = 'darkred', lwd = 1.5), window.mids, window.kmer.counts)
	
	# Now plot all chromosomes with, the position of each of the four detox gene regions and Ace1
	par(mar = c(2,4,0,0.8), mgp = c(2, 0.7, 0), ...) 
	add.chromosomes(gaps = gaps, cs = cs, ce = ce)
	if (!missing(filename))
		dev.off()
	invisible(list(snp.table, window.hist))
}

#plot.snps.along.genome(top.snps, window.size = 100000, family = 'Arial')
plot.snps.along.genome(top.snps, window.size = 100000, filename = 'plotted_snp_density_nosibs.png', family = 'Arial')
# Almost all on Ace1. Bit on 2La. 

# Now write a function to draw the kmers on a heatmap
snp.heatmap <- function (x, snp.pos, filename, filter.fail.marker, chrom.col = NULL, legend = T, revC = F, scale = c("row", "column", "none"), na.rm = TRUE, segment.lwd = 0.1, ...) {
	x <- copy(x)
    scale <- if (missing(scale)) 
        "none"
    else match.arg(scale)
	di <- dim(x)
	if (di[1] != di[2]) stop('input matrix should be square')
    nr <- di[1]
	# Make sure the chromosomes are in the order we want
	snp.pos$Chrom <- factor(as.character(snp.pos$Chrom), levels = c('2R', '2L', '3R', '3L', 'X'))
	snp.order <- order(snp.pos$Chrom, snp.pos$Pos)
	snp.pos <- snp.pos[snp.order, ]
	x <- x[snp.order, snp.order]
	# Scale if required
    if (scale == "row") {
        x <- sweep(x, 1, rowMeans(x, na.rm = na.rm), check.margin = FALSE)
        sx <- apply(x, 1, sd, na.rm = na.rm)
        x <- sweep(x, 1, sx, "/", check.margin = FALSE)
    }
    else if (scale == "column") {
        x <- sweep(x, 2, colMeans(x, na.rm = na.rm), check.margin = FALSE)
        sx <- apply(x, 2, sd, na.rm = na.rm)
        x <- sweep(x, 2, sx, "/", check.margin = FALSE)
    }
	if (!missing(filename)){
		if (legend)
			file.width = 7.2
		else
			file.width = 6.5
		file.height = 8
		if (grepl('\\.eps', filename))
			postscript(filename, width = file.width, height = file.height, horizontal = F, onefile = FALSE, paper = "special")
		else if (grepl('\\.png', filename))
			png(filename, width = file.width, height = file.height, units = 'in', res = 600)
		else if (grepl('\\.tif', filename))
			tiff(filename, width = file.width, height = file.height, units = 'in', res = 600)
	}
	lhei <- c(4,0.8)
	if (legend){
		lmat <- matrix(c(1,2,3,3), 2, 2, byrow = T)
		lwid <- c(10,1)
		layout(lmat, heights = lhei, widths = lwid)
	}
	else {
		lmat <- matrix(1:2, 2, 1)
		layout(lmat, heights = lhei)
	}
	side.margin <- ifelse (missing(filter.fail.marker), 0.3, 1)
    par(mar = c(0.1, side.margin, 0.3, 0), xpd = NA, ...)
    if (revC)
        x <- x[, nr:1]
    image(1:nr, 1:nr, x, xlim = 0.5 + c(0, nr), ylim = 0.5 + c(0, nr), 
	      axes = FALSE, xlab = "", ylab = "", breaks = palbreaks, col = palette)
	coords.1 <- par('usr')	
	xrange.1 <- coords.1[2] - coords.1[1]
	# Add markers for filter-failing SNPs if required
	if (!missing(filter.fail.marker)){
		filter.fail.marker <- filter.fail.marker[snp.order]
		if (revC)
			filter.fail.marker <- filter.fail.marker[nr:1]
		points(rep(coords.1[1] - xrange.1/100, sum(filter.fail.marker)), which(filter.fail.marker), pch = -9658, cex = 0.6)
	}
	
	# Plot the legend if necessary
	if (legend){
		par(mar = c(0,0,0,0))
		plot(0:1, 0:1, type = 'n', bty = 'n', xlab = '', ylab = '', xaxt = 'n', yaxt = 'n')
		image(seq(0.1, 0.7, length.out = length(palette)), seq(0.5, 0.9, length.out = length(palette) + 1), matrix(palbreaks[2:length(palbreaks)] - diff(palbreaks)/2, length(palette), length(palette), byrow = T), col = palette, breaks = palbreaks, add = T)
		text(rep(0.85, 4), c(0.94, 0.9, 0.7, 0.5), c('r', '1', '0', '-1'), cex = c(1.4, 1, 1, 1))
	}
	
	# Plot the chromosomes
	# get a vector of end points for each chromosome (ie: cumulative sizes + gaps)
	gaps <- 5000000
	ce <- cumsum(chrom.size + c(0, 0, gaps, 0, gaps))
	# get a vector of start points for each chromosome (ie: cumulative sizes + gaps)
	cs <- ce - chrom.size
	top.margin <- 2.5
	par(mar = c(2, side.margin, top.margin, 0), xpd = NA)
	add.chromosomes(gaps = gaps, cs = cs, ce = ce, point.cex = 1, gene.cex = 0.7, chrom.cex = 1.2)
	coords.2 <- par('usr')
	xrange.2 <- coords.2[2] - coords.2[1]
	# Write a function to convert old x coords to new x coords
	f <- if (legend) 1.1035 else 1
	convert <- function(x) (((x - coords.1[1])/xrange.1) * xrange.2 + coords.2[1])/f
	# And work out the top end of the current plot
	top.end <- par('cxy')[2] * par('cex') * par('lheight') * top.margin + coords.2[4]
	# Add arrows linking the genom picture to the correct positions on the heatmap. 
	# Get a column that shows the position along the whole genome ,rather than just on the chromosome
	if (!('Genome.pos' %in% names(snp.pos)))
		snp.pos$Genome.pos <- snp.pos$Pos + cs[as.character(snp.pos$Chrom)]
	# For every snp, draw a line from the matrix to the chromosome.
	chrom.col <- if (!is.null(chrom.col)) chrom.col else c('2R' = 'black', '2L' = 'black', '3R' = 'black', '3L' = 'black', 'X' = 'black')
	segments(snp.pos$Genome.pos, 1.1, convert(1:nrow(x)), top.end, lwd = segment.lwd, col = chrom.col[snp.pos$Chrom])
	if (!missing(filename))
		dev.off()
}

# Get correlation matrix after taking the residuals accounting for phenotype.
top.snp.residuals <- t(apply(top.snps[, ..sample.names], 1, get.residuals))
top.snp.residuals.cor.matrix <- corSimMat(top.snp.residuals)

#x11(); snp.heatmap(top.snp.residuals.cor.matrix, top.snps[, c('Chrom', 'Pos')], filter.fail.marker = top.snps$hwe <= 0.05, chrom.col = chrom.colours, family = 'Arial', segment.lwd = 0.6)
snp.heatmap(top.snp.residuals.cor.matrix, top.snps[, c('Chrom', 'Pos')], filter.fail.marker = top.snps$hwe <= 0.05, chrom.col = chrom.colours, filename = 'top_snp_heatmap_nosibs.png', family = 'Arial', segment.lwd = 0.6)

plot.significance.along.genome <- function(significance.value, snp.pos, filename, chrom.col = NULL, lwd = 0.3, ...){
	# Sanity check
	if (length(significance.value) != nrow(snp.pos))
		stop('Number of points in significance.value is different to rows in snp.pos.')
	# Make sure the chromosomes are in the order we want
	snp.pos$Chrom <- factor(as.character(snp.pos$Chrom), levels = c('2R', '2L', '3R', '3L', 'X'))
	snp.order <- order(snp.pos$Chrom, snp.pos$Pos)
	snp.pos <- snp.pos[snp.order, ]
	significance.value <- significance.value[snp.order]
	if (!missing(filename)){
		file.width = 6.5
		file.height = 5
		if (grepl('\\.eps', filename))
			postscript(filename, width = file.width, height = file.height, horizontal = F, onefile = FALSE, paper = "special")
		else if (grepl('\\.png', filename))
			png(filename, width = file.width, height = file.height, units = 'in', res = 600)
		else if (grepl('\\.tif', filename))
			tiff(filename, width = file.width, height = file.height, units = 'in', res = 600)
	}
	lhei <- c(2,0.8)
	lmat <- matrix(1:2, 2, 1)
	layout(lmat, heights = lhei)
    par(mar = c(0, 3, 1, 0), mgp = c(2, 0.6, 0), xpd = NA, ...)
	chrom.col <- if (!is.null(chrom.col)) chrom.col else c('2R' = 'black', '2L' = 'black', '3R' = 'black', '3L' = 'black', 'X' = 'black')
	point.cex = 0.5
	plot(-log10(significance.value), bty = 'n', xlab = '', ylab = '-log(P)', xaxt = 'n', pch = 21, cex = point.cex, col = 'grey50', bg = chrom.col[snp.pos$Chrom], lwd = point.cex/2)
	coords.1 <- par('usr')	
	xrange.1 <- coords.1[2] - coords.1[1]
	# Plot the chromosomes
	# get a vector of end points for each chromosome (ie: cumulative sizes + gaps)
	gaps <- 5000000
	ce <- cumsum(chrom.size + c(0, 0, gaps, 0, gaps))
	# get a vector of start points for each chromosome (ie: cumulative sizes + gaps)
	cs <- ce - chrom.size
	top.margin <- 2
	par(mar = c(2, 3, top.margin, 0), xpd = NA, ...)
	add.chromosomes(gaps = gaps, cs = cs, ce = ce, point.cex = 0.7, gene.cex = 0.7, chrom.cex = 1.2)
	coords.2 <- par('usr')
	xrange.2 <- coords.2[2] - coords.2[1]
	# Write a function to convert old x coords to new x coords
	f <- 1
	convert <- function(x) (((x - coords.1[1])/xrange.1) * xrange.2 + coords.2[1])/f
	# And work out the top end of the current plot
	top.end <- par('cxy')[2] * par('cex') * par('lheight') * top.margin + coords.2[4]
	# Now add arrows linking the plot to the correct positions on the heatmap. 
	# Get a column that shows the position along the whole genome ,rather than just on the chromosome
	if (!('Genome.pos' %in% names(snp.pos)))
		snp.pos$Genome.pos <- snp.pos$Pos + cs[as.character(snp.pos$Chrom)]
	# For every snp, draw a line from the plot to the chromosome.
	segments(snp.pos$Genome.pos, 1.15, convert(1:length(significance.value)), top.end, lwd = lwd, col = chrom.col[snp.pos$Chrom])
	if (!missing(filename))
		dev.off()
}

#x11(); plot.significance.along.genome(top.snps$logregP, top.snps[, c('Chrom', 'Pos')], chrom.col = chrom.colours, lwd = 0.3, family = 'Arial')
plot.significance.along.genome(top.snps$logregP, top.snps[, c('Chrom', 'Pos')], chrom.col = chrom.colours, lwd = 0.3, filename = 'sig_along_genome_nosibs.png', family = 'Arial')

# Let's look for clumps of SNPs. 
top.snps$approx.pos <- paste(top.snps$Chrom, round(top.snps$Pos / 100000), sep = ':')
snp.clumps <- table(top.snps$approx.pos)
snp.clumps <- sort(snp.clumps[snp.clumps > 2], decreasing = T)
top.snp.clumps <- top.snps[approx.pos %in% names(snp.clumps), c(.N, lapply(.SD, mean)), by = approx.pos, .SDcols = sample.names]
# Now get various stats on snp clumps rather than individual snps. 
top.snp.clumps[, Chrom :=  sub(':.*', '', approx.pos)]

# We can have a bit of a closer look at possible SNPs (or SNP groups) of interest
# Other than Ace1, there are no clumps

# Print out a quick summary:
cat('\nThere were no clumps with more than 10 SNPs.\n')

# Now lets split data by Ace1
greater.Dup1.region <- c(3000000, 4000000)
sig.snps$In.Dup1 <- with(sig.snps, Chrom == '2R' & Pos > greater.Dup1.region[1] & Pos < greater.Dup1.region[2])
top.snps$In.Dup1 <- with(top.snps, Chrom == '2R' & Pos > greater.Dup1.region[1] & Pos < greater.Dup1.region[2])

# Let's try just looking at non-Ace1 SNPs
non.ace1.snps <- top.snps[(!In.Dup1), ]
#plot.significance.along.genome(non.ace1.snps$logregP, non.ace1.snps[, c('Chrom', 'Pos')], chrom.col = chrom.colours, lwd = 0.3, family = 'Arial')
plot.significance.along.genome(non.ace1.snps$logregP, non.ace1.snps[, c('Chrom', 'Pos')], chrom.col = chrom.colours, lwd = 0.3, filename = 'sig_along_genome_nonAce1.png', family = 'Arial')
# We just get a bit of a 2La signal and everything else is mostly scattered. 

# Let's plot the significance just along Ace1
ace1.snps <- top.snps[(In.Dup1), ]
gene.positions <- fread('~/data/ML/GAARD/data/gene_regions.csv', col.names = c('AGAP', 'Chrom', 'start', 'end'), key = 'AGAP')
Dup.region.genes <- paste('AGAP00', c(1352:1358, 1360:1368), sep = '')
Dup.gene.positions <- gene.positions[Dup.region.genes, c('AGAP', 'start', 'end')]
Dup.gene.positions$Name <- ''
ace1.gene <- 'AGAP001356'
Dup.gene.positions[ace1.gene, 'Name'] <- 'Ace1'
Del1.FM <- matrix(c(3501850,3598750,3502150,3599050), 2, 2)
Del2.FM <- matrix(c(3539300,3573450,3539600,3573750), 2, 2)
Del3.FM <- matrix(c(3535850,3618700,3536150,3619000), 2, 2)
Dup.coords <- data.table(Dup.name = c('Dup1', 'Del1', 'Del2', 'Del3'),
                         start = c(3436926, Del1.FM[1, 2], Del2.FM[1, 2], Del3.FM[1, 2]), 
                         end = c(3639836, Del1.FM[2,1], Del2.FM[2,1], Del3.FM[2,1]),
                         key = 'Dup.name')

# Now draw the ace1 gene with its exons
draw.gene.model <- function(gene.positions, dup.positions, del.positions = NULL, y = 0, text.cex = 0.8){
	# There should only be one dup
	if (nrow(dup.positions) != 1)
		stop('There should only be one and only one set of dup positions.')
	minn <- min(gene.positions[, 'start'])
	maxx <- max(gene.positions[, 'end'])
	full.length <- maxx - minn
	gr <- par('usr')
	#lines(c(minn - (full.length / 40), maxx + (full.length / 40)), c(y, y), lwd = 2, col = 'grey20')
	lines(c(gr[1], gr[2]), c(y, y), lwd = 2, col = 'grey20')
	v.adj = (gr[4] - gr[3])/20
	h.adj = (gr[2] - gr[1])/50
	# Draw the genes
	ace1.col <- 'tomato'
	other.col <- 'grey60'
	gene.colors <- c(other.col, ace1.col)[(gene.positions$Name == 'Ace1') + 1]
	rect(as.matrix(gene.positions[,'start']), y - v.adj, as.matrix(gene.positions[,'end']), y + v.adj, col = gene.colors, border = NA)
	#text(apply(gene.positions[,c('start', 'end')], 1, mean), y + 3*v.adj, gene.positions$Name, adj = 0.5)
	# Add the dup
	d.v.adj = (gr[4] - gr[3])/30
	y.pos <- y - (5*v.adj)
	rect(as.matrix(dup.positions[,'start']), y.pos - d.v.adj, as.matrix(dup.positions[,'end']), y.pos + d.v.adj, col = c('deepskyblue4'), border = NA)
	text(as.matrix(dup.positions[,'start']) - h.adj, y.pos, dup.positions$Dup.name, adj = 1, cex = text.cex)
	# Add any dels
	if ((!is.null(del.positions)) & (nrow(del.positions) > 0)){
		y.pos <- y.pos - (5*d.v.adj) * 1:nrow(del.positions)
		rect(as.matrix(del.positions[,'start']), y.pos - d.v.adj, as.matrix(del.positions[,'end']), y.pos + d.v.adj, col = c('burlywood3'), border = NA)
		text(as.matrix(del.positions[,'start']) - h.adj, y.pos , del.positions$Dup.name, adj = 1, cex = text.cex)
	}
}

plot.significance.along.ace1 <- function(significance.value, significance.direction, snp.pos, filename, chrom.col = NULL, lwd = 0.3, ...){
	# Sanity checks
	if (length(significance.value) != nrow(snp.pos))
		stop('Number of points in significance.value is different to rows in snp.pos.')
	if (length(significance.value) != length(significance.direction))
		stop('Number of points in significance.value is different to significance.direction.')
	if (any(snp.pos$Chrom != '2R'))
		stop('Some snps are not on Chromosome 2R.')
	# Make sure the SNPs are in order
	snp.order <- order(snp.pos$Pos)
	snp.pos <- snp.pos[snp.order, ]
	significance.value <- significance.value[snp.order]
	significance.direction <- significance.direction[snp.order]
	# Set up the plotting output file
	if (!missing(filename)){
		file.width = 6.5
		file.height = 5
		if (grepl('\\.eps', filename))
			postscript(filename, width = file.width, height = file.height, horizontal = F, onefile = FALSE, paper = "special")
		else if (grepl('\\.png', filename))
			png(filename, width = file.width, height = file.height, units = 'in', res = 600)
		else if (grepl('\\.tif', filename))
			tiff(filename, width = file.width, height = file.height, units = 'in', res = 600)
	}
	lhei <- c(2,1.2)
	lmat <- matrix(1:2, 2, 1)
	# Color-code the points according to whether it is in the Dup, in a Del, or in Ace1 gene. And shape code whether
	# its a positive or negative relationship
	ace1.col <- 'tomato'
	dup1.col <- 'deepskyblue4'
	del.col <- 'burlywood3'
	snp.colours <- rep('black', nrow(snp.pos))
	snp.colours[snp.pos$Pos >= as.numeric(Dup.coords['Dup1', 'start']) & snp.pos$Pos < as.numeric(Dup.coords['Dup1', 'end'])] <- dup1.col
	snp.colours[snp.pos$Pos >= min(Dup.coords['Del1', 'start']) & snp.pos$Pos < max(Dup.coords['Del1', 'end'])] <- del.col
	snp.colours[snp.pos$Pos >= as.numeric(Dup.gene.positions[ace1.gene, 'start']) & snp.pos$Pos < as.numeric(Dup.gene.positions[ace1.gene, 'end'])] <- ace1.col
	snp.pch <- c(20, 17)[(significance.direction == 1)+1]
	# Create the layout
	layout(lmat, heights = lhei)
    par(mar = c(0, 3, 1, 1), xpd = NA, mgp = c(2, 0.6, 0), tcl = -0.3, ...)
	plot(-log10(significance.value), bty = 'n', xlab = '', xaxt = 'n', cex = 0.7, ylab = '-log(P)', col = snp.colours, pch = snp.pch)
	coords.1 <- par('usr')	
	xrange.1 <- coords.1[2] - coords.1[1]
	# Plot the Dup region
	top.margin <- 2
	par(mar = c(2.5, 3, top.margin, 1), xpd = NA, mgp = c(1.3, 0.3, 0), ...)
	plot(c(min(snp.pos$Pos), max(snp.pos$Pos)), c(0, 1.6), type = 'n', bty = 'n', xlab = 'Position on 2R', ylab = '', yaxt = 'n', cex.axis = 0.8, tcl = -0.25, cex.lab = 0.8)
	draw.gene.model(Dup.gene.positions, Dup.coords['Dup1', ], Dup.coords[paste('Del', 1:3, sep = ''), ], 1.5)
	coords.2 <- par('usr')
	xrange.2 <- coords.2[2] - coords.2[1]
	# Write a function to convert old x coords to new x coords
	f <- 1
	convert <- function(x) (((x - coords.1[1])/xrange.1) * xrange.2 + coords.2[1])/f
	# And work out the top end of the current plot
	top.end <- par('cxy')[2] * par('cex') * par('lheight') * top.margin + coords.2[4]
	# For every snp, draw a line from the plot to the chromosome.
	chrom.col <- if (!is.null(chrom.col)) chrom.col else c('2R' = 'black', '2L' = 'black', '3R' = 'black', '3L' = 'black', 'X' = 'black')
	segments(snp.pos$Pos, 1.62, convert(1:length(significance.value)), top.end, lwd = lwd, col = snp.colours)
	if (!missing(filename))
		dev.off()
}

#plot.significance.along.ace1(ace1.snps$logregP, ace1.snps$logreg.direction, ace1.snps[, c('Chrom', 'Pos')], family = 'Arial')
plot.significance.along.ace1(ace1.snps$logregP, ace1.snps$logreg.direction, ace1.snps[, c('Chrom', 'Pos')], filename = 'sig_along_ace1.png', family = 'Arial')

# Let's load up the Ace1 Dup calls. 
CNV.diagnostic.read.calls.fn <- '~/data/CNV_GAARD/target_regions_CNVs.csv'
cnvs.allsamples <- fread(CNV.diagnostic.read.calls.fn)
colnames(cnvs.allsamples)[1] <- 'Sample.name'
setkey(cnvs.allsamples, Sample.name)
ace1.Dup1 <- cnvs.allsamples[sample.names, 'Ace1_Dup1']
ace1.Del1 <- cnvs.allsamples[sample.names, 'Ace1_Del1']

# For each SNP, get the correlation coefficient with Dup1 and Del1
snp.cor <- function(x, y){
	x[phenotypes == 'alive'] = x[phenotypes == 'alive'] - mean(x[phenotypes == 'alive'])
	x[phenotypes == 'dead'] = x[phenotypes == 'dead'] - mean(x[phenotypes == 'dead'])
	y = unlist(y)
	y[phenotypes == 'alive'] = y[phenotypes == 'alive'] - mean(y[phenotypes == 'alive'])
	y[phenotypes == 'dead'] = y[phenotypes == 'dead'] - mean(y[phenotypes == 'dead'])
	cor(x, y)
}

top.snps[, Dup1.cor :=  apply(.SD, 1, snp.cor, ace1.Dup1), .SDcols = sample.names]
top.snps[, Del1.cor :=  apply(.SD, 1, snp.cor, ace1.Del1), .SDcols = sample.names]
ace1.snps[, Dup1.cor :=  apply(.SD, 1, snp.cor, ace1.Dup1), .SDcols = sample.names]
ace1.snps[, Del1.cor :=  apply(.SD, 1, snp.cor, ace1.Del1), .SDcols = sample.names]
non.ace1.snps[, Dup1.cor :=  apply(.SD, 1, snp.cor, ace1.Dup1), .SDcols = sample.names]
non.ace1.snps[, Del1.cor :=  apply(.SD, 1, snp.cor, ace1.Del1), .SDcols = sample.names]

#plot.significance.along.ace1(10^(-abs(ace1.snps$Dup1.cor)), ace1.snps$logreg.direction, ace1.snps[, c('Chrom', 'Pos')])
#plot.significance.along.genome(10^(-abs(top.snps$Dup1.cor)), top.snps[, c('Chrom', 'Pos')], chrom.col = chrom.colours, lwd = 0.3)
#plot.significance.along.ace1(10^(-abs(ace1.snps$Del1.cor)), ace1.snps$logreg.direction, ace1.snps[, c('Chrom', 'Pos')])
#plot.significance.along.genome(10^(-abs(top.snps$Del1.cor)), top.snps[, c('Chrom', 'Pos')], chrom.col = chrom.colours, lwd = 0.3)
# Some of the non-Ace1 SNPs have quite high correlation with Ace1, so we will need to be a bit careful 
# with those. 

# Let's plot the correlation with Ace1 against P-value
png('Dup1_correlation_vs_significance.png', width = 8, height = 5, units = 'in', res = 300)
par(mfrow = c(1,2), family = 'Arial')
plot(abs(top.snps$Dup1.cor), -log(top.snps$logregP), pch = 21, col = 'grey50', bg = chrom.colours[top.snps$Chrom])
plot(abs(non.ace1.snps$Dup1.cor), -log(non.ace1.snps$logregP), pch = 21, col = 'grey50', bg = chrom.colours[non.ace1.snps$Chrom])
dev.off()
# There aren't any majorly significant SNPs that aren't correlated with Ace1. 

save.image('classical_analysis_nosibs.Rdata')

png('manhattan_plot_nosibs.png', width = 6, height = 6, res = 300, units = 'in')
par(family = 'Arial')
manhattan(snp.table[, .(numChrom, Pos, logregP, ID)], chr = 'numChrom', bp = 'Pos', p = 'logregP', snp = 'ID', suggestiveline = F, genomewideline = -log10(0.01), chrlabs = chrom.order, col = c('blue4', 'orange3'), ylab = '-log10(P)', annotatePval = 0.05, main = 'Logistic regression P, no sibs')
dev.off()



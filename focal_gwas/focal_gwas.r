library(data.table)
library(RColorBrewer)
library(lme4)
library(stringr)
library(Biostrings)
library(magrittr)
library(ape)

arg.values <- commandArgs(trailingOnly=T)
cat('Arguments:', arg.values, '\n', sep = '\n')

study.pop <- arg.values[1]
filter.regions <- arg.values[-1]

# Allow the common filtered regions to be defined by name
filter.regions[filter.regions == 'Ace1'] <- '2R:3100000-3900000'
filter.regions[filter.regions == '2La'] <- '2L:20524000-42165600'

chrom.sizes <- c('2R' = 61545105, '2L' = 49364325, '3R' = 53200684, '3L' = 41963435, 'X' = 24393108)

dir.create(study.pop)

# Write a function to filter out a region from a "filtered.tables" entry
remove.region <- function(filtered.table, regions){
	chroms <- str_extract(regions, '^[23]?[LRX](?=:)')
	start <- str_extract(regions, '(?<=:)\\d+(?=-)')
	end <- str_extract(regions, '(?<=-)\\d+$')
	for (i in 1:length(chroms)){
		filtered.table <- filtered.table[window.chrom != chroms[i] | 
										 window.pos > end[i] |
										 window.pos < start[i]]
	}
	filtered.table
}


# We first need to load the full SNP table to work out which SNPs belong in which window
snp.table.fn <- '~/data/ML/GAARD_SNP/filtered_snp_tables/filtered_snp_tables.Rdata'
load(snp.table.fn)

snp.table <- snp.tables[[study.pop]][[1]][, .(Chrom, Pos)]
rm(snp.tables)

# Load the results of Fst-based window selection
load('../randomisations/Fst/fst_filtered_windows.Rdata')

filtered.table <- filtered.tables[[study.pop]]
rm(filtered.tables)

# Now remove the windows inside 2La, since there are so many that we wouldn't know where to start
# Inversion coordinates taken from Sharakov et al 2006 (PNAS)
for (region in filter.regions)
	filtered.table <- remove.region(filtered.table, region)

# Get the number of windows that will be found on each chromosome
chrom.snp.counts <- tapply(snp.table$Chrom, snp.table$Chrom, length)
window.num <- ceiling(chrom.snp.counts / window.size)
window.end <- cumsum(window.num)
window.start <- window.end - window.num + 1
snp.table$window <- unlist(sapply(names(window.num), function(chrom) rep(window.start[chrom]:window.end[chrom], each = window.size, length.out = chrom.snp.counts[chrom])))
window.pos <- snp.table[, .(window.chrom = unique(Chrom), window.pos = floor(median(Pos))), by = window][, window := NULL]
window.pos$window <- 1:nrow(window.pos)

# Now find out which SNPs are found in each of the significant windows
get.window.start.end <- function(chrom, pos, window.pos){
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

focal.window.ranges <- mapply(get.window.start.end,
                              filtered.table$window.chrom, 
                              filtered.table$window.pos, 
                              MoreArgs = list(window.pos = window.pos),
                              SIMPLIFY = F
                             ) %>%
                             rbindlist() %>%
							 setkey(name)

# Load the sib groups information
sib.groups <- fread('../NGSrelate/full_relatedness/sib_group_table.csv', sep = '\t')
sib.groups[, pop := paste(location, species, insecticide, sep = '_')]
samples.to.remove <- sib.groups[keep == F & pop == study.pop, sample.name]

# We now need to load the SNP data including non-accessible sites. 
# Load SNP data
SNP.filename <- paste(
	'~/data/ML/GAARD_SNP', 
	study.pop, 
	'csvs/filtered_SNPs.csv', 
	sep = '/'
)
cat('Loading SNP data from ', SNP.filename, '.\n', sep = '')
all.snps <- fread(SNP.filename)
# Remove sibs
all.snps[, c(samples.to.remove) := NULL]
sample.names <- grep('WA-\\d{4}', colnames(all.snps), value = T)
# Focus on focal regions
extract.window.snps <- function(window.name, region.border.size = 10000){
	all.snps[Chrom == focal.window.ranges[window.name, chrom] & 
	         Pos >= (focal.window.ranges[window.name, start] - region.border.size) &
	         Pos <= (focal.window.ranges[window.name, end] + region.border.size)] %>%
	.[, Window.name := window.name] %>%
	# Record whether the SNP is in the focal window, on in the border
	.[, In.window := Pos >= focal.window.ranges[window.name, start] & Pos <= focal.window.ranges[window.name, end]]
}

# The threshold MAF that we will use to filter SNPs
MAF.thresh <- 5
# Function to calculate minor allele frequency
maf <- function(genotype){
	half.hap.num <- length(genotype)
	half.hap.num - abs(sum(genotype) - half.hap.num)
}

focal.snps <- lapply(focal.window.ranges$name, extract.window.snps) %>%
              rbindlist() %>%
			  # Remove SNPs with MAF smaller than threshold
              .[, MAF := apply(.SD, 1, maf), .SDcols = sample.names] %>%
              .[MAF >= MAF.thresh, ]


rm(all.snps)
gc()

# Load phenotype data
phenotype.filename <- '../data/combined/sample_phenotypes.csv'
cat('Loading phenotype data from ', phenotype.filename, '.\n', sep = '')
# The object sample.names already excludes dicarded sibs, so applying here makes the phenotypes
# match the samples from the snp table, as well as their order. 
phenotype.table <- fread(phenotype.filename, key = 'specimen')[sample.names, ]
phenotypes <- setNames(as.factor(phenotype.table[sample.names, phenotype]),
                       nm = sample.names)

# Logistic regression
logreg.test <- function(genotype, phenotype){
	genotype <- as.numeric(genotype)
	null.model <- glm(phenotype ~ 1, family = 'binomial')
	model1 <- glm(phenotype ~ genotype, family = 'binomial')
	pval <- anova(model1, null.model, test = 'Chisq')[['Pr(>Chi)']][2]
	direction <- sign(model1$coefficients['genotype'])
	data.table(logregP = pval, direction = direction)
}

significance.test <- apply(focal.snps[, ..sample.names], 
                           1, 
                           logreg.test,
                           phenotype = phenotypes
                     ) %>%
                     rbindlist %>%
                     cbind(focal.snps[, .(Chrom, Pos, Window.name, In.window)], .)

# Set the NA P values to 1
significance.test$logregP[is.na(significance.test$logregP)] <- 1

# Keep only the significant SNPs with p < 0.01 and turn them into a VCF
p.thresh <- 0.01
sig.snps <- significance.test[logregP < p.thresh, ]

sample.names.string <- paste(sample.names, collapse = ',')
study.id.table <- fread('../data/study_ids.csv', key = 'population')
study.id <- study.id.table[study.pop, study_id]

get.allele.counts <- function(window.name, region.border.size = 10000){
	chrom <- focal.window.ranges[window.name, chrom]
	window.start <- max(0, focal.window.ranges[window.name, start] - region.border.size)
	window.end <- min(chrom.sizes[chrom], focal.window.ranges[window.name, end] + region.border.size)
	window.range <- paste(window.start, window.end, sep = '-')
	allele.counts.filename <- paste(study.pop, '/', window.name, '_allele_counts.csv', sep = '')
	system(paste('python extract_zarr_region.py ~/scratch/VObs_GAARD',
	                     study.id,
	                     chrom,
	                     window.range, 
	                     allele.counts.filename,
	                     sample.names.string, 
	                     sep = ' '))
	allele.counts <- fread(allele.counts.filename, sep = '\t')
	system(paste('rm', allele.counts.filename))
	colnames(allele.counts)[1] <- 'pos'
	setkey(allele.counts, pos)
	allele.counts
}

allele.counts.tables <- lapply(setNames(nm = focal.window.ranges$name), get.allele.counts)

# Create a vcf for SNP Eff
sig.snps.ac <- sig.snps[, ..allele.counts.tables[[unique(Window.name)]][J(Pos)], by = Window.name]
alt.counts.columns <- paste('counts', 1:3, sep = '')
alt.allele.columns <- paste('allele', 1:3, sep = '')
# as.numeric is needed because by taking the whole of the table's row as a vector, the numbers get turned into 
# characters
sig.snps.ac[, chrom := sub(':.*', '', Window.name)]
sig.snps.ac[, alts := apply(.SD, 1, function(x) paste(x[..alt.allele.columns][as.numeric(x[..alt.counts.columns]) > 0], collapse = ','))]
sig.snps.ac[, major.alt := apply(.SD, 1, function(x) paste(x[..alt.allele.columns][which.max(x[..alt.counts.columns])]))]
sig.snps.ac[, info := ifelse(apply(.SD, 1, function(x) sum(x > 0)) > 1, 'Multiallelic', ''), .SDcols = alt.counts.columns]
vcf <- sig.snps.ac[, .(chrom, pos, '.', ref = allele0, alt = alts, '.', '.', info)]
vcf.file <- paste(study.pop, '/temp.vcf', sep = '')
vcf.output.file <- paste(study.pop, '/focal_windows.vcf', sep = '')
write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO', file = vcf.file)
fwrite(vcf, vcf.file, sep = '\t', append = T, quote = F)

# Now run SNPEff on that. 
snp.eff.command <- paste(
	'java -jar ~/software/snpEff/snpEff.jar eff Anopheles_gambiae', 
	'-noStats', 
	vcf.file, 
	'>', 
	vcf.output.file
)
cat('Running snpEff on significant SNPs, using command', snp.eff.command, '.\n')
system(snp.eff.command)
system(paste('rm ', study.pop, '/temp.vcf', sep = ''))

# Now load up the new vcf table with the snpEff annotations.
annotated.vcf <- fread(vcf.output.file)

# Join all of the desired information into a single table
sig.snps.annotated <- cbind(
	sig.snps, 
	sig.snps.ac[, .(ref = allele0, alts, major.alt)], 
	annotated.vcf[, .(snpEff = INFO)]
)
sig.snps.annotated[, genes := sapply(snpEff, function(str) paste(unique(unlist(str_match_all(str, 'AGAP\\d{6}'))), sep = ','))]

# Write a function to pull out all the effect types from a given INFO field
pull.out.effect.types <- function(info){
	all.infos <- strsplit(info, ',')[[1]]
	all.effect.types <- all.infos %>%
	                    sapply(function(x) strsplit(x, '\\|')[[1]][[2]]) %>%
	                    paste(collapse = ',')
	all.effect.types
}
sig.snps.annotated[, effects := sapply(snpEff, pull.out.effect.types)]

# We will need some code that plots a given focal window, along with the (labelled) genes contained 
# therein.
gff.path <- '~/data/ML/GAARD_SNP/data/Anopheles-gambiae-PEST_BASEFEATURES_AgamP4.9.gff3'
gff <- as.data.table(read.gff(gff.path, GFF3 = T))

lighten.col <- function(color, lightness, alpha = alpha){
	col.rgb <- col2rgb(color)/255
	rgb(t(1-(1-col.rgb)*lightness), alpha = alpha)
}

draw.gene.model <- function(gene.positions, exon.positions, y = 0, text.cex = 0.5){
	gr <- par('usr')
	lwd = 2
	lines(c(gr[1], gr[2]), c(y, y), lwd = lwd, col = 'grey20')
	if (nrow(gene.positions) > 0){
		gene.positions <- gene.positions[order(start)]
		gene.thickness <- (gr[4] - gr[3])/10
		# Draw the genes
		v.adj = ifelse(gene.positions$strand == '+', 0, -gene.thickness)
		rect(
			gene.positions[,start], 
			y + v.adj, 
			gene.positions[,end], 
			y + gene.thickness + v.adj, 
			col = 'grey95',
			border = 'black',
			lwd = lwd,
			xpd = F
		)
		# Draw the exons
		exon.v.adj = ifelse(exon.positions$strand == '+', 0, -gene.thickness)
		rect(
			exon.positions[,start], 
			y + exon.v.adj, 
			exon.positions[,end], 
			y + gene.thickness + exon.v.adj, 
			col = 'black',
			border = NA,
			xpd = F
		)
		# Label the genes
		text(
			apply(gene.positions[,.(start, end)], 1, mean), 
			y - 1.6*gene.thickness, 
			gene.positions$gene.name, 
			adj = 1, 
			srt = 45,
			cex = text.cex
		)
	}
}

plot.focal.window <- function(window.name, filename, sig.snps.only = T, true.pos = T, region.border.size = 10000, segment.lwd = 0.2, ...){
	chrom <- focal.window.ranges[window.name, chrom]
	start.pos <- max(0, focal.window.ranges[window.name, start] - region.border.size)
	end.pos <- min(chrom.sizes[chrom], focal.window.ranges[window.name, end] + region.border.size)
	
	# Get the list of genes present in the window
	gene.gff <- gff[seqid == chrom & 
	                start < end.pos & 
	                end > start.pos &
	                type == 'gene', ]
	
	gene.gff[, gene.id := unique(str_extract(attributes, 'AGAP\\d{6}'))]
	gene.gff[, gene.name := str_extract(attributes, '(?<=Name=)[^;]+') %>%
	                        str_to_title]
	gene.gff$gene.name[is.na(gene.gff$gene.name)] <- gene.gff$gene.id[is.na(gene.gff$gene.name)]
	# Now the list of Exons
	exon.gff <- gff[seqid == chrom & 
	                start < end.pos & 
	                end > start.pos &
	                type == 'exon', ]
	
	# Set up the plotting output file
	if (!missing(filename)){
		file.width = 6.5
		file.height = 5
		if (grepl('\\.eps', filename))
			postscript(filename, width = file.width, height = file.height, horizontal = F, onefile = FALSE, paper = "special")
		else if (grepl('\\.png', filename))
			png(filename, width = file.width, height = file.height, units = 'in', res = 200)
		else if (grepl('\\.tif', filename))
			tiff(filename, width = file.width, height = file.height, units = 'in', res = 200)
	}
	
	lhei <- c(2,0.7)
	lmat <- matrix(1:2, 2, 1)
	if (sig.snps.only){
		these.snps <- sig.snps.annotated[Window.name == window.name][order(Pos), ]
		snp.colours <- ifelse(grepl('missense', these.snps$effects), 'red', 'cornflowerblue')
	}
	else {
		these.snps <- significance.test[Window.name == window.name][order(Pos), ]
		snp.colours <- ifelse(these.snps$logregP < p.thresh, 'red', 'cornflowerblue')
	}
	
	snp.colours[!these.snps$In.window] <- lighten.col(snp.colours[!these.snps$In.window], 0.5)
	
	snp.pch <- c(21, 24)[(these.snps$direction == 1)+1]
	# Create the layout
	layout(lmat, heights = lhei)
    par(mar = c(0, 3, 1, 1), xpd = NA, mgp = c(2, 0.6, 0), tcl = -0.3, ...)
	if (true.pos){
		plot(
			these.snps$Pos,
			-log10(these.snps$logregP), 
			bty = 'n', 
			xlim = c(start.pos, end.pos),
			xlab = '', ylab = '-log(P)', xaxt = 'n', 
			cex = 0.9, col = 'grey50', bg = snp.colours, 
			pch = snp.pch, main = window.name
		)
		abline(v = c(focal.window.ranges[window.name, start], 
		             focal.window.ranges[window.name, end]), 
		       col = lighten.col('grey60', 1, 0.5),
			   xpd = F
		)
	}
	else {
		plot(
			-log10(these.snps$logregP), 
			bty = 'n', 
			xlab = '', ylab = '-log(P)', xaxt = 'n', 
			cex = 0.7, col = snp.colours,
			pch = snp.pch, main = window.name
		)
		abline.1 <- min(which(these.snps$In.window)) - 0.5
		abline.2 <- max(which(these.snps$In.window)) + 0.5
		abline(v = c(abline.1, abline.2), col = lighten.col('grey60', 1, 0.5), xpd = F)
	}
	coords.1 <- par('usr')	
	xrange.1 <- coords.1[2] - coords.1[1]
	# Plot the genes
	top.margin <- 0.8
	par(mar = c(2.5, 3, top.margin, 1), xpd = NA, mgp = c(1.3, 0.3, 0), ...)
	plot(
		c(start.pos, end.pos), 
		c(0, 0.6), 
		type = 'n', bty = 'n', 
		xlab = paste('Position on', chrom), ylab = '', yaxt = 'n', 
		cex.axis = 0.8, tcl = -0.25, cex.lab = 0.8
	)
	draw.gene.model(gene.gff, exon.gff, 0.5)
	if (!true.pos){
		coords.2 <- par('usr')
		xrange.2 <- coords.2[2] - coords.2[1]
		# Write a function to convert old x coords to new x coords
		f <- 1
		convert <- function(x) (((x - coords.1[1])/xrange.1) * xrange.2 + coords.2[1])/f
		# And work out the top end of the current plot
		top.end <- par('cxy')[2] * par('cex') * par('lheight') * top.margin + coords.2[4]
		# For every snp, draw a line from the plot to the chromosome.
		segments(these.snps$Pos, 0.62, convert(1:nrow(these.snps)), top.end, lwd = segment.lwd, col = snp.colours)
	}
	if (!missing(filename))
		dev.off()
}

mapply(
	plot.focal.window, 
	focal.window.ranges$name, 
	# We replace ':' with '_' in the filename, because bash seems to struggle with "X:"
	paste(
		study.pop, 
		'/', 
		sub(':', '_', focal.window.ranges$name), 
		'.png', 
		sep = ''
	),
	MoreArgs = list(family = 'Arial')
)

system(paste('rm', vcf.output.file))
fwrite(focal.window.ranges, file = paste(study.pop, '/focal_window_ranges.csv', sep = ''), sep = '\t')
fwrite(sig.snps.annotated, file = paste(study.pop, '/sig_snps_annotated.csv', sep = ''), sep = '\t')


library(data.table)
library(stringr)
library(fdrtool)
library(magrittr)

# Load the data
locations <- c('Avrankou', 'Baguida', 'Korle-Bu', 'Madina', 'Obuasi')

fst.tables <- list()
for (l in locations){
	cat('Loading data for', l, '\n')
	load(paste('fst_randomisations_', l, '.Rdata', sep = ''))
	for (pop in names(fst.table))
		fst.tables[[pop]] <- cbind(window.pos[[pop]], fst.table[[pop]])
}

phenotype.column <- 'phenotype.moving.Fst'
randomisation.columns <- grep('^r\\d{5}\\.', colnames(fst.tables[[1]]), value = T)
populations <- names(fst.tables)

# Calculate a P-value for each window based on the number of randomisations larger than the observed value
for (pop in populations){
	fst.tables[[pop]][, window.p := rowSums(apply(.SD, 2, function(x) x >= get(phenotype.column))) / ncol(.SD),
	                    .SDcols = randomisation.columns]
	fst.tables[[pop]][, log.p := -log10(window.p + 1/(length(..randomisation.columns)))]
	fst.tables[[pop]][, fdr := fdrtool(window.p, 'pvalue', plot = F)$qval]
}


# Write a function to add transparency to a colour
add.transparency <- function(col, prop.alpha){
	rgb.val <- col2rgb(col)
	new.rgb <- rgb(rgb.val[1], rgb.val[2], rgb.val[3], max = 255, alpha = (1 - prop.alpha) * 255)
	new.rgb
}


chrom.sizes <- c('2R' = 61545105, '2L' = 49364325, '3R' = 53200684, '3L' = 41963435, 'X' = 24393108)

# Write a function to draw the chromosomes on an existing plotting device
add.chromosomes <- function(gaps, cs, ce, gene.cex = 0.9, gene.col = 'grey20', point.cex = 1.2, point.col = 'grey30', chrom.col = NULL, chrom.cex = 1.4, chrom.offset = 0){
	plot(c(cs[1], ce[5]), c(-6.5,1.3), xlim = c(cs[1] + gaps/2, ce[5] - gaps/2), type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
	# show the kdr region
	kdr.region.mean <- cs['2L'] + mean(c(2358158, 2431617))
	points(kdr.region.mean, 0, pch = 19, cex = point.cex, col = point.col)
	text(kdr.region.mean, -1.7, 'Vgsc', srt = 45, adj = 1.1, xpd = T, cex = gene.cex, col = gene.col, font = 2)
	# show the Rdl region
	rdl.region.mean <- cs['2L'] + mean(c(25363652, 25434556))
	points(rdl.region.mean, 0, pch = 19, cex = point.cex, col = point.col)
	text(rdl.region.mean, -1.7, 'Rdl', srt = 45, adj = 1.1, xpd = T, cex = gene.cex, col = gene.col, font = 2)
	# show the Carboxylesterase region
	coeae.region.mean <- cs['2L'] + mean(c(28548433, 28550748))
	points(coeae.region.mean, 0, pch = 19, cex = point.cex, col = point.col)
	text(coeae.region.mean, -1.7, 'Coeae2f', srt = 45, adj = 1.1, xpd = T, cex = gene.cex, col = gene.col, font = 2)
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


# Function to plot windowed Fst P values 
plot.randomised.fst <- function(fst.table, filename = NULL, num.randomisations = NULL, fdr.cutoff = 0.01, plot.title = '', gaps = 5000000){
	# get vectors of start and end points for each chromosome (ie: cumulative sizes + gaps)
	ce <- cumsum(chrom.sizes + c(0, 0, gaps, 0, gaps))
	cs <- ce - chrom.sizes
	# Create the plot
	if (!missing(filename)){
		file.width = 6.5
		file.height = 3.5
		if (grepl('\\.eps', filename))
			postscript(filename, width = file.width, height = file.height, horizontal = F, onefile = FALSE, paper = "special")
		else if (grepl('\\.png', filename))
			png(filename, width = file.width, height = file.height, units = 'in', res = 600)
		else if (grepl('\\.tif', filename))
			tiff(filename, width = file.width, height = file.height, units = 'in', res = 600)
	}
	layout(matrix(c(rep(1,4), rep(2,4), rep(3,2)), nrow = 10, ncol = 1))
	colours <- c(fst = add.transparency('orangered3', 0.4),
                 randomisations = 'grey50',
                 p = 'blue')
	par(mar = c(0,4,1,2), mgp = c(2, 0.7, 0), family = 'Arial', xpd = NA) 
	num.randomisations <- ifelse(is.null(num.randomisations), 
	                             length(randomisation.columns),
	                             num.randomisations)
	# This odd way of getting a sequence is to make sure we get the right outcome if num.randomisations == 0
	r.columns <- randomisation.columns[seq(1, num.randomisations, length.out = num.randomisations)]
	fst.columns <- c(phenotype.column, r.columns)
	max.y <- max(c(max(fst.table[, ..fst.columns]), 0.05))
	min.y <- min(fst.table[, ..fst.columns])
	# Create the empty plot.
	plot(c(cs[1], ce[5]), c(min.y, max.y), xlim = c(cs[1] + gaps/2, ce[5] - gaps/2), type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', main = plot.title, cex.main = 1.2)
	
	# Get plotting positions for windows
	fst.table$genome.pos <- fst.table$window.pos + cs[as.character(fst.table$window.chrom)]
	# Add randomised data
	sapply(r.columns, function(x) fst.table[, lines(genome.pos, get(x), col = colours['randomisations'], lwd = 0.4), by = window.chrom])
	# Add true data
	fst.table[, lines(genome.pos, phenotype.moving.Fst, col = colours['fst'], lwd = 1.2), by = window.chrom]
	# Add y axis
	fst.step.size <- ifelse(max.y > 0.2, 0.1, 0.05)
	axis(2, at = seq(0, max.y, fst.step.size))
	mtext('Fst', 2, 2, cex = 0.8)
	
	# Plot the P value data
	par(mar = c(0,4,1,2), mgp = c(2, 0.7, 0), family = 'Arial', xpd = NA) 
#	max.p <- max(fst.table$log.p)
#	min.p <- min(fst.table$log.p)
	# Instead of observing max.p and min.p, let's just use their theoretical extremes. 
	min.p <- 0
	max.p <- -log10(1/length(randomisation.columns))
	plot(c(cs[1], ce[5]), c(min.p, max.p), xlim = c(cs[1] + gaps/2, ce[5] - gaps/2), type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
	fst.table[, lines(genome.pos, log.p, col = colours['p'], lwd = 1.2), by = window.chrom]
	fst.table[fdr <= get('fdr.cutoff'), points(genome.pos, rep(..max.p, nrow(.SD)), col = 2, pch = 19, cex = 0.5)]
	p.cutoff.table <- unique(fst.table[, .(log.p, fdr)])
	# We define the p-cutoff and slightly more than the largest -log10(p) that was non-sig (that way we have a value
	# even for the datasets where nothing was significant). 
	p.cutoff <- max(p.cutoff.table[fdr > get('fdr.cutoff'), log.p]) + 0.05
	abline(h = p.cutoff, col = 2, xpd = F)
	# Add y axis
	p.step.size <- ifelse(max.p > 2, 1, 0.5)
	axis(2, at = seq(0, max.p, p.step.size))
	mtext('-log10(P)', 2, 2, cex = 0.8)
	
	# Now plot all chromosomes with, the position of each of the four detox gene regions and Ace1
	par(mar = c(1,4,0,2), mgp = c(2, 0.7, 0)) 
	add.chromosomes(gaps = gaps, cs = cs, ce = ce, gene.cex = 0.7, point.cex = 1, chrom.offset = -1.2, chrom.cex = 1.2)
	if (!missing(filename))
		dev.off()
}

save.image('Fst_p_values.Rdata')

# Draw all the plots
for (pop in names(fst.tables)){
	fdr.cutoff <- ifelse(grepl('PM', pop), 0.01, 0.05)
	fn <- paste(pop, '_fdr_plot.png', sep = '')
	title <- paste(gsub('_', ' ', pop), 'FDR cutoff =', fdr.cutoff)
	plot.randomised.fst(fst.tables[[pop]], filename = fn, plot.title = title, num.randomisations = 1000)
}

# Instead of calculating FDR, let's find some arbitrary Fst cut-off to do a fist pass control on the peaks, 
# and then use raw P-values for those peaks. 
# We could try a cutoff of Fst = 0.02
# We can also try a cutoff of Fst = 0.03
# Or we can try something like "take the 99th centile of the data and double its distance from the 50th centile. 
# In the end, after testing the effectiveness of the different cutoffs on the various datasets, I settled for 
# using the negative part of the distribution (to the left of the mode) to determine was the positive part of the
# distribution "should" look like, and used this to find extremems. 

# Use the negative part of the distribution to tell us how far to look in the positive part. 
first.pass.peaks <- function(fst.vector, return.thresh = F){
	fst.modal.bin <- which.max(table(cut(fst.vector, breaks = seq(min(fst.vector), max(fst.vector), 0.0005))))
	fst.mode <- mean(as.numeric(strsplit(sub(']$', '', sub('^\\(', '', names(fst.modal.bin))), ',')[[1]]))
	neg.diff <- fst.mode - min(fst.vector)
	thresh <- fst.mode + neg.diff*2.5
	if (return.thresh)
		thresh
	else
		fst.vector > thresh
}

for (pop in names(fst.tables))
	fst.tables[[pop]][, filter := first.pass.peaks(phenotype.moving.Fst)]

# Create filters that join all the filter-passing windows for all insecticidies. Since the different populations
# have different window positions, we need to find a way to get around this. We might need to do nearest windows.

find.nearest.windows <- function(pos.1, pos.2){
	distances <- abs(outer(pos.1, pos.2, '-'))
	nearest.windows <- apply(distances, 1, which.min)
	nearest.windows
}

# First, get a combined list of all filter passing windows
delta.pops <- grep('Delta', names(fst.tables), value = T)
delta.all.filtered.windows <- fst.tables[delta.pops] %>%
                              lapply(function(D) D[filter == T, .(window.chrom, window.pos)]) %>%
                              rbindlist %>%
                              split(.$window.chrom)
for (pop in delta.pops){
	# Create a new filter called joined.filter, which is true for every window that is the nearest window to one
	# of the filter-passing windows from any of the deltamethrin datasets. 
	fst.tables[[pop]][, joined.filter := ..delta.all.filtered.windows[[unique(window.chrom)]]$window.pos %>%
	                                     find.nearest.windows(window.pos) %>%
	                                     replace(logical(nrow(.SD)), ., T), 
	                    by = window.chrom]
}

# Do the same for PM
pm.pops <- grep('PM', names(fst.tables), value = T)
pm.all.filtered.windows <- fst.tables[pm.pops] %>%
                              lapply(function(D) D[filter == T, .(window.chrom, window.pos)]) %>%
                              rbindlist %>%
                              split(.$window.chrom)
for (pop in pm.pops){
	# Create a new filter called joined.filter, which is true for every window that is the nearest window to one
	# of the filter-passing windows from any of the pmmethrin datasets. 
	fst.tables[[pop]][, joined.filter := ..pm.all.filtered.windows[[unique(window.chrom)]]$window.pos %>%
	                                     find.nearest.windows(window.pos) %>%
	                                     replace(logical(nrow(.SD)), ., T), 
	                    by = window.chrom]
}

# Show how these filters affect the peaks
plot.fst.filters <- function(fst.table, filename = NULL, plot.title = '', gaps = 5000000, p.thresh = c(0.01, 0.001), filter.name = 'filter'){
	# Create the plot
	if (!missing(filename)){
		file.width = 6.5
		file.height = 3.5
		if (grepl('\\.eps', filename))
			postscript(filename, width = file.width, height = file.height, horizontal = F, onefile = FALSE, paper = "special")
		else if (grepl('\\.png', filename))
			png(filename, width = file.width, height = file.height, units = 'in', res = 600)
		else if (grepl('\\.tif', filename))
			tiff(filename, width = file.width, height = file.height, units = 'in', res = 600)
	}
	# Get vectors of start and end points for each chromosome (ie: cumulative sizes + gaps)
	ce <- cumsum(chrom.sizes + c(0, 0, gaps, 0, gaps))
	cs <- ce - chrom.sizes
	layout(matrix(c(rep(1,4),rep(2,1)), nrow = 5, ncol = 1))
	par(mar = c(0,4,1,2), mgp = c(2, 0.7, 0), family = 'Arial', xpd = NA) 
	num.randomisations <- ifelse(is.null(num.randomisations), 
	                             length(randomisation.columns),
	                             num.randomisations)
	# This odd way of getting a sequence is to make sure we get the right outcome if num.randomisations == 0
	r.columns <- randomisation.columns[seq(1, num.randomisations, length.out = num.randomisations)]
	fst.columns <- c(phenotype.column, r.columns)
	max.y <- max(c(max(fst.table[, ..fst.columns]), 0.05))
	min.y <- min(fst.table[, ..fst.columns])
	# Create the empty plot.
	plot(c(cs[1], ce[5]), c(min.y, max.y), xlim = c(cs[1] + gaps/2, ce[5] - gaps/2), type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', main = plot.title, cex.main = 1.2)
	
	# Get plotting positions for windows
	fst.table$genome.pos <- fst.table$window.pos + cs[as.character(fst.table$window.chrom)]
	# Add true data
	fst.table[, lines(genome.pos, phenotype.moving.Fst, col = 'orangered3', lwd = 1.2), by = window.chrom]
	# Add y axis
	fst.step.size <- ifelse(max.y > 0.2, 0.1, 0.05)
	axis(2, at = seq(0, max.y, fst.step.size))
	mtext('Fst', 2, 2, cex = 0.8)
	# Add peaks that pass filtering. Colour according to whether they have significant p-value
	if (length(p.thresh) == 1)
		p.thresh <- c(p.thresh, -1)
	p.colours <- c('green', 'orchid3', 'blue')
	filter.pass <- fst.table[[filter.name]]
	fst.table[filter.pass, points(genome.pos, phenotype.moving.Fst, 
	                              pch = 21,
	                              bg = p.colours[(window.p < p.thresh[1]) + (window.p < p.thresh[2]) + 1], 
	                              col = 'grey40', cex = 1.1, lwd = .5)
	         ]
	# For the abline, we always use the population's own filter. 
	abline(h = min(fst.table[(filter), phenotype.moving.Fst]), col = 'blue', xpd = F)
	
	# Now plot all chromosomes with, the position of each of the four detox gene regions and Ace1
	par(mar = c(1,4,0,2), mgp = c(2, 0.7, 0)) 
	add.chromosomes(gaps = gaps, cs = cs, ce = ce, gene.cex = 0.7, point.cex = 1, chrom.offset = -1.2, chrom.cex = 1.2)
	
	if (!missing(filename))
		dev.off()
}

# Having explored The possibility of using joint filters, it doesn't do all that much good. Just a couple of 
# "sensible" new windows picked up, but only with p < 0.01 (rather than p < 0.001), in which case we also get
# some "non-sensible" new windows alongside. So we stick to within-population filters. 
for (pop in names(fst.tables)){
	plot.fst.filters(fst.tables[[pop]], filename = paste(pop, 'peak_filter_plot.png', sep = '_'), plot.title = pop)
	#plot.fst.filters(fst.tables[[pop]], plot.title = pop, filter.name = 'joined.filter')
}

# Next we want to use the randomised data to find which of those peaks we believe. 
# For each population, pull out the peaks. 
filtered.tables <- lapply(fst.tables, function(D) D[filter == T & window.p < 0.001, ])

save.image('fst_filtered_windows.Rdata')



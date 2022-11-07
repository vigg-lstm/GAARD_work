library(stringi)
library(data.table)
library(magrittr)
library(future.apply)
plan(tweak(multisession, workers = 20))

# Load all of the H12 tables
h12.filenames <- list.files('H12_outputs', '*.tsv', full.names = T)

study.pops <- c('Avrankou.coluzzii.Delta',
                'Baguida.gambiae.Delta',
                'Baguida.gambiae.PM',
                'Korle-Bu.coluzzii.Delta',
                'Korle-Bu.coluzzii.PM',
                'Madina.gambiae.Delta',
                'Madina.gambiae.PM',
                'Obuasi.gambiae.Delta',
                'Obuasi.gambiae.PM')

# Get the randomisation ids. 
randomisation.ids <- unique(stri_extract_first_regex(h12.filenames, '\\d{4}(?=\\.)')) %>%
                     .[!is.na(.)]
names(randomisation.ids) <- paste('r', randomisation.ids, sep = '')

# A function that looks for positive peaks by identifying windows more extreme than twice the 
# 99% centile:
find.peaks <- function(values, centile = 0.99, multiplier = 2){
	thresh <- quantile(values, centile) * multiplier
	values > thresh
}

# A function to load and combine all data for a given randomisation id
load.and.combine.data <- function(rand.id, pop, find.peaks = F, print.id = F){
	if (print.id)
		cat('\t', rand.id, '\n', sep = '')
	
	alive.filenames <- grep(paste(pop, '[_.]alive', rand.id, '\\.', sep = ''), h12.filenames, value = T) %>%
	                   setNames(., stri_extract_first_regex(., '[23LRX]+(?=\\.tsv)'))
	dead.filenames <- grep(paste(pop, '[_.]dead', rand.id, '\\.', sep = ''), h12.filenames, value = T) %>%
	                  setNames(., stri_extract_first_regex(., '[23LRX]+(?=\\.tsv)'))
	
	if (any(names(dead.filenames) != names(alive.filenames)))
		stop(paste('Different chromosomes have been found for dead and alive samplesets for randomisation', rand.id))
	
	alive.data <- names(alive.filenames) %>%
	              lapply(function(chrom){ 
	                  fread(alive.filenames[chrom], sep = '\t') %>%
	                  .[, .(chromosome = ..chrom, midpoint = midpoint, H12 = h12)]
				  }) %>%
				  rbindlist()
	dead.data <- names(dead.filenames) %>%
	              lapply(function(chrom){ 
	                  fread(dead.filenames[chrom], sep = '\t') %>%
	                  .[, .(chromosome = ..chrom, midpoint = midpoint, H12 = h12)]
				  }) %>%
				  rbindlist()
	# Check that the midpoints are identical
	if (!identical(alive.data$midpoint, dead.data$midpoint))
		stop('Window positions do no match between alive and dead samples')
	diff.data <- alive.data[, .(chromosome = chromosome, 
	                            midpoint = midpoint, 
	                            H12.difference = H12 - ..dead.data[, H12])
	]
	
	output <- list(alive = alive.data, dead = dead.data, diff = diff.data)
	
	# Look for peaks in the diff data
	if (find.peaks)
		output$diff[, is.peak := find.peaks(H12.difference)]
	
	output
}

h12.tables <- list()
for (pop in study.pops){
	cat('Subtracting dead from alive H12 data for', pop, '\n')
	cat('\tTrue data\n')
	h12.tables[[pop]] <- load.and.combine.data('', pop = pop, find.peaks = T)
	cat('\tRandomised data\n')
	# Load the randomised data sets
	for (r in names(randomisation.ids)){
		this.randomisation <- load.and.combine.data(randomisation.ids[r], pop = pop, print.id = T)
		# The midpoints need to be the same as in the true data. 
		if (!identical(this.randomisation$diff$midpoint, h12.tables[[pop]]$diff$midpoint))
			stop('Window positions different between true data and randomisations')
		h12.tables[[pop]]$alive[, c(r) := ..this.randomisation$alive[, H12]]
		h12.tables[[pop]]$dead[, c(r) := ..this.randomisation$dead[, H12]]
		h12.tables[[pop]]$diff[, c(r) := ..this.randomisation$diff[, H12.difference]]
	}
	cat('\tCalculating p-value of peaks based on randomised data.\n')
	h12.tables[[pop]]$diff[is.peak == T, pval := apply(.SD - H12.difference, 
	                                                   1, 
	                                                   function(x) sum(x > 0)/length(x)
	                                             ), 
	                                     .SDcols = names(randomisation.ids)
	] 
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

# Function to plot the H12 data
plot.h12.diff <- function(h12.table, 
                          filename = NULL, 
                          num.randomisations = NULL, 
                          p.thresh = 0.01, 
                          plot.title = '', 
                          gaps = 5000000, 
                          filter.name = 'is.peak'){
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
	colours <- c(h12 = add.transparency('orangered3', 0.2),
                 randomisations = add.transparency('grey50', 0.8))
	par(mar = c(0,4,1,2), mgp = c(2, 0.7, 0), family = 'Arial', xpd = NA) 
	num.randomisations <- ifelse(is.null(num.randomisations), 
	                             length(randomisation.ids),
	                             num.randomisations)
	# This odd way of getting a sequence is to make sure we get the right outcome if num.randomisations == 0
	r.columns <- names(randomisation.ids)[seq(1, num.randomisations, length.out = num.randomisations)]
	h12.columns <- c('H12.difference', r.columns)
	max.y <- max(c(max(h12.table[, ..h12.columns]), 0.05))
	min.y <- min(h12.table[, ..h12.columns])
	# Create the empty plot.
	plot(c(cs[1], ce[5]), c(min.y, max.y), xlim = c(cs[1] + gaps/2, ce[5] - gaps/2), type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', main = plot.title, cex.main = 1.2)
	
	# Get plotting positions for windows
	h12.table$genome.pos <- h12.table$midpoint + cs[as.character(h12.table$chromosome)]
	# Add randomised data
	sapply(r.columns, function(x) h12.table[, lines(genome.pos, get(x), col = colours['randomisations'], lwd = 0.4), by = chromosome])
	# Add true data
	h12.table[, lines(genome.pos, H12.difference, col = colours['h12'], lwd = 1.2), by = chromosome]
	# Add y axis
	h12.step.size <- ifelse(max.y > 0.2, 0.1, 0.05)
	axis(2, at = seq(0, max.y, h12.step.size))
	mtext('H12 difference', 2, 2, cex = 0.8)
	# Add peaks that pass filtering. Colour according to whether they have significant p-value
	if (length(p.thresh) == 1)
		p.thresh <- c(p.thresh, -1)
	p.colours <- c('orchid3', 'green', 'blue')
	filter.pass <- h12.table[[filter.name]]
	h12.table[filter.pass, points(genome.pos, H12.difference, 
	                              pch = 21,
	                              bg = p.colours[(pval < p.thresh[1]) + (pval < p.thresh[2]) + 1], 
	                              col = colours['randomisations'], 
	                              cex = 1.1, lwd = .5)
	]
	
	# Now plot all chromosomes with, the position of each of the four detox gene regions and Ace1
	par(mar = c(1,4,0,2), mgp = c(2, 0.7, 0)) 
	add.chromosomes(gaps = gaps, cs = cs, ce = ce, gene.cex = 0.7, point.cex = 1, chrom.offset = -1.2, chrom.cex = 1.2)
	
	if (!missing(filename))
		dev.off()
}

p.threshold = 0.01

for (pop in names(h12.tables))
	plot.h12.diff(h12.tables[[pop]]$diff, 
	              filename = paste(pop, 'peak_filter_plot.png', sep = '_'),
	              p.thresh = p.threshold,
	              plot.title = pop)

save.image('h12_filtered_windows.Rdata')



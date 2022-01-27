library(data.table)
library(stringr)

# We don't include Aboisso since we don't have PBS or Fst data for it
populations <- c('Avrankou_coluzzii_Delta', 
                 'Baguida_gambiae_Delta', 'Baguida_gambiae_PM', 
                 'Korle-Bu_coluzzii_Delta', 'Korle-Bu_coluzzii_PM',
                 'Madina_gambiae_Delta', 'Madina_gambiae_PM',
                 'Obuasi_gambiae_Delta', 'Obuasi_gambiae_PM')

chroms <- c('2L', '2R', '3L', '3R', 'X')

# Load all of the PBS results
pbs.files <- list.files('../PBS/tsv', '*.tsv', full.names = T)
names(pbs.files) <- sub('\\.tsv', '', sub('.*_', '', pbs.files))
pbs.list.by.chrom <- lapply(pbs.files, fread)

join.pbs.chroms <- function(population){
	population.name.split <- strsplit(population, '_')[[1]]
	table.names <- paste(population.name.split[2], population.name.split[1], population.name.split[3], chroms, sep = '.')
	tables <- pbs.list.by.chrom[table.names]
	names(tables) <- str_match(table.names, '(?<=\\.)[23]?[LRX]$')
	do.call(rbind, c(tables, idcol = 'chrom'))
}

pbs.list <- lapply(setNames(nm = populations), join.pbs.chroms)

# Load all of the G12 results
g12.files <- list.files('../G12/tsv', '*.tsv', full.names = T)
names(g12.files) <- sub('\\..*', '', sub('.*\\/', '', g12.files))
g12.list.by.chrom <- lapply(g12.files, fread)

join.h12.chroms <- function(population){
	population.name.split <- strsplit(population, '_')[[1]]
	table.names.alive <- paste(population.name.split[1], population.name.split[2], chroms, population.name.split[3], 'alive', sep = '_')
	table.names.dead <- paste(population.name.split[1], population.name.split[2], chroms, population.name.split[3], 'dead', sep = '_')
	tables.alive <- g12.list.by.chrom[table.names.alive]
	tables.dead <- g12.list.by.chrom[table.names.dead]
	names(tables.alive) <- str_match(table.names.alive, '(?<=_)[23]?[LRX](?=_)')
	names(tables.dead) <- str_match(table.names.dead, '(?<=_)[23]?[LRX](?=_)')
	list(alive = do.call(rbind, c(tables.alive, idcol = 'chrom')),
	     dead = do.call(rbind, c(tables.dead, idcol = 'chrom')))
}

g12.list <- lapply(setNames(nm = populations), join.h12.chroms)

# Load all of the Fst results
fst.files <- list.files('../Fst/tsv', '*.tsv', full.names = T)
names(fst.files) <- sub('\\..*', '', sub('.*\\/', '', fst.files))
fst.list <- lapply(fst.files, fread)

chrom.size <- c('2R' = 61545105, '2L' = 49364325, '3R' = 53200684, '3L' = 41963435, 'X' = 24393108)

# Write a function to add transparency to a colour
add.transparency <- function(col, prop.alpha){
	rgb.val <- col2rgb(col)
	new.rgb <- rgb(rgb.val[1], rgb.val[2], rgb.val[3], max = 255, alpha = (1 - prop.alpha) * 255)
	new.rgb
}

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


# Function to plot PBS, H12 (or G12) and Fst along genome.
plot.stats.along.genome <- function(pbs.table = NULL, h12.table = NULL, fst.table = NULL, filename = NULL, plot.title = '', gaps = 5000000, chrom.sizes = chrom.size){
	if (is.null(pbs.table) & is.null(h12.table) & is.null(fst.table))
		stop('At least one of pbs.table, h12.table or fst.table must be provided.')
	# get a vector of end points for each chromosome (ie: cumulative sizes + gaps)
	ce <- cumsum(chrom.sizes + c(0, 0, gaps, 0, gaps))
	# get a vector of start points for each chromosome (ie: cumulative sizes + gaps)
	cs <- ce - chrom.sizes
	# Create the plot
	if (!is.null(filename)){
		file.width = 6.5
		file.height = 4.5
		if (grepl('\\.eps', filename))
			postscript(filename, width = file.width, height = file.height, horizontal = F, onefile = FALSE, paper = "special")
		else if (grepl('\\.png', filename))
			png(filename, width = file.width, height = file.height, units = 'in', res = 600)
		else if (grepl('\\.tif', filename))
			tiff(filename, width = file.width, height = file.height, units = 'in', res = 600)
	}
	layout(matrix(c(rep(1,4), rep(2,4), rep(3,2)), nrow = 10, ncol = 1))
	par(mar = c(0,4,1,4), mgp = c(2, 0.7, 0), family = 'Arial', xpd = NA) 
	colours <- c(pbs = 'dodgerblue3', 
	             fst = 'darkred', 
	             h12.alive = 'purple', 
	             h12.dead = 'green2')
	# Create the empty plot
	max.y <- max(pbs.table$PBS, fst.table$moving.Fst)
	min.y <- min(pbs.table$PBS, fst.table$moving.Fst)
	plot(c(cs[1], ce[5]), c(min.y, max.y), xlim = c(cs[1] + gaps/2, ce[5] - gaps/2), type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', main = plot.title, cex.main = 1.5)
	
	# For each of the input tables, get a column that shows the position along the whole genome ,rather than just 
	# on the chromosome
	if (!is.null(fst.table)){
		# Get plotting positions for windows
		fst.table$genome.pos <- fst.table$window.pos + cs[as.character(fst.table$window.chrom)]
		max.fst <- max(c(max(fst.table$moving.Fst, na.rm = T), 0.05))
		fst.table$normed.fst <- max.y * fst.table$moving.Fst / max.fst
		# Add data
		fst.table[, lines(genome.pos, normed.fst, col = colours['fst'], lwd = 1.5), by = window.chrom]
		# Add y axis
		fst.step.size <- ifelse(max.fst > 0.2, 0.1, 0.05)
		axis(4, at = seq(0, max.y, max.y * fst.step.size/max.fst), labels = as.character(seq(0, max.fst, fst.step.size)), col = colours['fst'], col.axis = colours['fst'], col.ticks = colours['fst'])
		mtext('Fst', 4, 2, cex = 0.8, col = colours['fst'])
	}
	
	if (!is.null(pbs.table)){
		# Get plotting positions for windows
		pbs.table$genome.pos <- pbs.table$midpoint + cs[as.character(pbs.table$chrom)]
		max.pbs <- max(c(max(pbs.table$PBS, na.rm = T), 0.05))
		pbs.table$normed.pbs <- max.y * pbs.table$PBS / max.pbs
		# Add data
		pbs.table[, lines(genome.pos, normed.pbs, col = add.transparency(colours['pbs'], 0.3), lwd = 1.5), by = chrom]
		# Add y axis
		pbs.step.size <- ifelse(max.pbs > 0.2, 0.1, 0.05)
		axis(2, at = seq(0, max.y, max.y * pbs.step.size/max.pbs), labels = as.character(seq(0, max.pbs, pbs.step.size)), col = colours['pbs'], col.axis = colours['pbs'], col.ticks = colours['pbs'])
		mtext('PBS', 2, 2, cex = 0.8, col = colours['pbs'])
	}
	
	if (!is.null(h12.table)){
		par(mar = c(0,4,1,4), mgp = c(2, 0.7, 0), family = 'Arial', xpd = NA) 
		# Check whether we have G12 or H12
		stat12 <- colnames(h12.table$alive)[grepl('12', colnames(h12.table$alive))]
		max.h12 <- max(c(h12.table$alive[[stat12]], h12.table$dead[[stat12]], 0.05))
		plot(c(cs[1], ce[5]), c(0, max.h12), xlim = c(cs[1] + gaps/2, ce[5] - gaps/2), type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', cex.main = 1.5)
		# Get plotting positions for windows
		h12.table$alive$genome.pos <- h12.table$alive$midpoint + cs[as.character(h12.table$alive$chrom)]
		h12.table$dead$genome.pos <- h12.table$dead$midpoint + cs[as.character(h12.table$dead$chrom)]
		# Add data to plot
		h12.table$alive[, lines(genome.pos, get(stat12), col = colours['h12.alive'], lwd = 2.5), by = chrom]
		h12.table$dead[, lines(genome.pos, get(stat12), col = add.transparency(colours['h12.dead'], 0.4), lwd = 2.5), by = chrom]
		# Add y axis
		h12.step.size <- ifelse(max.h12 > 0.6, 0.2, 0.1)
		axis(2, at = seq(0, max.h12, h12.step.size))
		mtext(stat12, 2, 2, cex = 0.8)
		legend(ce[5], max.h12, c('alive', 'dead'), col = colours[c('h12.alive', 'h12.dead')], lwd = 3, xjust = 0, yjust = 1, bty = 'n')
	}
	
	# Now plot all chromosomes with, the position of each of the four detox gene regions and Ace1
	par(mar = c(1,4,1,4), mgp = c(2, 0.7, 0)) 
	add.chromosomes(gaps = gaps, cs = cs, ce = ce, chrom.offset = -1, chrom.cex = 1.3)
	if (!missing(filename))
		dev.off()
}

# Draw all the plots
sapply(names(fst.list), function(pop) plot.stats.along.genome(pbs.table = pbs.list[[pop]], h12.table = g12.list[[pop]], fst.table = fst.list[[pop]], filename = paste(pop, '_plot.png', sep = ''), plot.title = pop))


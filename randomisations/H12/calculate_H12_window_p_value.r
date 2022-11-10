library(stringi)
library(data.table)
library(magrittr)
library(future.apply)
plan(tweak(multisession, workers = 20))

# Load all of the H12 tables
h12.filenames <- list.files('H12_outputs', '*.tsv', full.names = T)

study.pops <- setNames(nm = c('Avrankou.coluzzii.Delta',
                              'Baguida.gambiae.Delta',
                              'Baguida.gambiae.PM',
                              'Korle-Bu.coluzzii.Delta',
                              'Korle-Bu.coluzzii.PM',
                              'Madina.gambiae.Delta',
                              'Madina.gambiae.PM',
                              'Obuasi.gambiae.Delta',
                              'Obuasi.gambiae.PM'))

# Get the randomisation ids. 
randomisation.ids <- readLines(h12.filenames[1], n = 1) %>%
                     {strsplit(., '\t')[[1]]} %>%
                     grep('r\\d{4}', ., value = T)

# A function that looks for positive peaks by identifying windows more extreme than thrice the 
# 95% centile:
find.peaks <- function(values, centile = 0.95, multiplier = 3){
	thresh <- quantile(values, centile) * multiplier
	values > thresh
}

# A function to load and combine all data for a given randomisation id
load.and.combine.data <- function(pop, find.peaks = T, calculate.P.values = T){
	
	alive.filenames <- grep(paste(pop, '\\.alive', sep = ''), h12.filenames, value = T) %>%
	                   setNames(., stri_extract_first_regex(., '[23LRX]+(?=\\.tsv)'))
	dead.filenames <- grep(paste(pop, '\\.dead', '\\.', sep = ''), h12.filenames, value = T) %>%
	                  setNames(., stri_extract_first_regex(., '[23LRX]+(?=\\.tsv)'))
	
	if (any(names(dead.filenames) != names(alive.filenames)))
		stop(paste('Different chromosomes have been found for dead and alive samplesets for randomisation', rand.id))
	
	alive.data <- names(alive.filenames) %>%
	              lapply(function(chrom){ 
	                  fread(alive.filenames[chrom], sep = '\t') %>%
	                  .[, chromosome := ..chrom] %>%
					  setnames('h12', 'H12') %>%
					  setcolorder(c('chromosome', 'startpoint', 'endpoint', 'midpoint'))
				  }) %>%
				  rbindlist()
	dead.data <- names(dead.filenames) %>%
	             lapply(function(chrom){ 
	                 fread(dead.filenames[chrom], sep = '\t') %>%
	                 .[, chromosome := ..chrom] %>%
					 setnames('h12', 'H12') %>%
					 setcolorder(c('chromosome', 'startpoint', 'endpoint', 'midpoint'))
				 }) %>%
				 rbindlist()
	# Check that the midpoints are identical
	if (!identical(alive.data$midpoint, dead.data$midpoint))
		stop('Window positions do no match between alive and dead samples')
	diff.data <- cbind(alive.data[, .(chromosome, startpoint, endpoint, midpoint)],
	                   alive.data[, c('H12', ..randomisation.ids)] - dead.data[, c('H12', ..randomisation.ids)]
	)
	
	# Look for peaks in the diff data
	if (find.peaks){
		diff.data[, is.peak := find.peaks(H12)]
		setcolorder(diff.data, c('chromosome', 'startpoint', 'endpoint', 'midpoint', 'H12', 'is.peak'))
	}
	
	if (calculate.P.values){
		if(find.peaks){
			diff.data[is.peak == T, pval := apply(.SD - H12, 
			                                      1, 
			                                      function(x) sum(x > 0)/length(x)
			                                ), 
			                                .SDcols = randomisation.ids
			] 
		}
		else {
			diff.data[, pval := apply(.SD - H12, 
			                          1, 
			                          function(x) sum(x > 0)/length(x)
			                    ), 
			                    .SDcols = randomisation.ids
			] 
		}
	}
	
	output <- list(alive = alive.data, dead = dead.data, diff = diff.data)
	output
}

cat('Loading data\n')
h12.tables <- lapply(study.pops, function(pop) {cat(pop, '\n'); load.and.combine.data(pop)})

source('../../shared_functions/R_plotting.r')

chrom.sizes <- c('2R' = 61545105, '2L' = 49364325, '3R' = 53200684, '3L' = 41963435, 'X' = 24393108)

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
	colours <- c(h12 = add.transparency('limegreen', 0.2),
                 randomisations = add.transparency('grey50', 0.85))
	par(mar = c(0,4,1,2), mgp = c(2, 0.7, 0), family = 'Arial', xpd = NA) 
	num.randomisations <- ifelse(is.null(num.randomisations), 
	                             length(randomisation.ids),
	                             num.randomisations)
	# This odd way of getting a sequence is to make sure we get the right outcome if num.randomisations == 0
	r.columns <- randomisation.ids[seq(1, num.randomisations, length.out = num.randomisations)]
	h12.columns <- c('H12', r.columns)
	max.y <- max(c(max(h12.table[, ..h12.columns]), 0.05))
	min.y <- min(h12.table[, ..h12.columns])
	# Create the empty plot.
	plot(c(cs[1], ce[5]), c(min.y, max.y), xlim = c(cs[1] + gaps/2, ce[5] - gaps/2), type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', main = plot.title, cex.main = 1.2)
	
	# Get plotting positions for windows
	h12.table$genome.pos <- h12.table$midpoint + cs[as.character(h12.table$chromosome)]
	# Add randomised data
	sapply(r.columns, function(x) h12.table[, lines(genome.pos, get(x), col = colours['randomisations'], lwd = 0.4), by = chromosome])
	# Add true data
	h12.table[, lines(genome.pos, H12, col = colours['h12'], lwd = 1.2), by = chromosome]
	# Add y axis
	h12.step.size <- ifelse(max.y > 0.2, 0.1, 0.05)
	axis(2, at = seq(0, max.y, h12.step.size))
	mtext('H12 difference', 2, 2, cex = 0.8)
	# Add peaks that pass filtering. Colour according to whether they have significant p-value
	if (length(p.thresh) == 1)
		p.thresh <- c(p.thresh, -1)
	p.colours <- c('orchid3', 'green', 'blue')
	filter.pass <- h12.table[[filter.name]]
	h12.table[filter.pass, points(genome.pos, H12, 
	                              pch = 21,
	                              bg = p.colours[(pval < p.thresh[1]) + (pval < p.thresh[2]) + 1], 
	                              col = colours['randomisations'], 
	                              cex = 1.1, lwd = .5)
	]
	
	# Now plot all chromosomes with, the position of each of the four detox gene regions and Ace1
	par(mar = c(1,4,0,2), mgp = c(2, 0.7, 0)) 
	add.chromosomes(chrom.sizes, gaps = gaps, gene.cex = 0.7, point.cex = 1, chrom.offset = -1.2, chrom.cex = 1.2)
	
	if (!missing(filename))
		dev.off()
}

p.threshold = 0.01

for (pop in names(h12.tables))
	plot.h12.diff(h12.tables[[pop]]$diff, 
	              filename = paste(pop, 'peak_filter_plot.png', sep = '_'),
	              p.thresh = p.threshold,
	              plot.title = pop)

saveRDS(h12.tables, 'h12_filtered_windows.RDS')



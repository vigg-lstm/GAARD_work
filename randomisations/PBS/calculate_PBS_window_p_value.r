library(stringi)
library(data.table)
library(magrittr)
library(future.apply)
plan(tweak(multisession, workers = 20))

# Load all of the PBS tables
pbs.filenames <- list.files('PBS_outputs', '*.tsv', full.names = T)

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
randomisation.ids <- readLines(pbs.filenames[1], n = 1) %>%
                     {strsplit(., '\t')[[1]]} %>%
                     grep('r\\d{4}', ., value = T)

# A function that looks for positive peaks by identifying windows more extreme than thrice the 
# 95% centile of the non-2La data:
find.peaks <- function(pbs.table, centile = 0.95, multiplier = 3){
	thresh <- quantile(pbs.table[chromosome != '2L']$PBS, centile) * multiplier
	pbs.table$PBS > thresh
}

# A function to load and combine all data for a given randomisation id
load.and.combine.data <- function(pop, peak.function = NULL, calculate.P.values = T){
	
	filenames <- grep(pop, pbs.filenames, value = T) %>%
	             setNames(., stri_extract_first_regex(., '[23LRX]+(?=\\.tsv)'))
	
	pbs.data <- names(filenames) %>%
	            lapply(function(chrom){ 
	                fread(filenames[chrom], sep = '\t') %>%
	                .[, chromosome := ..chrom] %>%
					setnames('pbs', 'PBS') %>%
					setcolorder(c('chromosome', 'startpoint', 'endpoint', 'midpoint'))
				}) %>%
				rbindlist()
	
	# Look for peaks in the PBS data
	if (!is.null(peak.function))
		pbs.data[, is.peak := peak.function(.SD)]
	
	if (calculate.P.values){
		if (!is.null(peak.function)) {
			pbs.data[is.peak == T, pval := apply(.SD - PBS, 
			                                     1, 
			                                     function(x) sum(x > 0)/length(x)
			                               ), 
			                               .SDcols = randomisation.ids
			] 
		}
		else {
			pbs.data[, pval := apply(.SD - PBS, 
			                         1, 
			                         function(x) sum(x > 0)/length(x)
			                   ), 
			                   .SDcols = randomisation.ids
			] 
		}
	}
	
	pbs.data
}

cat('Loading data\n')
pbs.tables <- lapply(study.pops, function(pop) {cat(pop, '\n'); load.and.combine.data(pop, find.peaks)})

source('../../shared_functions/R_plotting.r')

chrom.sizes <- c('2R' = 61545105, '2L' = 49364325, '3R' = 53200684, '3L' = 41963435, 'X' = 24393108)

# Function to plot the PBS data
plot.pbs <- function(pbs.table, 
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
	colours <- c(pbs = add.transparency('royalblue1', 0.2),
                 randomisations = add.transparency('grey50', 0.8))
	par(mar = c(0,4,1,2), mgp = c(2, 0.7, 0), family = 'Arial', xpd = NA) 
	num.randomisations <- ifelse(is.null(num.randomisations), 
	                             length(randomisation.ids),
	                             num.randomisations)
	# This odd way of getting a sequence is to make sure we get the right outcome if num.randomisations == 0
	r.columns <- randomisation.ids[seq(1, num.randomisations, length.out = num.randomisations)]
	pbs.columns <- c('PBS', r.columns)
	max.y <- max(c(max(pbs.table[, ..pbs.columns]), 0.05))
	min.y <- min(pbs.table[, ..pbs.columns])
	# Create the empty plot.
	plot(c(cs[1], ce[5]), c(min.y, max.y), xlim = c(cs[1] + gaps/2, ce[5] - gaps/2), type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', main = plot.title, cex.main = 1.2)
	
	# Get plotting positions for windows
	pbs.table$genome.pos <- pbs.table$midpoint + cs[as.character(pbs.table$chromosome)]
	# Add randomised data
	sapply(r.columns, function(x) pbs.table[, lines(genome.pos, get(x), col = colours['randomisations'], lwd = 0.4), by = chromosome])
	# Add true data
	pbs.table[, lines(genome.pos, PBS, col = colours['pbs'], lwd = 1.2), by = chromosome]
	# Add y axis
	pbs.step.size <- ifelse(max.y > 0.2, 0.1, 0.05)
	axis(2, at = seq(0, max.y, pbs.step.size))
	mtext('PBS', 2, 2, cex = 0.8)
	# Add peaks that pass filtering. Colour according to whether they have significant p-value
	if (length(p.thresh) == 1)
		p.thresh <- c(p.thresh, -1)
	p.colours <- c('orchid3', 'green', 'blue')
	filter.pass <- pbs.table[[filter.name]]
	pbs.table[filter.pass, points(genome.pos, PBS, 
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

for (pop in names(pbs.tables))
	plot.pbs(pbs.tables[[pop]], 
	         filename = paste(pop, 'peak_filter_plot.png', sep = '_'),
	         p.thresh = p.threshold,
	         plot.title = pop)

saveRDS(pbs.tables, 'pbs_filtered_windows.RDS')



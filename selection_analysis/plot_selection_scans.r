library(stringi)
library(stringr)
library(data.table)
library(magrittr)
library(ape)

# Load all of the H12 tables
h12.filenames <- list.files('.', 'H12_.*.csv', full.names = T)

study.pops <- setNames(nm = c('Avrankou',
                              'Aboisso',
                              'Baguida',
                              'Korle-Bu',
                              'Madina',
                              'Obuasi'))

# A function to load and combine all data for a given randomisation id
load.and.combine.h12 <- function(pop){
	
	filenames <- grep(pop, h12.filenames, value = T) %>%
	             setNames(., stri_extract_first_regex(., '[23LRX]+(?=\\.csv)'))
	
	combined.data <- names(filenames) %>%
	                 lapply(function(chrom){ 
	                     fread(filenames[chrom], sep = '\t') %>%
	                     .[, chromosome := ..chrom] %>%
					     setnames('h12', 'H12') %>%
					     setcolorder(c('chromosome', 'midpoint'))
				     }) %>%
				     rbindlist()
	
	combined.data
}

cat('Loading H12 data\n')
h12.tables <- lapply(study.pops, function(pop) {cat(pop, '\n'); load.and.combine.h12(pop)})

source('../shared_functions/R_plotting.r')

chrom.sizes <- c('2R' = 61545105, '2L' = 49364325, '3R' = 53200684, '3L' = 41963435, 'X' = 24393108)

# Function to plot the H12 data
plot.h12 <- function(h12.table, 
                     filename = NULL, 
                     plot.title = '', 
                     gaps = 5000000, 
                     h12.col = 'darkgreen',
                     y.label = 'H12'){
	# Create the plot
	if (!missing(filename)){
		file.width = 6.5
		file.height = 2
		if (grepl('\\.eps', filename)){
			postscript(filename, width = file.width, height = file.height, horizontal = F, onefile = FALSE, paper = "special")
			par(mar = c(0,4,1,2), mgp = c(2, 0.7, 0), family = 'Arial', xpd = NA) 
		}
		else if (grepl('\\.png', filename)){
			png(filename, width = file.width, height = file.height, units = 'in', res = 600)
			par(mar = c(0,4,1,2), mgp = c(2, 0.7, 0), family = 'Arial', xpd = NA) 
		}
		else if (grepl('\\.tif', filename)){
			tiff(filename, width = file.width, height = file.height, units = 'in', res = 600)
			par(mar = c(0,4,1,2), mgp = c(2, 0.7, 0), family = 'Arial', xpd = NA) 
		}
		else if (grepl('\\.pdf', filename)){
			pdf(filename, width = file.width, height = file.height)
			par(mar = c(0,4,1,2), mgp = c(2, 0.7, 0), xpd = NA) 
		}
	}
	layout(matrix(c(rep(1,7),rep(2,3)), nrow = 10, ncol = 1))
	# Get vectors of start and end points for each chromosome (ie: cumulative sizes + gaps)
	ce <- cumsum(chrom.sizes + c(0, 0, gaps, 0, gaps))
	cs <- ce - chrom.sizes
	max.y <- max(c(max(h12.table[, 'H12']), 0.05))
	min.y <- min(h12.table[, 'H12'])
	# Create the empty plot.
	plot(c(cs[1], ce[5]), c(min.y, max.y), xlim = c(cs[1] + gaps/2, ce[5] - gaps/2), 
	     type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', 
	     main = plot.title, cex.main = 1.2)
	
	# Get plotting positions for windows
	h12.table$genome.pos <- h12.table$midpoint + cs[as.character(h12.table$chromosome)]
	# Add the data
	h12.table[, lines(genome.pos, H12, col = h12.col, lwd = 1), by = chromosome]
	# Add y axis
	h12.step.size <- ifelse(max.y > 0.2, 0.2, 0.05)
	axis(2, at = seq(0, max.y, h12.step.size))
	mtext(y.label, 2, 2, cex = 0.8)
	
	# Now plot all chromosomes with, the position of each of the four detox gene regions and Ace1
	par(mar = c(1,4,0,2), mgp = c(2, 0.7, 0)) 
	add.chromosomes(chrom.sizes, gaps = gaps, gene.cex = 0.6, point.cex = 0.6, chrom.offset = 0, chrom.cex = 0.8)
	
	if (!missing(filename))
		dev.off()
}

for (pop in names(h12.tables))
	plot.h12(h12.tables[[pop]], 
	         filename = paste(pop, 'H12.pdf', sep = '_'),
	         plot.title = pop)


# Now all the same for H1x
# Load all of the H12 tables
h1x.filenames <- list.files('.', 'H1x_.*.csv', full.names = T)

pop.pairs <- setNames(nm = c('Avrankou_vs_Aboisso',
                             'Avrankou_vs_Baguida',
                             'Avrankou_vs_Korle-Bu',
                             'Avrankou_vs_Madina',
                             'Avrankou_vs_Obuasi',
                             'Aboisso_vs_Baguida',
                             'Aboisso_vs_Korle-Bu',
                             'Aboisso_vs_Madina',
                             'Aboisso_vs_Obuasi',
                             'Baguida_vs_Korle-Bu',
                             'Baguida_vs_Madina',
                             'Baguida_vs_Obuasi',
                             'Korle-Bu_vs_Madina',
                             'Korle-Bu_vs_Obuasi',
                             'Madina_vs_Obuasi'))

# A function to load and combine all data for a given randomisation id
load.and.combine.h1x <- function(pop1, pop2){
	
	filenames <- grep(pop1, h1x.filenames, value = T) %>%
	             grep(pop2, ., value = T) %>%
	             setNames(., stri_extract_first_regex(., '[23LRX]+(?=\\.csv)'))
	
	combined.data <- names(filenames) %>%
	                 lapply(function(chrom){ 
	                     fread(filenames[chrom], sep = '\t') %>%
	                     .[, chromosome := ..chrom] %>%
					     setnames('h1x', 'H1x') %>%
					     setcolorder(c('chromosome', 'midpoint'))
				     }) %>%
				     rbindlist()
	
	combined.data
}

cat('Loading data\n')
h1x.tables <- lapply(strsplit(pop.pairs, '_vs_'),
                     function(pops) {print(pops); load.and.combine.h1x(pops[1], pops[2])}
)

# Function to plot the H12 data
plot.h1x <- function(h1x.table, 
                     filename = NULL, 
                     plot.title = '', 
                     gaps = 5000000, 
                     h1x.column = 'H1x',
                     h1x.colour = 'darkgreen',
                     y.label = 'H1X'){
	# Create the plot
	if (!missing(filename)){
		file.width = 6.5
		file.height = 2
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
	layout(matrix(c(rep(1,7),rep(2,3)), nrow = 10, ncol = 1))
	par(mar = c(0,4,1,2), mgp = c(2, 0.7, 0), family = 'Arial', xpd = NA) 
	max.y <- max(c(max(h1x.table[, ..h1x.column]), 0.05))
	min.y <- min(h1x.table[, ..h1x.column])
	# Create the empty plot.
	plot(c(cs[1], ce[5]), c(min.y, max.y), xlim = c(cs[1] + gaps/2, ce[5] - gaps/2), 
	     type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', 
	     main = plot.title, cex.main = 1.2)
	
	# Get plotting positions for windows
	h1x.table$genome.pos <- h1x.table$midpoint + cs[as.character(h1x.table$chromosome)]
	# Add the data
	h1x.table[, lines(genome.pos, get(h1x.column), col = h1x.colour, lwd = 1.2), by = chromosome]
	# Add y axis
	h1x.step.size <- ifelse(max.y > 0.2, 0.2, 0.05)
	axis(2, at = seq(0, max.y, h1x.step.size))
	mtext(y.label, 2, 2, cex = 0.8)
	
	# Now plot all chromosomes with, the position of each of the four detox gene regions and Ace1
	par(mar = c(1,4,0,2), mgp = c(2, 0.7, 0)) 
	add.chromosomes(chrom.sizes, gaps = gaps, gene.cex = 0.6, point.cex = 0.6, chrom.offset = 0, chrom.cex = 0.8)
	
	if (!missing(filename))
		dev.off()
}

for (pop.pair in names(h1x.tables))
	plot.h1x(h1x.tables[[pop.pair]], 
	         filename = paste(pop.pair, 'H1x.png', sep = '_'),
	         plot.title = pop.pair,
	         h1x.column = 'H1x',
	         h1x.colour = 'purple',
	         y.label = 'H1x')


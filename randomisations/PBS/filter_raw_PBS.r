library(stringi)
library(data.table)
library(magrittr)
library(ggplot2)

# Load all of the H12 tables
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

# A function to load and combine all data for a population
load.and.combine.data <- function(pop, peak.function = NULL, print.id = F){
	if (print.id)
		cat('\t', rand.id, '\n', sep = '')
	
	filenames <- grep(pop, pbs.filenames, value = T) %>%
	             setNames(., stri_extract_first_regex(., '[23LRX]+(?=\\.tsv)'))
	
	pbs.data <- names(filenames) %>%
	            lapply(function(chrom) {
	                       fread(filenames[chrom], sep = '\t') %>%
	            	       .[, .(chromosome = ..chrom, midpoint = midpoint, PBS = PBS)]
	                   }
	            ) %>%
	            rbindlist()
	
	# Look for peaks in the diff data
	if (!is.null(peak.function))
		pbs.data[, is.peak := peak.function(.SD)]
	
	pbs.data
}

pbs.tables <- lapply(study.pops, load.and.combine.data)

# Where should we set the threshold for identifying peaks. Populations differ a lot in the shape
# of their distributions at the extremes (strongly driven by whether or not there is a signal in
# 2La) but less so closer to the middle. We can look at how the distribution similarity changes 
# as we get further from the median:
quantiles.table <- sapply(study.pops, function(s) quantile(pbs.tables[[s]][chromosome != '2L']$PBS, seq(0,1,0.01))) %>%
                   data.table() %>%
                   melt(measure.vars = study.pops, variable.name = 'Study', value.name = 'PBS') %>%
				   .[, quant := rep(round(seq(0,1,0.01), 2), length(..study.pops))]
quantile.comparison <- ggplot(data = quantiles.table[quant > 0.5], aes(quant, PBS, colour = Study)) + 
                       geom_line() +
                       scale_y_continuous(trans='log2')
x11()
print(quantile.comparison)

# As long as we don't include 2L, we can use the 95% centile * 3, which seems to work quite well:
thresholds <- quantiles.table[quant == 0.95][, thresh := PBS*3]
for (s in study.pops){
	x11()
	par(mfrow = c(5,1), mar = c(2,2,1,1), mgp = c(1.2, 0.6, 0))
	pbs.tables[[s]][, {
		plot(midpoint, PBS, pch = 19, cex = 0.3, col = 'purple', main = ..s)
		abline(h = ..thresholds[Study == ..s, thresh], col = 'red')
	}, by = 'chromosome']
}
# That looks good

# So now reload all the data, looking for peaks

# A function that identifies the threshold to use as a cutoff to find peaks. 
find.peaks <- function(pbs.table, centile = 0.95, multiplier = 3){
	thresh <- quantile(pbs.table[chromosome != '2L']$PBS, centile) * multiplier
	pbs.table$PBS > thresh
}

pbs.tables <- lapply(study.pops, load.and.combine.data, find.peaks)

# Subset the tables to include only peak windows
pbs.tables.filtered <- lapply(study.pops, function(pop) pbs.tables[[pop]][is.peak == T])

# Save these to file. 
save.pbs.table <- function(pop){
	pbs.tables.filtered[[pop]][, fwrite(.SD, 
	                                    paste('PBS_outputs/PBS_', ..pop,
	                                          '_filtered_', .BY, '.tsv', sep = ''
	                                    ), 
	                                    sep = '\t'
	                             ),
	                             .SDcols = c('midpoint', 'PBS'),
	                             by = 'chromosome'
	]
}

lapply(study.pops, save.pbs.table)





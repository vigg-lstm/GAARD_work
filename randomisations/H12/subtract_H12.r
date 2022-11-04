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

# We will only work on randomisation ids up to 20 for now, for the practice run
randomisation.ids <- randomisation.ids[1:10]

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
# Somewhat convoluted way of getting the min and max y value.
all.h12s <- c('H12', names(randomisation.ids))
h12.limits <- list()
for (pop in study.pops){
	h12.limits[[pop]] <- h12.tables[[pop]] %>%
	                     {c(.$dead[, ..all.h12s], .$alive[, ..all.h12s])} %>%
				         unlist %>%
				         {c(min(.), max(.))}
}

for (pop in study.pops){
pop = 'Madina.gambiae.PM'
chrom <- '2R'
x11()
par(mfrow=c(2,1), mar = c(3, 3, 1.5, 0), mgp = c(2,0.8,0))
plot(c(1, chrom.sizes[chrom]), h12.limits[[pop]], type = 'n', xlab = chrom, ylab = 'H12', main = pop)
# plot the randomised data in grey
h12.tables[[pop]]$dead[chromosome == chrom, 
                       lapply(names(.SD), function(R) lines(midpoint, 
                                                           get(R),
                                                           col = add.transparency('grey50', 0.7),
                                                           lwd = 1))
                   ]
# add the true data
h12.tables[[pop]]$dead[chromosome == chrom, lines(midpoint, H12, col = add.transparency('red', 0.7), lwd = 2)]

plot(c(1, chrom.sizes[chrom]), h12.limits[[pop]], type = 'n', xlab = chrom, ylab = 'H12')
# plot the randomised data in grey
h12.tables[[pop]]$alive[chromosome == chrom, 
                       lapply(names(.SD), function(R) lines(midpoint, 
                                                           get(R),
                                                           col = add.transparency('grey50', 0.7),
                                                           lwd = 1))
                   ]
# add the true data
h12.tables[[pop]]$alive[chromosome == chrom, lines(midpoint, H12, col = add.transparency('blue', 0.7), lwd = 2)]
}

h12.diff.limits <- list()
all.h12.diffs <- c('H12.difference', names(randomisation.ids))
for (pop in study.pops){
	h12.diff.limits[[pop]] <- h12.tables[[pop]] %>%
	                          {.$diff[, ..all.h12.diffs]} %>%
				              unlist %>%
				              {c(min(.), max(.))}
}

for (pop in study.pops){
	pop = 'Madina.gambiae.PM'
	chrom <- '2R'
	x11()
	plot(c(1, chrom.sizes[chrom]), h12.diff.limits[[pop]], type = 'n', xlab = chrom, ylab = 'H12 difference', main = pop)
	# plot the randomised data in grey
	h12.tables[[pop]]$diff[chromosome == chrom, 
						   lapply(names(.SD), function(R) lines(midpoint, 
															    get(R),
															    col = add.transparency('grey50', 0.7),
															    lwd = 1))
					      ]
	# add the true data
	h12.tables[[pop]]$diff[chromosome == chrom, lines(midpoint, H12.difference, col = add.transparency('red', 0.7), lwd = 2)]

}

###
#h12.tables[[pop]][[1]]$dead[chromosome == '2R', lines(midpoint, H12, col = add.transparency('red', 0.4), lwd = 3)]
#h12.tables[[pop]][[1]]$alive[chromosome == '2R', lines(midpoint, H12, col = add.transparency('green', 0.4), lwd = 3)]
##
x11()
par(mfrow = c(1,3))
h12.tables[[pop]][[1]]$alive[chromosome == '2R', hist(H12, breaks = 100)]
h12.tables[[pop]][[1]]$dead[chromosome == '2R', hist(H12, breaks = 100)]
h12.tables[[pop]][[1]]$diff[chromosome == '2R', hist(H12.difference, breaks = 100)]
#
x11()
par(mfrow = c(2,1))
h12.tables[[pop]][[1]]$diff[, hist(H12.difference, breaks = 1000)]
h12.tables[[pop]][[1]]$diff[, hist(H12.difference, breaks = 1000, ylim = c(0,10))]



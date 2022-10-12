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
randomisation.ids <- unique(stri_extract_first_regex(h12.filenames, '\\d{4}(?=\\.)'))
names(randomisation.ids) <- paste('r', randomisation.ids, sep = '')

# We will only work on randomisation ids up to 20 for now, for the practice run
randomisation.ids <- randomisation.ids[1:20]

# We have the problem that the genomic windows for H12 are not the same for the dead and the alive. There
# are a few ways in which we could get around this, and I don't think it will make much difference which 
# one we choose. 
# We will do the following. First, we identify which of the dead / alive datasets has the most windows. 
# Then we take that dataset and, for each window, we identify the closest window from the opposite dataset. 
# We then subtract dead from alive at all those windows

# First, we write a function that takes two vectors of table positions and, for each window in the first 
# vector, finds the nearest window in the second.
find.nearest.windows <- function(pos.1, pos.2){
	distances <- abs(outer(pos.1, pos.2, '-'))
	nearest.windows <- apply(distances, 1, which.min)
	nearest.windows
}

# Now a function that applies that to take the difference between nearest windows. 
get.h12.difference <- function(alive.table, dead.table){
	if (nrow(alive.table) > nrow(dead.table)){
		nearest.window.dead.h12 <- dead.table[find.nearest.windows(alive.table$midpoint, dead.table$midpoint), H12]
		diff.table <- alive.table[, .(chromosome = chromosome, midpoint = midpoint, H12.difference = H12 - ..nearest.window.dead.h12)]
	}
	else {
		nearest.window.alive.h12 <- alive.table[find.nearest.windows(dead.table$midpoint, alive.table$midpoint), H12]
		diff.table <- dead.table[, .(chromosome = chromosome, midpoint = midpoint, H12.difference = ..nearest.window.alive.h12 - H12)]
	}
	diff.table
}

# Or we do a linear interpolation. 
find.nearest.windows.either.side <- function(pos.1, pos.2){
	wins <- t(sapply(pos.1, function(x) c(max(which(pos.2 <= x)), min(which(pos.2 > x)))))
	wins
}

get.h12.interpolation <- function(reference.table, interpol.table, h12.column = 'H12'){
	# Reference table is the one that we will use to decide which genetic positions to use.
	# Interpol table is the one for which we will interpolate values between positions in 
	# order to match Reference table. 
	# Get the nearest windows in interpol.table to the positions in reference.table
	nearest.win <- find.nearest.windows.either.side(reference.table$midpoint, interpol.table$midpoint)
	nearest.win.val <- interpol.table[[h12.column]][nearest.win] %>%
	                   matrix(nrow(nearest.win), 2)
	nearest.win.pos <- interpol.table$midpoint[nearest.win] %>%
	                   matrix(nrow(nearest.win), 2)
	# For every position in the reference table, find out where it lies relative to the windows
	# either side of it in the interpol table. 
	interpolation.fraction <- (reference.table$midpoint - nearest.win.pos[,1]) / (nearest.win.pos[,2] - nearest.win.pos[,1])
	# Calculate the interpolated values in the interpol table
	interpolation.value <- nearest.win.val[,1] + interpolation.fraction * (nearest.win.val[,2] - nearest.win.val[,1])
	# Where there is an NA, it's because the window position in reference.table was smaller than the 
	# smallest (or larger than the largest) window in interpol.table, meaning that it was impossible
	# to interpolate. So we just assign the closest window value instead. 
	interpolation.na <- is.na(interpolation.value)
	single.nearest.win <- find.nearest.windows(reference.table[interpolation.na, midpoint], interpol.table[, midpoint])
	interpolation.value[interpolation.na] <- interpol.table[single.nearest.win][[h12.column]]
	interpolation.value
}

get.h12.difference <- function(alive.table, dead.table){
	if (nrow(alive.table) > nrow(dead.table)){
		interpolated.dead.h12 <- get.h12.interpolation(alive.table, dead.table)
		diff.table <- alive.table[, .(chromosome = chromosome, midpoint = midpoint, H12.difference = H12 - ..interpolated.dead.h12)]
	}
	else {
		interpolated.alive.h12 <- get.h12.interpolation(dead.table, alive.table)
		diff.table <- dead.table[, .(chromosome = chromosome, midpoint = midpoint, H12.difference = ..interpolated.alive.h12 - H12)]
	}
	diff.table
}



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
	              lapply(function(chrom) fread(alive.filenames[chrom], sep = '\t')[, .(chromosome = ..chrom, midpoint = midpoint, H12 = H12)])
	dead.data <- names(dead.filenames) %>%
	             lapply(function(chrom) fread(dead.filenames[chrom], sep = '\t')[, .(chromosome = ..chrom, midpoint = midpoint, H12 = H12)])
	diff.data <- 1:length(alive.data) %>%
				 lapply(function(i) get.h12.difference(alive.data[[i]], dead.data[[i]]))
	
	output <- list(alive = rbindlist(alive.data), dead = rbindlist(dead.data), diff = rbindlist(diff.data))
	
	# Look for peaks in the diff data
	if (find.peaks)
		output$diff[, is.peak := find.peaks(H12.difference)]
	
	output
}

h12.tables <- list()
for (pop in study.pops){
	h12.tables[[pop]] <- list()
	cat('Subtracting dead from alive H12 data for', pop, '\n')
	cat('\tTrue data\n')
	h12.tables[[pop]][['true']] <- load.and.combine.data('', pop = pop, find.peaks = T)
	cat('\tRandomised data\n')
	# For each of the randomised data sets, add the interpolated values to the true data set
	for (r in names(randomisation.ids)){
		h12.tables[[pop]][[r]] <- load.and.combine.data(randomisation.ids[r], pop = pop, print.id = T)
#		h12.tables[[pop]]$true$diff[[r]] <- 0
#		for (chrom in c('2L', '2R', '3L', '3R', 'X')){
#			h12.tables[[pop]]$true$diff[chromosome == chrom][[r]] <- get.h12.interpolation(h12.tables[[pop]]$true$diff[chromosome == chrom], h12.tables[[pop]]$rando[[r]]$diff[chromosome == chrom], h12.column = 'H12.difference')
#		}
		#h12.tables[[pop]]$true$diff[, c(r) := get.h12.interpolation(.SD, 
		h12.tables[[pop]]$true$diff[is.peak == T, c(r) := get.h12.interpolation(.SD, 
		                                                                        ..h12.tables[[pop]][[r]]$diff[chromosome == .BY[[1]]], 
		                                                                        h12.column = 'H12.difference'
		                                                  ), 
		                                          by = chromosome,
		                                          .SDcols = c('chromosome', 'midpoint', 'H12.difference')
		]
	}
	cat('\tCalculating p-value of peaks based on randomised data.\n')
	h12.tables[[pop]]$true$diff[is.peak == T, pval := apply(.SD - H12.difference, 
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

# Now plot the data. Since we don't have all the chromosomes, we won't do the fancy plot with the chromosome
# underneath. Let's plot 2R and X separately.
chrom.sizes <- c('2R' = 61545105, '2L' = 49364325, '3R' = 53200684, '3L' = 41963435, 'X' = 24393108)
# Somewhat convoluted way of getting the min and max y value.
h12.limits <- list()
for (pop in study.pops){
	h12.limits[[pop]] <- h12.tables[[pop]] %>%
				         lapply(function(L) L$diff$H12.difference) %>%
				         unlist %>%
				         {c(min(.), max(.))}
}

pop == 'Obuasi.gambiae.Delta'
x11()
plot(c(1, chrom.sizes['2R']), h12.limits[[pop]], type = 'n', xlab = '2R position', ylab = 'H12')
# plot the randomised data in grey
sapply(h12.tables[[pop]], function(L) L$diff[chromosome == '2R', lines(midpoint, H12.difference, col = add.transparency('grey50', 0.7), lwd = 1)])
# add the true data
#true.h12.tables$diff[chromosome == '2R', lines(midpoint, H12.difference, col = add.transparency('red', 0.7), lwd = 3)]
###
#h12.tables[[pop]][[1]]$dead[chromosome == '2R', lines(midpoint, H12, col = add.transparency('red', 0.4), lwd = 3)]
#h12.tables[[pop]][[1]]$alive[chromosome == '2R', lines(midpoint, H12, col = add.transparency('green', 0.4), lwd = 3)]
h12.tables[[pop]][[1]]$diff[chromosome == '2R', lines(midpoint, H12.difference, col = add.transparency('green', 0.4), lwd = 3)]
h12.tables[[pop]][[1]]$diff[chromosome == '2R' & is.peak == T, points(midpoint, H12.difference, col = add.transparency('red', 0.4), pch = 19, cex = 0.7)]
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



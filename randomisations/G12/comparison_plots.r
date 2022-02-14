library(stringi)
library(data.table)
library(magrittr)

# Load all of the G12 tables
g12.filenames <- list.files('results_for_eric', '*.tsv', full.names = T)

# There are lots of randomisations that don't have 2L data. So let's just kick out the 2L ones for now
g12.filenames <- g12.filenames[!grepl('2L', g12.filenames)]

# For now, all the files are from Madina gambiae Delta, so we don't need to worry about splitting by dataset. 
# We just need to split by randomisation id. 
randomisation.ids <- unique(stri_extract_first_regex(g12.filenames, '(?<=\\.)\\d+(?=\\.)'))

# We will only work on randomisation ids 1 and 2 for now, since later ones have less data
randomisation.ids <- randomisation.ids[1:2]

# We have the problem that the genomic windows for G12 are not the same for the dead and the alive. There
# are a few ways in which we could get around this, and I don't think it will make much difference which 
# one we choose. 
# As a first attempt, we will do the following. First, we identify which of the dead / alive datasets has 
# the fewest windows. We take that dataset and, for each window, we identify the closest window from the 
# opposite dataset. And then we subtract dead from alive. Let's go. 
# First, we write a function that takes two vectors of table positions and, for each window in the first 
# vector, finds the nearest window in the second.
find.nearest.windows <- function(pos.1, pos.2){
	distances <- abs(outer(pos.1, pos.2, '-'))
	nearest.windows <- apply(distances, 1, which.min)
	nearest.windows
}

# Now a function that applies that to take the difference between nearest windows
get.g12.difference <- function(alive.table, dead.table){
	if (nrow(alive.table) < nrow(dead.table)){
		nearest.window.dead.g12 <- dead.table[find.nearest.windows(alive.table$midpoint, dead.table$midpoint), G12]
		diff.table <- alive.table[, .(chromosome = chromosome, midpoint = midpoint, G12.difference = G12 - ..nearest.window.dead.g12)]
	}
	else {
		nearest.window.alive.g12 <- alive.table[find.nearest.windows(dead.table$midpoint, alive.table$midpoint), G12]
		diff.table <- dead.table[, .(chromosome = chromosome, midpoint = midpoint, G12.difference = ..nearest.window.alive.g12 - G12)]
	}
	diff.table
}

# Write a function to load and combine all data for a given randomisation id
load.and.combine.data <- function(rand.id){
	alive.filenames <- grep(paste('alive.', rand.id, sep = ''), g12.filenames, value = T) %>%
	                   setNames(., stri_extract_first_regex(., '(?<=\\d\\.)[23LRX]+(?=\\.tsv)'))
	dead.filenames <- grep(paste('dead.', rand.id, sep = ''), g12.filenames, value = T) %>%
	                  setNames(., stri_extract_first_regex(., '(?<=\\d\\.)[23LRX]+(?=\\.tsv)'))
	
	if (any(names(dead.filenames) != names(alive.filenames)))
		stop(paste('Different chromosomes have been found for dead and alive samplesets for randomisation', rand.id))

	alive.data <- names(alive.filenames) %>%
	              lapply(function(chrom) fread(alive.filenames[chrom], sep = '\t')[, .(chromosome = ..chrom, midpoint = midpoint, G12 = G12)])
	dead.data <- names(dead.filenames) %>%
	             lapply(function(chrom) fread(dead.filenames[chrom], sep = '\t')[, .(chromosome = ..chrom, midpoint = midpoint, G12 = G12)])
	diff.data <- 1:length(alive.data) %>%
				 lapply(function(i) get.g12.difference(alive.data[[i]], dead.data[[i]]))
	
	list(alive = rbindlist(alive.data), dead = rbindlist(dead.data), diff = rbindlist(diff.data))
}
	
g12.tables <- lapply(randomisation.ids, load.and.combine.data)

# We now have a single difference for each randomisation. Let's add the real data
true.alive.filenames <- list.files('../../results/G12/tsv', 'Madina_gambiae.*Delta_alive', full.names = T) %>%
                        setNames(stri_extract_first_regex(., '(?<=gambiae_)[23LRX]+(?=_[DeltaPM]+_)'))
true.dead.filenames <- list.files('../../results/G12/tsv', 'Madina_gambiae.*Delta_dead', full.names = T) %>%
                        setNames(stri_extract_first_regex(., '(?<=gambiae_)[23LRX]+(?=_[DeltaPM]+_)'))
# Since we only had the data for 2R and X for the randomisations, we'll use the same for the true data for now.
true.alive.filenames <- grep('_[2RX]+_', true.alive.filenames, value = T)
true.dead.filenames <- grep('_[2RX]+_', true.dead.filenames, value = T)


true.alive.data <- names(true.alive.filenames) %>%
			       lapply(function(chrom) fread(true.alive.filenames[chrom], sep = '\t')[, .(chromosome = ..chrom, midpoint = midpoint, G12 = G12)])
true.dead.data <- names(true.dead.filenames) %>%
			      lapply(function(chrom) fread(true.dead.filenames[chrom], sep = '\t')[, .(chromosome = ..chrom, midpoint = midpoint, G12 = G12)])
true.diff.data <- 1:length(true.alive.data) %>%
			      lapply(function(i) get.g12.difference(true.alive.data[[i]], true.dead.data[[i]]))
true.g12.tables <- list(alive = rbindlist(true.alive.data), dead = rbindlist(true.dead.data), diff = rbindlist(true.diff.data))

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
g12.limits <- c(g12.tables, list(true.g12.tables)) %>%
              lapply(function(L) L$diff$G12.difference) %>%
              unlist %>%
              {c(min(.), max(.))}

x11()
plot(c(1, chrom.sizes['2R']), g12.limits, type = 'n', xlab = '2R position', ylab = 'G12')
# plot the randomised data in grey
sapply(g12.tables, function(L) L$diff[chromosome == '2R', lines(midpoint, G12.difference, col = 'grey50', lwd = 3)])
# add the true data
true.g12.tables$diff[chromosome == '2R', lines(midpoint, G12.difference, col = add.transparency('red', 0.7), lwd = 5)]

# Now the X
x11()
plot(c(1, chrom.sizes['X']), g12.limits, type = 'n', xlab = 'X position', ylab = 'G12')
# plot the randomised data in grey
sapply(g12.tables, function(L) L$diff[chromosome == 'X', lines(midpoint, G12.difference, col = 'grey50', lwd = 3)])
# add the true data
true.g12.tables$diff[chromosome == 'X', lines(midpoint, G12.difference, col = add.transparency('red', 0.7), lwd = 5)]






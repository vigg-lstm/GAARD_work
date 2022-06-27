library(data.table)
library(magrittr)
library(stringr)

# Load meta data
meta <- fread('gaard_metadata.tsv', sep = '\t')

# Load the karyotype data
kary.2La <- fread('../karyotypes/gaard_karyotypes.tsv', sep = '\t') %>%
            .[inversion == '2La', ] %>%
            .[, .(partner_sample_id, karyotype = round(mean_genotype))] %>%
            setkey(partner_sample_id)

# Load relatedness results
ngs.relate.output.fn <- 'king.csv'
ngs.relate <- fread(ngs.relate.output.fn, sep = '\t')
# Change the name of the king column for easier typing
colnames(ngs.relate)[3] <- tolower(colnames(ngs.relate)[3])

# Change numerical indices to sample names. Indices are 0-indexed, so add +1
ngs.relate$a <- meta[ngs.relate$a + 1, 'partner_sample_id']
ngs.relate$b <- meta[ngs.relate$b + 1, 'partner_sample_id']

# Now set the metadata table key to be the GAARD sample name. 
setkey(meta, 'partner_sample_id')

ngs.relate$species.a <- meta[ngs.relate$a, species_gambiae_coluzzii]
ngs.relate$species.b <- meta[ngs.relate$b, species_gambiae_coluzzii]
ngs.relate$location.a <- meta[ngs.relate$a, location]
ngs.relate$location.b <- meta[ngs.relate$b, location]
ngs.relate$insecticide.a <- meta[ngs.relate$a, insecticide]
ngs.relate$insecticide.b <- meta[ngs.relate$b, insecticide]
ngs.relate$phenotype.a <- meta[ngs.relate$a, phenotype]
ngs.relate$phenotype.b <- meta[ngs.relate$b, phenotype]
ngs.relate$karyotype.a <- kary.2La[ngs.relate$a, karyotype]
ngs.relate$karyotype.b <- kary.2La[ngs.relate$b, karyotype]

gamgam <- ngs.relate[species.a == 'gambiae' & species.b == 'gambiae', ]
colcol <- ngs.relate[species.a == 'coluzzii' & species.b == 'coluzzii', ]
gamcol <- ngs.relate[(species.a == 'gambiae' & species.b == 'coluzzii') | 
                     (species.a == 'coluzzii' & species.b == 'gambiae'), ]

# Overall histogram split by species (each species gets its own histogram)
png('king_histograms.png', width = 6, height = 2.2, units = 'in', res = 300)
par(mfrow = c(1,3), cex = 0.4, lwd = 0.6, mar = c(3,2,2,0), oma = c(0,1,0,0), mgp = c(2,0.8,0), xpd = NA)
hist(gamgam$king, ylim = c(0,25000), col = rgb(0,0,1,0.5), breaks = 100, main = 'all pairs', xlab = 'King')
hist(colcol$king, col = rgb(1,0,0,0.5), add = T, breaks = 50)
hist(gamcol$king, col =rgb(1,0,1,0.5), add = T, breaks = 50)
# Within gamgam, show the proportion of each bar that are from same location (unlike histogram
# above, here the coloured bar is a subset of the larger histogram
hist(gamgam$king, col = rgb(0,0,0.8), breaks = 100, main = 'gambiae pairs', xlab = 'King', ylab = '')
with(gamgam[location.a == location.b], hist(king, col = rgb(0.6,0.6,1), add = T, breaks = 100))
# And again for coluzzii
hist(colcol$king, col = rgb(0.8,0,0), breaks = 100, main = 'coluzzii pairs', xlab = 'King', ylab = '')
with(colcol[location.a == location.b], hist(king, col = rgb(1,0.6,0.6), add = T, breaks = 100))
dev.off()

# Let's draw the overall histogram, separately for each species pair type, colour-coded by karyotype
png('king_histograms_karyotype.png', width = 6, height = 2.2, units = 'in', res = 300)
par(mfrow = c(1,3), cex = 0.4, lwd = 0.6, mar = c(3,2,2,0), oma = c(0,1,0,0), mgp = c(2,0.8,0), xpd = NA)
hist(gamcol$king, col = rgb(0.8,0,0.8), breaks = 100, main = 'x-species', xlab = 'King')
with(gamcol[karyotype.a == karyotype.b], hist(king, col = rgb(0.8,0.6,0.8), add = T, breaks = 100))
# Within gamgam, show the proportion of each bar that are from same location (unlike histogram
# above, here the coloured bar is a subset of the larger histogram
hist(gamgam$king, col = rgb(0,0,0.8), breaks = 100, main = 'gambiae pairs', xlab = 'King', ylab = '')
with(gamgam[karyotype.a == karyotype.b], hist(king, col = rgb(0.6,0.6,1), add = T, breaks = 100))
# And again for coluzzii
hist(colcol$king, col = rgb(0.8,0,0), breaks = 100, main = 'coluzzii pairs', xlab = 'King', ylab = '')
with(colcol[karyotype.a == karyotype.b], hist(king, col = rgb(1,0.6,0.6), add = T, breaks = 100))
dev.off()


# Now zoom in to try and see sib clusters. 
png('king_histograms_zoomed.png', width = 4, height = 2.2, units = 'in', res = 300)
par(mfrow = c(1,2), cex = 0.4, lwd = 0.6, mar = c(3,2,2,0), oma = c(0,1,0,0), mgp = c(2,0.8,0))
hist(gamgam$king, col = rgb(0,0,0.8), breaks = 100, main = 'gambiae pairs', xlab = 'King', ylab = '', ylim = c(0,100))
mtext('frequency', 2, 0, T, cex = 0.45)
with(gamgam[location.a == location.b], hist(king, col = rgb(0.6,0.6,1), add = T, breaks = 100))
hist(colcol$king, col = rgb(0.8,0,0), breaks = 100, main = 'coluzzii pairs', xlab = 'King', ylab = '', ylim = c(0,100))
with(colcol[location.a == location.b], hist(king, col = rgb(1,0.6,0.6), add = T, breaks = 100))
dev.off()

# We still have a peak of negatively related individual
find.sib.groups <- function(sib.table, verbose = F){
	# First, create a list where each entry is a sib pair declared by the sib.table
	pair.list <- lapply(split(sib.table[, .(a, b)], 1:nrow(sib.table)), unlist)
	# Now merge any two entries in that list that share a sample
	for (i in length(pair.list):2){
		shared.samples <- which(sapply(pair.list[1:(i-1)], function(x) any(pair.list[[i]] %in% x)))
		if (length(shared.samples) >= 1){
			pair.list[[shared.samples[1]]] <- unique(c(pair.list[[shared.samples[1]]], pair.list[[i]]))
			pair.list[[i]] <- NULL
		}
	}
	sib.groups <- lapply(pair.list, function(group) sib.table[a %in% group | b %in% group, ])
	cat(length(sib.groups), 'sib groups were found.\n')
	
	cross.study.sibs.per.group <- sapply(sib.groups, function(sib) sum(sib$location.a != sib$location.b))
	population.bridge <- cross.study.sibs.per.group > 0
	if (any(population.bridge))
		cat('\t', sum(population.bridge), ' sib groups bridged across populations.')
	total.cross.study.pairs <- sum(cross.study.sibs.per.group)
	
	# Now go through and deal with groups that don't have coherent sets of sibs
	missing.pair.tables <- list()
	total.missing.pairs <- 0
	for (sib in sib.groups){
		members <- with(sib, unique(unlist(c(a, b))))
		expected.pairs <- combn(members, 2)
		num.missing <- ncol(expected.pairs) - nrow(sib)
		if (num.missing > 0){
			total.missing.pairs <- total.missing.pairs + num.missing
			# Get the table of missing pairs from the full table
			missing.pairs <- c()
			for (j in 1:ncol(expected.pairs)){
				this.pair <- expected.pairs[, j]
				pair.rows <- with(sib, (a %in% this.pair) & (b %in% this.pair))
				# If none of those values are true, then this is a missing pair and needs to be found in the full table
				if (!any(pair.rows))
					missing.pairs <- c(missing.pairs, which(with(ngs.relate, (a %in% this.pair) & (b %in% this.pair))))
			}
			missing.pair.table <- ngs.relate[missing.pairs, ]
			if (verbose)
				cat('\t', num.missing, ' missing pairs were found in group containing ', 
					paste(members, collapse = ', '), 
					' with min, mean and max KING scores of ', 
					paste(round(summary(missing.pair.table$king)[c('Min.', 'Mean', 'Max.')], 2), collapse  =', '),
					' respectively.\n', sep = '')
			missing.pair.tables <- c(missing.pair.tables, list(missing.pair.table))
		}
		else {
			if (verbose)
				cat('\tThere were no missing pairs.\n')
		}
	}
	total.pairs <- sum(sapply(sib.groups, nrow))
	list(sib.groups = sib.groups, 
	     total.pairs = total.pairs, 
	     total.missing.pairs = total.missing.pairs, 
		 missing.pair.ratio = total.missing.pairs / total.pairs,
	     total.cross.study.pairs = total.cross.study.pairs, 
	     missing.pair.tables = missing.pair.tables, 
	     cross.study.sibs.per.group = cross.study.sibs.per.group)
}

test.threshold <- function(king.thresh, pairs.table, verbose = F){
	cat('\nRunning find.sib.groups with a king threshold of', king.thresh, '\n')
	sib.pairs <- pairs.table[king >= king.thresh, ]
	invisible(find.sib.groups(sib.pairs, verbose = verbose))
}

thresholds <- seq(0.15, 0.35, 0.005)
threshold.exploration.gam <- lapply(thresholds, test.threshold, gamgam)
gam.num.cross.pairs <- sapply(threshold.exploration.gam, function(x) x$total.cross.study.pairs)
gam.total.missing.pairs <- sapply(threshold.exploration.gam, function(x) x$total.missing.pairs)
gam.total.pairs <- sapply(threshold.exploration.gam, function(x) x$total.pairs)
gam.missing.pair.ratio <- sapply(threshold.exploration.gam, function(x) x$missing.pair.ratio)
x11()
par(mfrow = c(1,2))
plot(thresholds, gam.total.pairs, pch = 19, col = 'black', ylim = c(0, max(c(gam.total.pairs, gam.total.missing.pairs))), main = 'gambiae pairs')
points(thresholds, gam.total.missing.pairs, pch = 19, col = 'blue')
points(thresholds, gam.num.cross.pairs, pch = 19, col = 'magenta')
plot(thresholds, gam.missing.pair.ratio, pch = 19, col = rgb(0.5,0.5,1))
thresholds[order(gam.missing.pair.ratio)]
#
threshold.exploration.col <- lapply(thresholds, test.threshold, colcol)
col.num.cross.pairs <- sapply(threshold.exploration.col, function(x) x$total.cross.study.pairs)
col.total.missing.pairs <- sapply(threshold.exploration.col, function(x) x$total.missing.pairs)
col.total.pairs <- sapply(threshold.exploration.col, function(x) x$total.pairs)
col.missing.pair.ratio <- sapply(threshold.exploration.col, function(x) x$missing.pair.ratio)
x11()
par(mfrow = c(1,2))
plot(thresholds, col.total.pairs, pch = 19, col = 'black', ylim = c(0, max(c(col.total.pairs, col.total.missing.pairs))), main = 'coluzzii pairs')
points(thresholds, col.total.missing.pairs, pch = 19, col = 'red')
points(thresholds, col.num.cross.pairs, pch = 19, col = 'magenta')
plot(thresholds, col.missing.pair.ratio, pch = 19, col = rgb(1,0.5,0.5))
thresholds[order(col.missing.pair.ratio)]
# What if we combine gamgam and colcol, so we just have a table of same-species pairs .
threshold.exploration <- lapply(thresholds, test.threshold, rbind(gamgam, colcol))
num.cross.pairs <- sapply(threshold.exploration, function(x) x$total.cross.study.pairs)
total.missing.pairs <- sapply(threshold.exploration, function(x) x$total.missing.pairs)
total.pairs <- sapply(threshold.exploration, function(x) x$total.pairs)
missing.pair.ratio <- sapply(threshold.exploration, function(x) x$missing.pair.ratio)
x11()
par(mfrow = c(1,2))
plot(thresholds, total.pairs, pch = 19, col = 'black', ylim = c(0, max(c(total.pairs, total.missing.pairs))), main = 'same species pairs')
points(thresholds, total.missing.pairs, pch = 19, col = 'red')
points(thresholds, num.cross.pairs, pch = 19, col = 'magenta')
plot(thresholds, missing.pair.ratio, pch = 19, col = rgb(1,0.5,0.5))
thresholds[order(missing.pair.ratio)]

# So we settle on a threshold of 0.195.
cat('We choose a king threshold of', thresholds[which.min(missing.pair.ratio)], 'as the best threshold\n.')
best.threshold.results <- threshold.exploration[[which.min(missing.pair.ratio)]]

# Create a list of sib groups, which includes the experiment that each individual was taken from
create.sib.group.table <- function(sib.pair.table, sib.group.id){
	sib.group.table <- data.table(sample.name = c(sib.pair.table$a, sib.pair.table$b), 
	                              location = c(sib.pair.table$location.a, sib.pair.table$location.b),
	                              species = c(sib.pair.table$species.a, sib.pair.table$species.b),
	                              insecticide = c(sib.pair.table$insecticide.a, sib.pair.table$insecticide.b),
	                              sib.group.id = sib.group.id
	                             ) %>%
	                   {.[!duplicated(.)]} %>%
	                   {.[order(sample.name)]}
	sib.group.table
}

sib.group.table <- mapply(create.sib.group.table, 
                          best.threshold.results$sib.groups, 
                          1:length(best.threshold.results$sib.groups),
                          SIMPLIFY = F) %>%
                   rbindlist()

set.seed(42)
sib.group.table[, keep := sample(c(T, rep(F, nrow(.SD) -1)), nrow(.SD)), by = c('location', 'insecticide', 'sib.group.id')]

# One sample in the table is a male. Remove it and it's sib
male.sibs <- intersect(sib.group.table$sample.name, meta[sex_call == 'M', partner_sample_id])
cat('\nThe following samples from the sib.group.table were male:\n')
cat('\n\t', male.sibs, '\n')
cat('\nThis sample had just one sib:\n')
male.sib.group.id <- sib.group.table[sample.name %in% male.sibs, sib.group.id]
print(setdiff(sib.group.table[sib.group.id == male.sib.group.id, sample.name], male.sibs))
cat('Removing this sib group from the table.\n')
sib.group.table <- sib.group.table[sib.group.id != male.sib.group.id]
# Adjust the sib group ids so as to not have a gap:
sib.group.table[sib.group.id > male.sib.group.id, sib.group.id := sib.group.id - 1]


cat('\nWe have removed', table(sib.group.table$keep)['FALSE'], 'samples as sibs, broken down as follows:\n\n')
meta[sib.group.table[!(keep), sample.name]] %>%
with(., table(location, phenotype, insecticide)) %>%
print

fwrite(sib.group.table, file = 'sib_group_table.csv', sep = '\t')


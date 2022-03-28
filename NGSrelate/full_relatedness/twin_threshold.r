library(data.table)

ngs.relate.output.fn <- 'ag3_gaard1244.3L.ngsRelate.tsv'
useful.columns <- c('specimen.x', 'location.x', 'species.x', 'specimen.y', 'location.y', 'species.y', 'KING')
ngs.relate <- fread(ngs.relate.output.fn, sep = '\t')[, ..useful.columns]
# Rename the column for each of typing
colnames(ngs.relate)[colnames(ngs.relate) == 'KING'] <- 'king'

find.twin.groups <- function(twin.table){
	# First, create a list where each entry is a twin pair declared by the twin.table
	pair.list <- lapply(split(twin.table[, .(specimen.x, specimen.y)], 1:nrow(twin.table)), unlist)
	# Now merge any two entries in that list that share a sample
	for (i in length(pair.list):2){
		shared.samples <- which(sapply(pair.list[1:(i-1)], function(x) any(pair.list[[i]] %in% x)))
		if (length(shared.samples) >= 1){
			pair.list[[shared.samples[1]]] <- unique(c(pair.list[[shared.samples[1]]], pair.list[[i]]))
			pair.list[[i]] <- NULL
		}
	}
	twin.groups <- lapply(pair.list, function(group) twin.table[specimen.x %in% group | specimen.y %in% group, ])
	cat(length(twin.groups), 'twin groups were found.\n')
	# Now go through and deal with groups that don't have coherent sets of twins
	cross.study.twin.groups <- sapply(twin.groups, function(tw) length(unique(c(tw$location.x, tw$location.y))))
	
	population.bridge <- cross.study.twin.groups > 1
	if (any(population.bridge))
		cat('\t', sum(population.bridge), ' twin groups bridged across populations.')
	
	missing.pair.tables <- list()
	for (tw in twin.groups){
		members <- with(tw, unique(unlist(c(specimen.x, specimen.y))))
		expected.pairs <- combn(members, 2)
		num.missing <- ncol(expected.pairs) - nrow(tw)
		if (num.missing > 0){
			# Get the table of missing pairs from the full table
			missing.pairs <- c()
			for (j in 1:ncol(expected.pairs)){
				this.pair <- expected.pairs[, j]
				pair.rows <- with(tw, (specimen.x %in% this.pair) & (specimen.y %in% this.pair))
				# If none of those values are true, then this is a missing pair and needs to be found in the full table
				if (!any(pair.rows))
					missing.pairs <- c(missing.pairs, which(with(ngs.relate, (specimen.x %in% this.pair) & (specimen.y %in% this.pair))))
			}
			missing.pair.tables <- c(missing.pair.tables, list(ngs.relate[missing.pairs, ]))
			cat('\t', num.missing, ' missing pairs were found in group containing ', paste(members, collapse = ','), 
			    ' with max, min and mean KING scores of ', summary(missing.pair.tables$king)[c('Min.', 'Mean', 'Max.')],
				' respectively.\n', sep = '')
		}
		else {
			cat('\tThere were no missing pairs.\n')
		}
	}
	list(twin.groups, missing.pair.tables, cross.study.twin.groups)
}

test.threshold <- function(king.thresh){
	cat('\nRunning find.twin.groups with a king threshold of', king.thresh, '\n')
	twin.pairs <- ngs.relate[king >= king.thresh, ]
	find.twin.groups(twin.pairs)
}
	

thresholds <- seq(0.3, 0.5, 0.02)
threshold.exploration <- lapply(thresholds, test.threshold)


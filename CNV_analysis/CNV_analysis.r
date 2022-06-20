library(data.table)
library(magrittr)
library(stringr)

study.ids.table <- fread('../data/study_ids.csv')[population != 'Aboisso_gambiae_PM']
study.ids.table[, study_id := paste('v3.2_', study_id, sep = '')]
meta <- fread('../data/combined/all_samples.samples.meta.csv', key = 'sample_id')
meta.sample.names <- meta$partner_sample_id
phen <- fread('../data/combined/sample_phenotypes.csv', key = 'specimen')[sort(meta.sample.names), ]

load.target.cnv.table <- function(study.id){
	target.cnv.table <- fread(paste('Ag1000G_CNV_data', study.id, 'target_regions_analysis/focal_region_CNV_table.csv', sep = '/'))
	these.sample.names <- meta[target.cnv.table$V1, partner_sample_id]
	target.cnv.table$sample.id <- these.sample.names
	target.cnv.table$V1 <- NULL
	target.cnv.table
}

target.CNV.table <- lapply(unique(study.ids.table$study_id), load.target.cnv.table) %>%
                    rbindlist() %>%
                    setkey(sample.id)

Dup.names <- grep('_D[upel]{2}', colnames(target.CNV.table), value = T)
Dup.families <- sub('_.*', '', Dup.names)
Dup.clusters <- split(Dup.names, Dup.families)

known.Dups <- Dup.names[!grepl('Dup0', Dup.names)]
sample.names <- target.CNV.table$sample.id
good.var.sample.names <- target.CNV.table[High.var.sample == F]$sample.id
# The "Any" category is where any of the known CNVs are present
Dups.any <- c('Ace1_Any', 'Cyp6aap_Any', 'Cyp6mz_Any', 'Gstue_Any', 'Cyp9k1_Any')
for (D in Dups.any){
	this.region <- sub('_Any', '', D)
	these.cnvs <- known.Dups[grepl(this.region, known.Dups)]
	target.CNV.table[[D]] <- apply(subset(target.CNV.table[, ..these.cnvs]), 1, any)
}

target.CNV.table$phenotype <- as.factor(phen[target.CNV.table$sample.id, phenotype])

# Get the CNV count per population
population.target.numCNVs <- aggregate(target.CNV.table[, c(Dup.names, Dups.any), with = F], 
                                       phen[sample.names, .(location, species)],
                                       sum
)
# Get the CNV frequency per population
population.target.CNVs <- aggregate(target.CNV.table[, c(Dup.names, Dups.any), with = F], 
                                    phen[sample.names, .(location, species)], 
                                    function(x) sum(x) / length(x)
)
# Get the size of each population-year combination
population.popsize <- interaction(phen[sample.names, .(location, species)]) %>%
                      droplevels %>%
                      table
# Now again using only samples with good variance
population.target.goodvar.numCNVs <- aggregate(target.CNV.table[good.var.sample.names, c(Dup.names, Dups.any), with = F], 
                                               phen[good.var.sample.names, .(location, species)],
                                               sum
)
population.target.goodvar.CNVs <- aggregate(target.CNV.table[good.var.sample.names, c(Dup.names, Dups.any), with = F], 
                                            phen[good.var.sample.names, .(location, species)], 
                                            function(x) sum(x) / length(x)
)
population.goodvar.popsize <- interaction(phen[good.var.sample.names, .(location, species)]) %>%
                      droplevels %>%
                      table

# Function to spin a table of by-population values and fill in missing genes with zeros.
spin.and.fill <- function(x, keyword) {
	out <- as.data.frame(t(x[, c(Dup.names, Dups.any)]))
	out[is.na(out)] <- 0
	colnames(out) <- paste(x$location, x$species, keyword, sep = '.')
	out
}

focal.region.cnv.frequencies <- cbind(spin.and.fill(population.target.numCNVs, 'all'),
                                      spin.and.fill(population.target.CNVs, 'all.prop'),
                                      spin.and.fill(population.target.goodvar.numCNVs, 'goodvar'),
                                      spin.and.fill(population.target.goodvar.CNVs, 'goodvar.prop')
)


# Now get the modal CNVs by gene
gene.table <- fread('Ag1000G_CNV_data/gene_annotation_fullgenetable.csv', key = 'Gene_stable_ID', check.names = T)

load.modal.cnv.table <- function(study.id){
	modal.cnv.table <- fread(paste('Ag1000G_CNV_data', study.id, 'modal_CNVs/modal_copy_number_gambcolu.csv', sep = '/'))
	these.sample.names <- meta[modal.cnv.table$V1, partner_sample_id]
	modal.cnv.table$V1 <- these.sample.names
	colnames(modal.cnv.table)[1] <- 'sample.id'
	modal.cnv.table
}

modal.copy.number <- lapply(unique(study.ids.table$study_id), load.modal.cnv.table) %>%
                     rbindlist() %>%
                     .[high_variance == F, -c('sex_call', 'high_variance')]

modal.CNV.table <- copy(modal.copy.number) %>%
                   .[, colnames(.)[-1] := lapply(.SD, `>`, 0), .SDcols = colnames(.)[-1]]

detox.genes <- c('Ace1',
                 paste('Cyp6aa', 1:2, sep = ''),
                 paste('Cyp6p', 1:5, sep = ''),
                 paste('Gste', 1:8, sep = ''),
                 'Cyp9k1')

detox.gene.conversion <- gene.table[Gene.name %in% toupper(detox.genes), 
                                    .(Gene.id = Gene_stable_ID, Gene.name = str_to_title(Gene.name))]

# We shrink the modal copy number table to be just the genes we are interested in. These are the detox 
# genes and also some of the ones in the Ace1 deletion in order to get the copy number of that. 
ace1.del.genes <- c('AGAP001358', 'AGAP001360', 'AGAP001361', 'AGAP001362', 'AGAP001363', 'AGAP001364', 'AGAP001365', 'AGAP001366')
genes.of.interest <- c(detox.gene.conversion$Gene.id, ace1.del.genes)
modal.copy.number <- modal.copy.number[, c('sample.id', genes.of.interest), with = F] %>%
                     merge(phen[good.var.sample.names],
                           by.x = 'sample.id', by.y = 'specimen',
                           all = F
                     ) %>%
                     .[, phenotype := as.factor(phenotype)] %>%
                     setnames(detox.gene.conversion$Gene.id, detox.gene.conversion$Gene.name) %>%
                     setkey(insecticide, location)

# For a smaller table, we only keep genes where a CNV was present. 
no.cnv <- names(which(apply(modal.CNV.table[, -1], 2, sum, na.rm = T) < 1))
modal.CNV.table <- modal.CNV.table[, -c(..no.cnv)]

# Calculate the frequency of CNVs in each population by year
population.modal.CNVs <- aggregate(modal.CNV.table[, -c('sample.id')],
                                   phen[modal.CNV.table$sample.id, .(location, species)], 
                                   function(x) sum(x) / length(x)
)

# A function that will lighten the colour by a given proportion
lighten.col <- function(color, lightness){
	col.rgb <- col2rgb(color)/255
	red <- 1-(1-col.rgb[1])*lightness
	blue <- 1-(1-col.rgb[2])*lightness
	green <- 1-(1-col.rgb[3])*lightness
	rgb(cbind(red, blue, green))
}

# A function that will round to x significant digits, unless the number is small, then just take
# one
signif2 <- function(x, digits, threshold = 0.1){
	small.x <- x < threshold 
	x[small.x] <- signif(x[small.x], 1)
	x[!small.x] <- signif(x[!small.x], digits)
	x
}


target.contable <- function(Dups, include.highvar = T, shrink = T, include.plot = T, 
                            add.sample.size = F, star.small.samples = T, 
                            species.colours = c(coluzzii = 'coral3', gambiae = 'dodgerblue3'), 
                            text.cell.cex = 0.65, pop.cex = 0.5, dup.cex = 0.5, ...){
	if (include.highvar)
		cnvs.by.population <- population.target.CNVs
	else 
		cnvs.by.population <- population.target.goodvar.CNVs
	population <- droplevels(interaction(cnvs.by.population[, c('location', 'species')]))
	output.table.rownames <- levels(population)
	output.table.colnames <- Dups
	output.table <- matrix(NA, length(output.table.rownames), length(Dups), dimnames = list(output.table.rownames, Dups))
	for (ro in output.table.rownames){
		which.row <- which(population == ro)
		if (length(which.row) == 1)
			output.table[ro, ] <- as.numeric(cnvs.by.population[which.row, Dups])
		else
			output.table[ro, ] <- NA
	}
	if (include.plot){
		if (include.highvar){
			cnvs.by.population <- population.target.CNVs
			popsize <- population.popsize[output.table.rownames]
		}
		else {
			cnvs.by.population <- population.target.goodvar.CNVs
			popsize <- population.goodvar.popsize[output.table.rownames]
		}
		#
		if ((length(names(species.colours)) == 0) | !any(names(species.colours) %in% c('arabiensis', 'gambiae', 'coluzzii', 'intermediate')))
			names(species.colours) <- c('arabiensis', 'gambiae', 'coluzzii', 'intermediate')
		par(mar = c(0,3,3,0), ...)
		plot(c(0,ncol(output.table)), -c(0,nrow(output.table)), type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
		x.matrix <- matrix(1:ncol(output.table), nrow(output.table), ncol(output.table), byrow = T)
		y.matrix <- matrix(1:nrow(output.table), nrow(output.table), ncol(output.table))
		col.matrix <- matrix(rgb(0.8, 0.8, 0.8), nrow(output.table), ncol(output.table))
		text.col.matrix <- c('black', 'white')[(output.table > 0.5) + 1]
		if (add.sample.size){
			text.cell.matrix <- sub('.*\\.', '.', paste(signif2(output.table, 2), ' (', popsize, ')', sep = ''))
			text.cell.matrix[grepl('NA', text.cell.matrix)] <- NA
		}
		else if (star.small.samples){
			text.cell.matrix <- sub('.*\\.', '.', paste(signif2(output.table, 2), c('', ' *', ' **')[(popsize < 5) + (popsize == 1) + 1], sep = ''))
			text.cell.matrix[grepl('NA', text.cell.matrix)] <- NA
		}
		else
			text.cell.matrix <- signif2(output.table, 2)
		for (i in 1:nrow(output.table)){
			this.row <- output.table[i, ]
			if (any(!is.na(this.row))){
				this.species <- sub('^.*\\.', '', rownames(output.table)[i])
				col.matrix[i, !is.na(this.row)] <- lighten.col(species.colours[this.species], sqrt(this.row[!is.na(this.row)]))
			}
		}
		rect(x.matrix - 1, 1 - y.matrix, x.matrix, - y.matrix, col = col.matrix, border = 'white', lwd = text.cell.cex)
		text(x.matrix - 0.5, 0.5 - y.matrix, text.cell.matrix, cex = text.cell.cex, col = text.col.matrix, font = 2)
		text(x.matrix[1,] - 0.5, (par('usr')[4] - par('usr')[3])/100, srt = 45, adj = c(0,0), colnames(output.table), xpd = NA, cex = dup.cex)
		country.names <- sub('\\.[^.]*$', '', rownames(output.table))
		species <- sub('.*\\.', '', rownames(output.table))
		rect(x.matrix[1,1]-1, 1 - y.matrix, par('usr')[1] - (par('usr')[2] - par('usr')[1])/5, -y.matrix, col = sapply(species.colours[species], lighten.col, 0.4), border = NA, xpd = NA)
		text(x.matrix[1,1]-1, 0.5 - y.matrix[,1], paste(country.names, '  \n', species, '  ', sep = ''), adj = 1, cex = pop.cex, xpd = NA)
	}
	if (shrink){
		# We keep the rows where at least one entry is > 0. We then remove the columns where all the remaining
		# rows are NA.
		which.rows <- apply(output.table, 1, sum, na.rm = T) > 0
		output.table <- output.table[which.rows, , drop = F]
		which.columns <- apply(output.table, 2, function(x) all(is.na(x))) == F
		output.table <- output.table[, which.columns, drop = F]
	}
	invisible(output.table)
}

modal.contable <- function(genes, shrink = T, include.plot = T, add.sample.size = F, 
                           star.small.samples = T, use.agaps = F, 
                           species.colours = c(coluzzii = 'coral3', gambiae = 'dodgerblue3'), 
                           text.cell.cex = 0.65, pop.cex = 0.5, gene.cex = 0.5, ...){
	cnvs.by.population <- population.modal.CNVs
	#
	if (genes[1] %in% gene.table$Gene_stable_ID)
		gene.names <- gene.table[genes, Gene.name]
	else {
		gene.names <- character()
		for (i in 1:length(genes)){
			gene <- genes[i]
			if (toupper(gene) %in% toupper(gene.table$Gene.name)){
				gene.name <- gene.table$Gene.name[toupper(gene.table$Gene.name) == toupper(gene)]
				if (length(gene.name) != 1) 
					stop('There should be only one gene matching the gene name.')
				gene.names[i] <- str_to_title(gene.name)
				genes[i] <- gene.table[Gene.name == gene.name, Gene_stable_ID]
			}
			else
				stop('Could not find the requested gene in gene.table.')
		}
	}
	population <- droplevels(interaction(cnvs.by.population[, c('location', 'species')]))
	output.table.rownames <- levels(population)
	output.table.colnames <- if(use.agaps) genes else gene.names
	output.table <- matrix(NA, length(output.table.rownames), length(output.table.colnames), dimnames = list(output.table.rownames, output.table.colnames))
	gene.has.cnv <- genes %in% colnames(cnvs.by.population)
	for (ro in output.table.rownames){
		which.row <- which(population == ro)
		if (length(which.row) != 1)
			output.table[ro, ] <- NA
		else {
			for (i in 1:ncol(output.table)){
				co <- output.table.colnames[i]
				if (gene.has.cnv[i])
					output.table[ro, co] <- as.numeric(cnvs.by.population[which.row, genes[i]])
				else
					output.table[ro, co] <- 0
			}
		}
	}
	if (include.plot){
		if (is.null(names(species.colours)) | !any(names(species.colours) %in% c('gambiae', 'coluzzii')))
			names(species.colours) <- c('gambiae', 'coluzzii')
		par(mar = c(0,3,3,0), ...)
		plot(c(0,ncol(output.table)), -c(0,nrow(output.table)), type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
		x.matrix <- matrix(1:ncol(output.table), nrow(output.table), ncol(output.table), byrow = T)
		y.matrix <- matrix(1:nrow(output.table), nrow(output.table), ncol(output.table))
		col.matrix <- matrix(rgb(0.8, 0.8, 0.8), nrow(output.table), ncol(output.table))
		popsize <- population.goodvar.popsize[output.table.rownames]
		text.col.matrix <- c('black', 'white')[(output.table > 0.5) + 1]
		if (add.sample.size){
			text.cell.matrix <- sub('.*\\.', '.', paste(signif2(output.table, 2), ' (', popsize, ')', sep = ''))
			text.cell.matrix[grepl('NA', text.cell.matrix)] <- NA
		}
		else if (star.small.samples){
			text.cell.matrix <- sub('.*\\.', '.', paste(signif2(output.table, 2), c('', ' *', ' **')[(popsize < 5) + (popsize == 1) + 1], sep = ''))
			text.cell.matrix[grepl('NA', text.cell.matrix)] <- NA
		}
		else {
			text.cell.matrix <- signif2(output.table, 2)
		}
		for (i in 1:nrow(output.table)){
			this.row <- output.table[i, ]
			if (any(!is.na(this.row))){
				this.species <- sub('^.*\\.', '', rownames(output.table)[i])
				col.matrix[i, !is.na(this.row)] <- lighten.col(species.colours[this.species], sqrt(this.row[!is.na(this.row)]))
			}
		}
		rect(x.matrix - 1, 1 - y.matrix, x.matrix, - y.matrix, col = col.matrix, border = 'white', lwd = text.cell.cex)
		text(x.matrix - 0.5, 0.5 - y.matrix, text.cell.matrix, cex = text.cell.cex, col = text.col.matrix, font = 2)
		text(x.matrix[1,] - 0.5, (par('usr')[4] - par('usr')[3])/100, srt = 45, adj = c(0,0), colnames(output.table), xpd = NA, cex = gene.cex)
		country.names <- sub('\\.[^.]*$', '', rownames(output.table))
		species <- sub('.*\\.', '', rownames(output.table))
		rect(x.matrix[1,1]-1, 1 - y.matrix, par('usr')[1] - (par('usr')[2] - par('usr')[1])/5, -y.matrix, col = sapply(species.colours[species], lighten.col, 0.4), border = NA, xpd = NA)
		text(x.matrix[1,1]-1, 0.5 - y.matrix[,1], paste(country.names, '  \n', species, '  ', sep = ''), adj = 1, cex = pop.cex, xpd = NA)
	}
	if (shrink){
		# We keep the rows where at least one entry is > 0. We then remove the columns where all the remaining
		# rows are NA.
		which.rows <- apply(output.table, 1, sum, na.rm = T) > 0
		output.table <- output.table[which.rows, , drop = F]
		which.columns <- apply(output.table, 2, function(x) all(is.na(x))) == F
		output.table <- output.table[, which.columns, drop = F]
	}
	invisible(output.table)
}


contable <- function(x, ...){
	if (grepl('_Dup', x[1]) | grepl('_Del', x[1]) | grepl('_Any', x[1]))
		target.contable(x, ...)
	else
		modal.contable(x, ...)
}

png('detox_gene_consensus_CNVs.png', width = 4, height = 1, units = 'in', res = 300)
par(family = 'Arial')
contable(detox.genes, 
         text.cell.cex = 0.35,
         pop.cex = 0.35,
         gene.cex = 0.35,
         mai = c(0,0.17,0.18,0)
)
dev.off()

png('Ace1_diagnostic_read_CNVs.png', width = 2, height = 1.3, units = 'in', res = 300)
par(family = 'Arial')
contable(Dup.clusters$Ace1, 
         text.cell.cex = 0.35,
         pop.cex = 0.35,
         dup.cex = 0.35,
         mai = c(0,0.23,0.23,0.05)
)
dev.off()

png('Cyp6aap_diagnostic_read_CNVs.png', width = 6, height = 1.4, units = 'in', res = 300)
par(family = 'Arial')
contable(Dup.clusters$Cyp6aap, 
         text.cell.cex = 0.35,
         pop.cex = 0.35,
         dup.cex = 0.35,
         mai = c(0,0.09,0.33,0.05)
)
dev.off()

png('Cyp6mz_diagnostic_read_CNVs.png', width = 2.5, height = 1.4, units = 'in', res = 300)
par(family = 'Arial')
contable(Dup.clusters$Cyp6mz, 
         text.cell.cex = 0.35,
         pop.cex = 0.35,
         dup.cex = 0.35,
         mai = c(0,0.21,0.33,0.18)
)
dev.off()

png('Gste_diagnostic_read_CNVs.png', width = 3, height = 1.4, units = 'in', res = 300)
par(family = 'Arial')
contable(Dup.clusters$Gstue, 
         text.cell.cex = 0.35,
         pop.cex = 0.35,
         dup.cex = 0.35,
         mai = c(0,0.21,0.3,0.18)
)
dev.off()

png('Cyp9k1_diagnostic_read_CNVs.png', width = 6, height = 1.4, units = 'in', res = 300)
par(family = 'Arial')
contable(Dup.clusters$Cyp9k1, 
         text.cell.cex = 0.35,
         pop.cex = 0.35,
         dup.cex = 0.35,
         mai = c(0,0.09,0.33,0.05)
)
dev.off()


# Now correlate modal copy number in cyp9k1 with phenotype
cat('\nModal copy number in Cyp9k1:\n')
for (insecticide in c('Delta', 'PM')){
	for (location in c('Avrankou', 'Baguida', 'Korle-Bu', 'Madina', 'Obuasi')){
		if (insecticide == 'PM' & location == 'Avrankou')
			next
		cat('\n\n\t', location, insecticide, 'Cyp9k1:\n')
		# For some reason, this doesn't work if I do the indexing directly. I have to create this object first
		index <- list(insecticide, location)
		print(modal.copy.number[index, table(phenotype, Cyp9k1)])
		print(modal.copy.number[index, drop1(glm(phenotype ~ Cyp9k1, family = 'binomial'), test = 'Chisq')])
	}
}

# None of those tests are significant.

# Let's try a glm that includes population as a fixed factor. We don't also include species, since the 
# effect of species will be perfectly captured by the population fixed factor (all populations contained
# only one species). 
cat('\n\nCyp9k1 combined Deltamethrin, location as fixed factor:\n')
print(modal.copy.number['Delta', drop1(glm(phenotype ~ Cyp9k1 + location, family = 'binomial'), test = 'Chisq')])
# That's non-significant. 
cat('\n\nCyp9k1 combined PM, location as fixed factor:\n')
print(modal.copy.number['PM', drop1(glm(phenotype ~ Cyp9k1 + location, family = 'binomial'), test = 'Chisq')])
# Also non-significant.

# Now again with Gste2. We only find dups at appreciable frequency in the coluzzii populations
cat('\nModal copy number in Gste2:\n')
for (insecticide in c('Delta', 'PM')){
	for (location in c('Avrankou', 'Korle-Bu')){
		if (insecticide == 'PM' & location == 'Avrankou')
			next
		cat('\n\n\t', location, insecticide, 'Gste2:\n')
		index <- list(insecticide, location)
		print(modal.copy.number[index, table(phenotype, Gste2)])
		print(modal.copy.number[index, drop1(glm(phenotype ~ Gste2, family = 'binomial'), test = 'Chisq')])
	}
}
# Nothing significant

# Cyp6aa1. Again, only at decent numbers in coluzzii
cat('\nModal copy number in Cyp6aa1:\n')
for (insecticide in c('Delta', 'PM')){
	for (location in c('Avrankou', 'Korle-Bu')){
		if (insecticide == 'PM' & location == 'Avrankou')
			next
		cat('\n\n\t', location, insecticide, 'Cyp6aa1:\n')
		index <- list(insecticide, location)
		print(modal.copy.number[index, table(phenotype, Cyp6aa1)])
		print(modal.copy.number[index, drop1(glm(phenotype ~ Cyp6aa1, family = 'binomial'), test = 'Chisq')])
	}
}
# Korle-Bu Delta significant. Avrankou Delta nearly significant

# Cyp6p3. Found at appreciable freqs in Korle-Bu and Madina
cat('\nModal copy number in Cyp6p3:\n')
for (insecticide in c('Delta', 'PM')){
	for (location in c('Korle-Bu', 'Madina')){
		cat('\n\n\t', location, insecticide, 'Cyp6p3:\n')
		index <- list(insecticide, location)
		print(modal.copy.number[index, table(phenotype, Cyp6p3)])
		print(modal.copy.number[index, drop1(glm(phenotype ~ Cyp6p3, family = 'binomial'), test = 'Chisq')])
	}
}
# Nothing significant

# Combine populations
cat('\n\nCyp6aa1 + Cyp6p3 combined Deltamethrin, location as fixed factor:\n')
print(modal.copy.number[.('Delta', c('Avrankou', 'Korle-Bu', 'Madina')), drop1(glm(phenotype ~ Cyp6aa1 + Cyp6p3 + location, family = 'binomial'), test = 'Chisq')])
# P3 non-significant, so drop it from the model
print(modal.copy.number[.('Delta', c('Avrankou', 'Korle-Bu', 'Madina')), drop1(glm(phenotype ~ Cyp6aa1 + location, family = 'binomial'), test = 'Chisq')])
# What's the direction?
cat('\n\t Coefficients (negative means positively association with resistance:\n')
print(modal.copy.number[.('Delta', c('Avrankou', 'Korle-Bu', 'Madina')), glm(phenotype ~ Cyp6aa1 + location, family = 'binomial')])
# That's significant for Cyp6aa1. 

cat('\n\nCyp6aa1 + Cyp6p3 combined PM, location as fixed factor:\n')
print(modal.copy.number[.('PM', c('Avrankou', 'Korle-Bu', 'Madina')), drop1(glm(phenotype ~ Cyp6aa1 + Cyp6p3 + location, family = 'binomial'), test = 'Chisq')])
# P3 least significant, so drop it from the model
print(modal.copy.number[.('PM', c('Avrankou', 'Korle-Bu', 'Madina')), drop1(glm(phenotype ~ Cyp6aa1 + location, family = 'binomial'), test = 'Chisq')])
# Non-significant

# Ace1.
cat('\nModal copy number in Ace1:\n')
for (insecticide in c('Delta', 'PM')){
	for (location in c('Baguida', 'Korle-Bu', 'Madina', 'Obuasi')){
		cat('\n\n\t', location, insecticide, 'Ace1:\n')
		index <- list(insecticide, location)
		print(modal.copy.number[index, table(phenotype, Ace1)])
		print(modal.copy.number[index, drop1(glm(phenotype ~ Ace1, family = 'binomial'), test = 'Chisq')])
	}
}
# PM significant in KB, Madina and Obuasi, but not in Baguida, despite variation in copy number. In fact, copy
# number overall is quite high in Baguida, and yet there were appreciable numbers of dead with high copy number
# (eg: 6). Exposure wasn't even that high (0.5x and changing exposure time). 

# Combine populations
cat('\n\nAce1 combined Deltamethrin, location as fixed factor:\n')
print(modal.copy.number['Delta', drop1(glm(phenotype ~ Ace1 + location, family = 'binomial'), test = 'Chisq')])
# That's non-significant for Delta, as expected
cat('\n\nAce1 combined PM, location as fixed factor:\n')
print(modal.copy.number['PM', drop1(glm(phenotype ~ Ace1 + location, family = 'binomial'), test = 'Chisq')])
# And highly significant for PM. 

# Add the presence of the Del to the model? Presence of the Del is of course correlated with copy number. But
# does it add anything to the model? 
modal.copy.number$Ace1_Del1_presence <- target.CNV.table[modal.copy.number$sample.id, Ace1_Del1]
cat('\n\nAce1 combined + Del1 presence, PM, location as fixed factor:\n')
print(modal.copy.number['PM', drop1(glm(phenotype ~ Ace1 + Ace1_Del1_presence + location, family = 'binomial'), test = 'Chisq')])
cat('\n\t Coefficients (negative means positively association with resistance:\n')
print(modal.copy.number['PM', glm(phenotype ~ Ace1 + Ace1_Del1_presence + location, family = 'binomial')])
# Having the deletion increases the level of resistance. 
# So even once the number of copies of Ace1 is taken into account, we are still seeing a significant effect of the
# number of deletions. 

# An alternative model would be one that counts the number of deletions. We can do this by taking the mean modal
# copy number of several genes present in the deletion and subtracting this from the copy number in Ace1. That's 
# the number of deletions. We leave out AGAP001357 because it's close to the breakpoint.
modal.copy.number[, Ace1_Del_cn := Ace1 - apply(.SD, 1, median), .SDcols = ace1.del.genes]
# There is a mismatch with the presence of Del1 because of the other deletions. So some samples have reduced 
# copy number and no Del1. 

cat('\n\nAce1 combined + Del1 copy number, PM, location as fixed factor:\n')
print(modal.copy.number['PM', drop1(glm(phenotype ~ Ace1 + Ace1_Del_cn + location, family = 'binomial'), test = 'Chisq')])

# It's now non-significant. It's at least partly because in the samples with a deletion copy number of 1, survival
# is higher in those that have Del1 than those that don't. 

# For Cyp6aap, let's look at each Dup in turn. We'll do that by population because there are different CNVs in
# different populations. In gambiae, there is very little diversity of Dups, so not really worth looking into 
# individual CNVs. 
# Avrankou:
cat('\n\nCyp6aap alleles in Avrankou:\n')
avrankou.delta.samples <- phen[location == 'Avrankou', specimen]
print(target.CNV.table[avrankou.delta.samples, drop1(glm(phenotype ~ Cyp6aap_Dup7 + Cyp6aap_Dup10, family = 'binomial'), test = 'Chisq')])
# Dup10 is least significant
print(target.CNV.table[avrankou.delta.samples, drop1(glm(phenotype ~ Cyp6aap_Dup7, family = 'binomial'), test = 'Chisq')])
# Dup7 also non significant 

# Korle-Bu:
korlebu.delta.samples <- phen[location == 'Korle-Bu' & insecticide == 'Delta', specimen]
print(target.CNV.table[korlebu.delta.samples, drop1(glm(phenotype ~ Cyp6aap_Dup10 + Cyp6aap_Dup17 + Cyp6aap_Dup18, family = 'binomial'), test = 'Chisq')])
# Dup17 is least significant
print(target.CNV.table[korlebu.delta.samples, drop1(glm(phenotype ~ Cyp6aap_Dup10 + Cyp6aap_Dup18, family = 'binomial'), test = 'Chisq')])
# Dup18 is non-significant. 
print(target.CNV.table[korlebu.delta.samples, drop1(glm(phenotype ~ Cyp6aap_Dup10, family = 'binomial'), test = 'Chisq')])
# Dup10 is significant
cat('\n\t Coefficients (negative means positively association with resistance:\n')
print(target.CNV.table[korlebu.delta.samples, glm(phenotype ~ Cyp6aap_Dup10 + Cyp6aap_Dup17 + Cyp6aap_Dup18, family = 'binomial')])
# Both Dup10 and Dup18 are positively associated with resistance, although Dup18 is non-sig. 



library(data.table)
library(magrittr)
library(stringr)
library(ape)

cyp9k1.region <- c(14000000, 17000000)

locations <- c('Avrankou', 'Baguida', 'Korle-Bu', 'Madina', 'Obuasi')

fst <- list()
for (loc in locations){
	load(paste('../../randomisations/Fst/fst_randomisations_', loc, '.Rdata', sep = ''))
	for (pop in names(fst.tables)){
		fst[[pop]] <- cbind(window.pos[[pop]], fst.tables[[pop]]) %>%
		              .[window.chrom == 'X' & window.pos > cyp9k1.region[1] & window.pos < cyp9k1.region[2]] %>%
		              .[, .(pos = window.pos, Fst = phenotype.moving.Fst, is.sig = F)]
	}
}

keep <- c('fst', 'cyp9k1.region')
rm(list = setdiff(ls(), keep))

fst.regions.fn <- '../../haplotypes/haplotype_significance_tests.csv'
fst.regions <- fread(fst.regions.fn) %>%
               .[, c('chrom', 'start', 'end') := tstrsplit(window, ':|-')] %>%
               .[(direction == 'alive' & logregP < 0.05)] %>%
               .[chrom == 'X' & as.numeric(end) > cyp9k1.region[1] & as.numeric(start) < cyp9k1.region[2]] %>%
               .[, .(population, start = as.numeric(start), end = as.numeric(end), pval = logregP)] %>%
               unique() 

for (i in seq_len(nrow(fst.regions))){
	pop <- fst.regions[i, population]
	start.point <- fst.regions[i, start]
	end.point <- fst.regions[i, end]
	fst[[pop]][pos > start.point & pos < end.point, is.sig := T]
}



# Get the list of genes present in the cyp9k1 region
gff.path <- '~/data/ML/GAARD_SNP/data/Anopheles-gambiae-PEST_BASEFEATURES_AgamP4.9.gff3'
gff <- as.data.table(read.gff(gff.path, GFF3 = T))
gene.gff <- gff[seqid == 'X' & 
				start < cyp9k1.region[2] & 
				end > cyp9k1.region[1] &
				type == 'gene', ]
gene.gff[, gene.id := unique(str_extract(attributes, 'AGAP\\d{6}'))]
gene.gff[, gene.name := str_extract(attributes, '(?<=Name=)[^;]+') %>%
						str_to_title]
gene.gff$gene.name[is.na(gene.gff$gene.name)] <- gene.gff$gene.id[is.na(gene.gff$gene.name)]
# Add the name of NADH dehydrogenase subunit 1 (it was labelled as "description" rather than "name"
# in the gff, which is why it got missed. 
gene.gff[gene.id == 'AGAP000849', gene.name := 'NADH dehydrogenase\n(ubiquinone) 1\nbeta subcomplex 1']
#
exon.gff <- gff[seqid == 'X' & 
				start < cyp9k1.region[2] & 
				end > cyp9k1.region[1] &
	            type == 'exon', ]

h12.fn <- '../../randomisations/H12/h12_filtered_windows.RDS'
h12 <- readRDS(h12.fn) %>%
       lapply(function(D) D$diff[chromosome == 'X' & midpoint > cyp9k1.region[1] & midpoint < cyp9k1.region[2], 
                                 .(pos = midpoint, start = startpoint, end = endpoint, H12 = H12, is.peak, pval)]
       )
names(h12) <- gsub('\\.', '_', names(h12))

pbs.fn <- '../../randomisations/PBS/pbs_filtered_windows.RDS'
pbs <- readRDS(pbs.fn) %>%
       lapply(function(D) D[chromosome == 'X' & midpoint > cyp9k1.region[1] & midpoint < cyp9k1.region[2], 
                            .(pos = midpoint, start = startpoint, end = endpoint, PBS = PBS, is.peak, pval)]
       )
names(pbs) <- gsub('\\.', '_', names(pbs))

source('../../shared_functions/R_plotting.r')

# Right, now let's make the composite plots. 
composite.plot <- function(pop, filename = NULL){
	fst.table <- fst[[pop]]
	h12.table <- h12[[pop]]
	pbs.table <- pbs[[pop]]
	if (!is.null(filename)){
		file.width = 6.5
		file.height = 4.5
		if (grepl('\\.eps', filename))
			postscript(filename, width = file.width, height = file.height, horizontal = F, onefile = FALSE, paper = "special")
		else if (grepl('\\.png', filename))
			png(filename, width = file.width, height = file.height, units = 'in', res = 300)
		else if (grepl('\\.tif', filename))
			tiff(filename, width = file.width, height = file.height, units = 'in', res = 300)
	}
	
	par(mar = c(0,7,1,4), mgp = c(2, 0.7, 0), family = 'Arial', xpd = NA) 
	colours <- c(pbs = 'royalblue1',
	             fst = 'orangered3', 
	             h12 = 'limegreen')
	# Create the empty plot
	max.y <- max(pbs.table$PBS, fst.table$Fst, h12.table$H12)
	min.y <- min(pbs.table$PBS, fst.table$Fst, h12.table$H12)
	min.plot.y <- min.y - (max.y - min.y)/5
	plot(c(cyp9k1.region[1], cyp9k1.region[2]), 
	     c(min.plot.y, max.y), 
	     type = 'n', bty = 'n', 
	     xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', 
	     main = pop, cex.main = 1.5
	)

	max.fst <- max(c(max(fst.table$Fst, na.rm = T), 0.05))
	fst.table$normed.fst <- max.y * fst.table$Fst / max.fst
	max.h12 <- max(c(max(h12.table$H12, na.rm = T), 0.05))
	h12.table$normed.h12 <- max.y * h12.table$H12 / max.h12
	max.pbs <- max(c(max(pbs.table$PBS, na.rm = T), 0.05))
	pbs.table$normed.pbs <- max.y * pbs.table$PBS / max.pbs

	# Add data
	fst.table[, lines(pos, normed.fst, col = colours['fst'], lwd = 1.5)]
	h12.table[, lines(unlist(transpose(.(start, end))), rep(normed.h12, each = 2), col = colours['h12'], lwd = 1.5)]
	#h12.table[, lines(pos, normed.h12, col = colours['h12'], lwd = 1.5)]
	pbs.table[, lines(unlist(transpose(.(start, end))), rep(normed.pbs, each = 2), col = colours['pbs'], lwd = 1.5)]

	# Highlight the significant regions
	usr <- par('usr')
	bead.coords <- usr[3] + (usr[4] - usr[3])/10
	h12.table[is.peak == T, text(apply(.SD[, .(start, end)], 1, median), y = normed.h12, '*', cex = 2, col = colours['h12'])]
	h12.table[is.peak == T & pval < 0.01][, rect(start, bead.coords, end, normed.h12, col = lighten.col(colours['h12'], 1, alpha = 0.2), border = NA)]
	pbs.table[is.peak == T, text(apply(.SD[, .(start, end)], 1, median), y = normed.pbs, '*', cex = 2, col = colours['pbs'])]
	pbs.table[is.peak == T & pval < 0.01][, rect(start, bead.coords, end, normed.pbs, col = lighten.col(colours['pbs'], 1, alpha = 0.2), border = NA)]

	# Add y axes
	fst.step.size <- ifelse(max.fst > 0.2, 0.1, 0.05)
	axis(2, at = seq(0, max.y, max.y * fst.step.size/max.fst), 
	     labels = as.character(seq(0, max.fst, fst.step.size)), 
	     col = colours['fst'], col.axis = colours['fst'], 
	     col.ticks = colours['fst'], lwd = 1.5
	)
	mtext(expression(F[ST]) , 2, 2, cex = 0.8, col = colours['fst'], font = 2)
	#
	h12.step.size <- ifelse(max.h12 > 0.2, 0.1, 0.05)
	axis(4, at = seq(0, max.y, max.y * h12.step.size/max.h12), 
	     labels = as.character(seq(0, max.h12, h12.step.size)), 
	     col = colours['h12'], col.axis = colours['h12'], 
	     col.ticks = colours['h12'], lwd = 1.5
	)
	mtext(expression(paste(Delta, 'H'[12])), 4, 2, cex = 0.8, col = colours['h12'], font = 2)
	#
	pbs.step.size <- ifelse(max.pbs > 0.2, 0.1, 0.05)
	axis(2, at = seq(0, max.y, max.y * pbs.step.size/max.pbs), 
	     labels = as.character(seq(0, max.pbs, pbs.step.size)), 
	     col = colours['pbs'], col.axis = colours['pbs'], 
	     col.ticks = colours['pbs'], line = 4, lwd = 1.5
	)
	mtext('PBS', 2, 6, cex = 0.8, col = colours['pbs'], font = 2)
	
	#par(mar = c(0,7,0,4), mgp = c(2, 0.7, 0), family = 'Arial', xpd = NA) 
#	plot(c(cyp9k1.region[1], cyp9k1.region[2]), 
#	     c(0,0), 
#	     type = 'n', bty = 'n', 
#	     xaxt = 'n', yaxt = 'n', xlab = '', ylab = ''
#	)

	draw.gene.model(
		cyp9k1.region[1], cyp9k1.region[2], 
		gene.gff, exon.gff, include.gene.names = F,
		y = bead.coords, gene.thickness.fraction = 100
	)
	# Now just draw the genes of interest again to get the name labels
	genes.to.label <- c('AGAP000818', 'AGAP000877', 'AGAP000849')
	draw.gene.model(
		cyp9k1.region[1], cyp9k1.region[2], 
		gene.gff[gene.id %in% genes.to.label], exon.gff[numeric(), ], 
		text.cex = 1, y = bead.coords, gene.thickness.fraction = 100
	)
	







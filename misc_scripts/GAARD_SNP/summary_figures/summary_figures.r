library(data.table)
library(magrittr)

sample.sets <- setNames(nm = c('Avrankou_coluzzii_Delta',
                               'Baguida_gambiae_Delta',
                               'Korle-Bu_coluzzii_Delta',
                               'Madina_gambiae_Delta',
                               'Obuasi_gambiae_Delta',
                               'Baguida_gambiae_PM',
                               'Korle-Bu_coluzzii_PM',
                               'Madina_gambiae_PM',
                               'Obuasi_gambiae_PM'))

chrom.sizes <- c('2R' = 61545105, '2L' = 49364325, '3R' = 53200684, '3L' = 41963435, 'X' = 24393108)

# Some plotting functions

# Write a function to draw the chromosomes on an existing plotting device
add.chroms <- function(gaps, cs, ce, gene.cex = 0.9, gene.col = 'grey20', point.cex = 1.2, point.col = 'grey30', chrom.col = NULL, chrom.cex = 1.4, chrom.offset = 0){
	plot(c(cs[1], ce[5]), c(-6.5,1.3), xlim = c(cs[1] + gaps/2, ce[5] - gaps/2), type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
	# show the kdr region
	kdr.region.mean <- cs['2L'] + mean(c(2358158, 2431617))
	points(kdr.region.mean, 0, pch = 19, cex = point.cex, col = point.col)
	text(kdr.region.mean, -1.7, 'Vgsc', srt = 45, adj = 1.1, xpd = T, cex = gene.cex, col = gene.col, font = 2)
	# show the Rdl region
	rdl.region.mean <- cs['2L'] + mean(c(25363652, 25434556))
	points(rdl.region.mean, 0, pch = 19, cex = point.cex, col = point.col)
	text(rdl.region.mean, -1.7, 'Rdl', srt = 45, adj = 1.1, xpd = T, cex = gene.cex, col = gene.col, font = 2)
	# show the Carboxylesterase region
	coeae.region.mean <- cs['2L'] + mean(c(28548433, 28550748))
	points(coeae.region.mean, 0, pch = 19, cex = point.cex, col = point.col)
	text(coeae.region.mean, -1.7, 'Coeae2f', srt = 45, adj = 1.1, xpd = T, cex = gene.cex, col = gene.col, font = 2)
	# show the Ace1 gene region
	ace1.region.mean <- cs['2R'] + mean(c(3484107, 3495790))
	points(ace1.region.mean, 0, pch = 19, cex = point.cex, col = point.col)
	text(ace1.region.mean, -1.7, 'Ace1', srt = 45, adj = 1.1, xpd = T, cex = gene.cex, col = gene.col, font = 2)
	# show the CYP6 region
	cyp6.region.mean <- cs['2R'] + mean(c(28463000, 28568000))
	points(cyp6.region.mean, 0, pch = 19, cex = point.cex, col = point.col)
	text(cyp6.region.mean, -1.7, 'Cyp6aa1-\nCyp6p2 ', srt = 45, adj = 1.1, xpd = T, cex = gene.cex, col = gene.col, font = 2)
	# show the Keap1 region (AGAP003645)
	keap1.region.mean <- cs['2R'] + mean(c(40927085, 40930289))
	points(keap1.region.mean, 0, pch = 19, cex = point.cex, col = point.col)
	text(keap1.region.mean, -1.7, 'Keap1', srt = 45, adj = 1, xpd = NA, cex = gene.cex, col = gene.col, font = 2)
	# show the Maf-S region (AGAP010405)
	mafs.region.mean <- cs['3L'] + mean(c(2785613, 2786954))
	points(mafs.region.mean, 0, pch = 19, cex = point.cex, col = point.col)
	text(mafs.region.mean, -1.7, 'Maf-S', srt = 45, adj = 1, xpd = NA, cex = gene.cex, col = gene.col, font = 2)
	# show the GST region
	gst.region.mean <- cs['3R'] + mean(c(28580000, 28605000))
	points(gst.region.mean, 0, pch = 19, cex = point.cex, col = point.col)
	text(gst.region.mean, -1.7, 'Gstu4-\nGste3 ', srt = 45, adj = 1.1, xpd = T, cex = gene.cex, col = gene.col, font = 2)
	# show the CYP6M2-Z1 region
	cyp6m2.z1.region.mean <- cs['3R'] + mean(c(6900000, 7030000))
	points(cyp6m2.z1.region.mean, 0, pch = 19, cex = point.cex, col = point.col)
	text(cyp6m2.z1.region.mean, -1.7, 'Cyp6m2-\nCyp6z1 ', srt = 45, adj = 1.1, xpd = T, cex = gene.cex, col = gene.col, font = 2)
	# show the CYP9K1 region
	cyp9k1.region.mean <- cs['X'] + mean(c(15222000, 15257000))
	points(cyp9k1.region.mean, 0, pch = 19, cex = point.cex, col = point.col)
	text(cyp9k1.region.mean, -1.7, 'Cyp9k1', srt = 45, adj = 1, xpd = T, cex = gene.cex, col = gene.col, font = 2)
	# Plot the outline of the chromosomes
	chrom.col <- if (!is.null(chrom.col)) chrom.col else c('2R' = 'black', '2L' = 'black', '3R' = 'black', '3L' = 'black', 'X' = 'black')
	chrom.y <- -8.5 - chrom.offset
	lines(c(ce['2R'], ce['2R'] - gaps/2, cs['2R'], cs['2R'], ce['2R'] - gaps/2, ce['2R']), 
		  c(-0.2, -1, -1, 1, 1, 0.2), lwd = 2, col = chrom.col['2R'])
	text((cs['2R'] + ce['2R'])/2, chrom.y, '2R', adj = 0.5, xpd = NA, cex = chrom.cex)
	lines(c(cs['2L'], cs['2L'] + gaps/2, ce['2L'], ce['2L'], cs['2L'] + gaps/2, cs['2L']), 
		  c(0.2, 1, 1, -1, -1, -0.2), lwd = 2, col = chrom.col['2L'])
	text((cs['2L'] + ce['2L'])/2, chrom.y, '2L', adj = 0.5, xpd = NA, cex = chrom.cex)
	lines(c(ce['3R'], ce['3R'] - gaps/2, cs['3R'], cs['3R'], ce['3R'] - gaps/2, ce['3R']), 
		  c(-0.2, -1, -1, 1, 1, 0.2), lwd = 2, col = chrom.col['3R'])
	text((cs['3R'] + ce['3R'])/2, chrom.y, '3R', adj = 0.5, xpd = NA, cex = chrom.cex)
	lines(c(cs['3L'], cs['3L'] + gaps/2, ce['3L'], ce['3L'], cs['3L'] + gaps/2, cs['3L']), 
		  c(0.2, 1, 1, -1, -1, -0.2), lwd = 2, col = chrom.col['3L'])
	text((cs['3L'] + ce['3L'])/2, chrom.y, '3L', adj = 0.5, xpd = NA, cex = chrom.cex)
	lines(c(cs['X'], cs['X'], ce['X'] - gaps/2, ce['X'], ce['X'], ce['X'] - gaps/2, cs['X']), 
		  c(-1, 1, 1, 0.2, -0.2, -1, -1), lwd = 2, col = chrom.col['X'])
	text((cs['X'] + ce['X'])/2, chrom.y, 'X', adj = 0.5, xpd = NA, cex = chrom.cex)
}

# Now write a function to draw the SNPs on a heatmap, marking the SNPs clumps. 
snp.heatmap.and.signif <- function (significance.value, 
                                    x, 
                                    snp.pos,
                                    clumps,
                                    filename,
                                    chrom.col = NULL,
                                    legend = T, 
                                    revC = F, 
                                    scale = c("row", "column", "none"), 
                                    na.rm = TRUE, 
                                    segment.lwd = 0.1,
                                    ...) {
	x <- copy(x)
    scale <- if (missing(scale)) 
        "none"
    else match.arg(scale)
	di <- dim(x)
	if (di[1] != di[2]) stop('input matrix should be square')
    nr <- di[1]
	# Make sure the chromosomes are in the order we want
	snp.pos$Chrom <- factor(as.character(snp.pos$Chrom), levels = c('2R', '2L', '3R', '3L', 'X'))
	snp.order <- order(snp.pos$Chrom, snp.pos$Pos)
	snp.pos <- snp.pos[snp.order, ]
	x <- x[snp.order, snp.order]
	significance.value <- significance.value[snp.order]
	# Scale if required
    if (scale == "row") {
        x <- sweep(x, 1, rowMeans(x, na.rm = na.rm), check.margin = FALSE)
        sx <- apply(x, 1, sd, na.rm = na.rm)
        x <- sweep(x, 1, sx, "/", check.margin = FALSE)
    }
    else if (scale == "column") {
        x <- sweep(x, 2, colMeans(x, na.rm = na.rm), check.margin = FALSE)
        sx <- apply(x, 2, sd, na.rm = na.rm)
        x <- sweep(x, 2, sx, "/", check.margin = FALSE)
    }
	if (!missing(filename)){
		if (legend)
			file.width = 7.2
		else
			file.width = 6.5
		file.height = 9
		if (grepl('\\.eps', filename))
			postscript(filename, width = file.width, height = file.height, horizontal = F, onefile = FALSE, paper = "special")
		else if (grepl('\\.png', filename))
			png(filename, width = file.width, height = file.height, units = 'in', res = 400)
		else if (grepl('\\.svg', filename))
			svg(filename, width = file.width, height = file.height)
		else if (grepl('\\.tif', filename))
			tiff(filename, width = file.width, height = file.height, units = 'in', res = 400)
	}
	lhei <- c(1,4,0.8)
	if (legend){
		lmat <- matrix(c(1,5,2,3,4,4), 3, 2, byrow = T)
		lwid <- c(10,1)
		layout(lmat, heights = lhei, widths = lwid)
	}
	else {
		lmat <- matrix(1:3, 3, 1)
		layout(lmat, heights = lhei)
	}
	side.margin <- 0.3
	# Plot the SNP significance
    par(mar = c(0, side.margin, 0.1, 0), mgp = c(2, 0.6, 0), xpd = NA, ...)
	chrom.col <- if (!is.null(chrom.col)) chrom.col else c('2R' = 'black', '2L' = 'black', '3R' = 'black', '3L' = 'black', 'X' = 'black')
	sig.point.cex = 0.6
	plot(-log10(significance.value), bty = 'n', xlab = '', ylab = '', xaxt = 'n', yaxt = 'n', pch = 21, cex = sig.point.cex, col = 'grey50', bg = chrom.col[snp.pos$Chrom], lwd = sig.point.cex/2, xaxs = 'i')
	axis(4, line = 0.5)
	mtext('log10(P)', 4, 2.5)
	# Plot the heatmap
	bottom.margin <- ifelse(missing(clumps), 0.1, 0.2)
    par(mar = c(bottom.margin, side.margin, 0.1, 0), xpd = NA, ...)
    if (revC)
        x <- x[, nr:1]
    image(1:nr, 1:nr, x, xlim = 0.5 + c(0, nr), ylim = 0.5 + c(0, nr), 
	      axes = FALSE, xlab = "", ylab = "", breaks = palbreaks, col = palette)
	coords.1 <- par('usr')	
	xrange.1 <- coords.1[2] - coords.1[1]
	# Add bars indicating the location of clumps
	if (!missing(clumps)){
		sapply(clumps, function(cc) {
						   clump.range <- which(snp.pos$Chrom == cc$chrom & snp.pos$Pos >= cc$start & snp.pos$Pos <= cc$end)
						   rect(
							   min(clump.range),
							   coords.1[1],
							   max(clump.range),
							   -(coords.1[1] + xrange.1/200),
							   col = 'grey30',
							   border = NA
						   )
					   }
		)
	}
	
	# Plot the legend if necessary
	if (legend){
		par(mar = c(0,0,0,0))
		plot(0:1, 0:1, type = 'n', bty = 'n', xlab = '', ylab = '', xaxt = 'n', yaxt = 'n')
		image(seq(0.1, 0.7, length.out = length(palette)), seq(0.5, 0.9, length.out = length(palette) + 1), matrix(palbreaks[2:length(palbreaks)] - diff(palbreaks)/2, length(palette), length(palette), byrow = T), col = palette, breaks = palbreaks, add = T)
		text(rep(0.85, 4), c(0.94, 0.9, 0.7, 0.5), c('r', '1', '0', '-1'), cex = c(1.4, 1, 1, 1))
	}
	
	# Plot the chromosomes
	# get a vector of end points for each chromosome (ie: cumulative sizes + gaps)
	gaps <- 5000000
	ce <- cumsum(chrom.size + c(0, 0, gaps, 0, gaps))
	# get a vector of start points for each chromosome (ie: cumulative sizes + gaps)
	cs <- ce - chrom.size
	top.margin <- 2.5
	par(mar = c(2, side.margin, top.margin, 0), xpd = NA)
	add.chroms(gaps = gaps, cs = cs, ce = ce, point.cex = 1, gene.cex = 0.9, chrom.cex = 1.2)
	coords.2 <- par('usr')
	xrange.2 <- coords.2[2] - coords.2[1]
	# Write a function to convert old x coords to new x coords
	f <- if (legend) 1.1035 else 1
	convert <- function(x) (((x - coords.1[1])/xrange.1) * xrange.2 + coords.2[1])/f
	# And work out the top end of the current plot
	top.end <- par('cxy')[2] * par('cex') * par('lheight') * top.margin + coords.2[4]
	# Add arrows linking the genome picture to the correct positions on the heatmap. 
	# Get a column that shows the position along the whole genome ,rather than just on the chromosome
	if (!('Genome.pos' %in% names(snp.pos)))
		snp.pos$Genome.pos <- snp.pos$Pos + cs[as.character(snp.pos$Chrom)]
	# For every snp, draw a line from the matrix to the chromosome.
	segments(snp.pos$Genome.pos, 1.1, convert(1:nrow(x)), top.end, lwd = segment.lwd, col = chrom.col[snp.pos$Chrom])
	if (!missing(filename))
		dev.off()
}

# As we to along, we will record all of the clumps, whether in 2La / Ace1 or not, and output them
# to a csv
regions.list <- list()
window.size <- 100000
for (ss in sample.sets){
	R.env <- paste('../', ss, 'classical_analysis_sibgroups/classical_analysis_nosibs.Rdata', sep = '/')
	load(R.env)
	# Before saving the .Rdata file, the top.snps table was re-ordered such that it is no longer in 
	# the same order as top.snp.residual.cor.matrix, which causes problems for the snp.heatmap.and.signif
	# function. So we put it back in the original order here.
	top.snps <- top.snps[order(Chrom, Pos)]
	
	# "snp.clumps" contains all clumps of any size (we want a minimum size of 10), while the object
	# "clumps" is the manually filtered clumps of sufficient size not in Ace1 or 2La where 2La 
	# dominates. 
	# First collate all regions with a clump of at least 10 SNPs
	region.ids <- names(snp.clumps)[snp.clumps >= 10]
	regions.list[[ss]] <- tstrsplit(region.ids, ':') %>%
	                       as.data.table() %>%
	                       .[, pos := as.numeric(V2) * window.size] %>%
	                       .[, .(sample.set = ss,
	                             chrom = V1, 
	                             start = pos - window.size/2, 
	                             end = pos + window.size/2
	                            ) 
	                        ]
	# Then make the summary plots of all the top SNPs, regardless of whether they are in a clump
	if (exists('clumps')){
		snp.heatmap.and.signif(top.snps$logregP, 
							   top.snp.residuals.cor.matrix, 
							   top.snps[, .(Chrom, Pos)],
							   clumps,
							   chrom.col = chrom.colours,
							   filename = paste(ss, 'top_snps_summary_plot.png', sep = '_'),
							   family = 'Arial',
							   segment.lwd = 0.4)
	}
	else {
		snp.heatmap.and.signif(top.snps$logregP, 
							   top.snp.residuals.cor.matrix, 
							   top.snps[, .(Chrom, Pos)],
							   chrom.col = chrom.colours,
							   filename = paste(ss, 'top_snps_summary_plot.png', sep = '_'),
							   family = 'Arial',
							   segment.lwd = 0.4)
	}
	rm(list = c('top.snps', 'top.snp.residuals.cor.matrix', 'clumps'))
}

regions <- rbindlist(regions.list)
fwrite(regions, 'classical_analysis_snp_clump_regions.csv', sep = '\t')



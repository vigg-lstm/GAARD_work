# A function to lighten the tone of the colour
lighten.col <- function(color, lightness, alpha = alpha){
	col.rgb <- col2rgb(color)/255
	rgb(t(1-(1-col.rgb)*lightness), alpha = alpha)
}

# Write a function to draw the chromosomes on an existing plotting device
add.chromosomes <- function(gaps, chrom.sizes, gene.cex = 0.9, gene.col = 'grey20', point.cex = 1.2, point.col = 'grey30', chrom.col = NULL, chrom.cex = 1.4, chrom.offset = 0){
	ce <- cumsum(chrom.sizes + c(0, 0, gaps, 0, gaps))
	cs <- ce - chrom.sizes
	plot(c(cs[1], ce[5]), c(-6.5,1.3), xlim = c(cs[1] + gaps/2, ce[5] - gaps/2), type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
	# show the kdr region
	kdr.region.mean <- cs['2L'] + mean(c(2358158, 2431617))
	points(kdr.region.mean, 0, pch = 19, cex = point.cex, col = point.col)
	text(kdr.region.mean, -1.7, 'Vgsc', srt = 45, adj = 1.1, xpd = NA, cex = gene.cex, col = gene.col, font = 2)
	# show the Rdl region
	rdl.region.mean <- cs['2L'] + mean(c(25363652, 25434556))
	points(rdl.region.mean, 0, pch = 19, cex = point.cex, col = point.col)
	text(rdl.region.mean, -1.7, 'Rdl', srt = 45, adj = 1.1, xpd = NA, cex = gene.cex, col = gene.col, font = 2)
	# show the Carboxylesterase region
	coeae.region.mean <- cs['2L'] + mean(c(28548433, 28550748))
	points(coeae.region.mean, 0, pch = 19, cex = point.cex, col = point.col)
	text(coeae.region.mean, -1.7, 'Coeae2f', srt = 45, adj = 1.1, xpd = NA, cex = gene.cex, col = gene.col, font = 2)
	# show the Ace1 gene region
	ace1.region.mean <- cs['2R'] + mean(c(3484107, 3495790))
	points(ace1.region.mean, 0, pch = 19, cex = point.cex, col = point.col)
	text(ace1.region.mean, -1.7, 'Ace1', srt = 45, adj = 1.1, xpd = NA, cex = gene.cex, col = gene.col, font = 2)
	# show the CYP6 region
	cyp6.region.mean <- cs['2R'] + mean(c(28463000, 28568000))
	points(cyp6.region.mean, 0, pch = 19, cex = point.cex, col = point.col)
	text(cyp6.region.mean, -1.7, 'Cyp6aa1-\nCyp6p2 ', srt = 45, adj = 1.1, xpd = NA, cex = gene.cex, col = gene.col, font = 2)
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
	text(gst.region.mean, -1.7, 'Gstu4-\nGste3 ', srt = 45, adj = 1.1, xpd = NA, cex = gene.cex, col = gene.col, font = 2)
	# show the CYP6M2-Z1 region
	cyp6m2.z1.region.mean <- cs['3R'] + mean(c(6900000, 7030000))
	points(cyp6m2.z1.region.mean, 0, pch = 19, cex = point.cex, col = point.col)
	text(cyp6m2.z1.region.mean, -1.7, 'Cyp6m2-\nCyp6z1 ', srt = 45, adj = 1.1, xpd = NA, cex = gene.cex, col = gene.col, font = 2)
	# show the CYP9K1 region
	cyp9k1.region.mean <- cs['X'] + mean(c(15222000, 15257000))
	points(cyp9k1.region.mean, 0, pch = 19, cex = point.cex, col = point.col)
	text(cyp9k1.region.mean, -1.7, 'Cyp9k1', srt = 45, adj = 1, xpd = NA, cex = gene.cex, col = gene.col, font = 2)
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


draw.gene.model <- function(start.pos, end.pos, gene.positions, exon.positions, include.gene.names = T, y = 0, text.cex = 0.5, gene.thickness.fraction = 15, lwd = 2){
	lines(c(start.pos, end.pos), c(y, y), lwd = lwd, col = 'grey20')
	if (nrow(gene.positions) > 0){
		gene.positions <- gene.positions[order(start)]
		gr <- par('usr')
		gene.thickness <- (gr[4] - gr[3])/gene.thickness.fraction
		# Draw the genes
		v.adj = ifelse(gene.positions$strand == '+', 0, -gene.thickness)
		rect(
			gene.positions[,start], 
			y + v.adj, 
			gene.positions[,end], 
			y + gene.thickness + v.adj, 
			col = 'grey95',
			border = 'black',
			lwd = lwd,
			xpd = F
		)
		# Draw the exons
		exon.v.adj = ifelse(exon.positions$strand == '+', 0, -gene.thickness)
		rect(
			exon.positions[,start], 
			y + exon.v.adj, 
			exon.positions[,end], 
			y + gene.thickness + exon.v.adj, 
			col = 'black',
			border = NA,
			xpd = F
		)
		# Label the genes
		if (include.gene.names){
			text(
				apply(gene.positions[,.(start, end)], 1, mean), 
				y - 1.6*gene.thickness, 
				gene.positions$gene.name, 
				adj = 1, 
				srt = 45,
				cex = text.cex
			)
		}
	}
}


library(data.table)
library(lme4)
library(stringr)
library(ape)

# Function to perform logistic regression for genotype-phenotype association.
logreg.test.random.effect <- function(genotype, phenotype, location){
	genotype <- as.numeric(genotype)
	null.model <- glm(phenotype ~ 1 + location, family = 'binomial')
	model1 <- glm(phenotype ~ genotype + location, family = 'binomial')
	pval <- anova(model1, null.model, test = 'Chisq')[['Pr(>Chi)']][2]
	direction <- levels(phenotype)[sign(model1$coefficients['genotype'])/2+1.5]
	data.table(logregP = pval, direction = direction)
}

haplotypes <- fread('Madina_Obuasi_delta_shared_Cyp6_cluster.csv', sep = '\t', key = 'sample_name')
sig.test <- haplotypes[, cbind(cluster.size = sum(Cyp6_main_cluster), logreg.test.random.effect(Cyp6_main_cluster, as.factor(phenotype), location))]

cat('\nSignificance test of shared haplotype combining Madina and Obuasi:\n')
print(sig.test)

# A function that will perform a fisher test on cluster allele frequencies and report whether the mutant at
# each allele is more frequent in the haplotype that outside of it, with a p.value that is smaller than a 
# threshold. This is effectively just checking that a SNP present in a haplotype is actually associated with
# that haplotype compared to the wt. 
test.fisher.p <- function(cluster_ref, cluster_alt, non_cluster_ref, non_cluster_alt, threshold = 0.01){
	freq.test <- (cluster_alt / cluster_ref) > (non_cluster_alt / non_cluster_ref)
	d <- c(cluster_ref, cluster_alt, non_cluster_ref, non_cluster_alt)
	p.test <- fisher.test(matrix(d, 2, 2))$p.value <= threshold
	freq.test & p.test
}

# Load the SNP data
csv.name = 'Madina_Obuasi_delta_shared_Cyp6_cluster_hap_sequence.csv'
vcf.name = 'Madina_Obuasi_delta_shared_Cyp6_cluster.vcf'
snp.table <- merge(
	fread(csv.name), 
	fread(vcf.name)[, -c('ID', 'QUAL', 'FILTER')], 
	by.x = c('Chrom', 'Pos', 'Ref', 'Alt'), 
	by.y = c('#CHROM', 'POS', 'REF', 'ALT')
)
# Note which snps are non-synonymous
snp.table[, non.synonymous := grepl('missense', INFO)]
# Note which snps have significantly different frequencies between clusters and non-clusters. 
snp.table[, isdiff := mapply(test.fisher.p, n_Ref, n_Alt, n_Ref_non_cluster, n_Alt_non_cluster)]

# Write a function to plot each cluster. 
# We want a bead diagram of the genes in the region, with coloured points along it showing SNPs. Also, cluster
# size should be indicated, as should p-value
gff.path <- '~/data/ML/GAARD_SNP/data/Anopheles-gambiae-PEST_BASEFEATURES_AgamP4.9.gff3'
gff <- as.data.table(read.gff(gff.path, GFF3 = T))

draw.gene.model <- function(start.pos, end.pos, gene.positions, exon.positions, include.gene.names = T, y = 0, text.cex = 0.5, gene.thickness.fraction = 15){
	lwd = 2
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

# A function to lighten the tone of the colour
lighten.col <- function(color, lightness, alpha = alpha){
	col.rgb <- col2rgb(color)/255
	rgb(t(1-(1-col.rgb)*lightness), alpha = alpha)
}

# A function to plot the haplotype clusters. The remove.background.snps argument controls whether SNPs
# where the alt alleles is not significantly more frequent in the cluster than in non-cluster haplotypes
# are removed from the plot or simply plotted in a lighter colour. 
plot.hap.cluster <- function(snp.table, sig.test, filename, point.cex, remove.background.snps = T, add = F){
	chrom <- unique(snp.table$Chrom)
	start.pos <- min(snp.table$Pos)
	end.pos <- max(snp.table$Pos)
	
	# Get the list of genes present in the window
	gene.gff <- gff[seqid == chrom & 
	                start < end.pos & 
	                end > start.pos &
	                type == 'gene', ]
	
	gene.gff[, gene.id := unique(str_extract(attributes, 'AGAP\\d{6}'))]
	gene.gff[, gene.name := str_extract(attributes, '(?<=Name=)[^;]+') %>%
	                        str_to_title]
	gene.gff$gene.name[is.na(gene.gff$gene.name)] <- gene.gff$gene.id[is.na(gene.gff$gene.name)]
	# Now the list of Exons
	exon.gff <- gff[seqid == chrom & 
	                start < end.pos & 
	                end > start.pos &
	                type == 'exon', ]
	
	plot.height = 1.5
	plot.width = 6.5
	if (!missing(filename)){
		if (grepl('\\.eps', filename))
			postscript(filename, width = plot.width, height = plot.height, horizontal = F, onefile = FALSE, paper = "special")
		else if (grepl('\\.png', filename))
			png(filename, width = plot.width, height = plot.height, units = 'in', res = 200)
		else if (grepl('\\.svg', filename))
			svg(filename, width = plot.width, height = plot.height)
		else if (grepl('\\.tif', filename))
			tiff(filename, width = plot.width, height = plot.height, units = 'in', res = 200)
	}
	else if (!add){
		x11(width = plot.width, height = plot.height)
	}
	
	layout(
		matrix(1:2, 2, 1), 
		heights = c(0.4, 0.1)
	)
	par(
		mar = c(0,0,0,0), 
		omi = c(0.4,0.6,0.175,0.175), 
		cex = 0.75,
		mex = 1.3,
		mgp = c(1,0.3,0),
		tcl = -0.2,
		xaxs = 'i',
		xpd = NA
	)
	plot(c(start.pos, end.pos), c(0, 0.3 + 0.2), type = 'n', bty = 'n', xaxt = 'n', xlab = '', yaxt = 'n', ylab = '')
	alt.column <- 'n_Alt'
	ref.column <- 'n_Ref'
	cluster.size <- sig.test[, cluster.size]
	cluster.p <- sig.test[, logregP]
	cluster.dir <- sig.test[, direction]
	y.pos <- 0.3 + 0.2
	draw.gene.model(
		start.pos, 
		end.pos, 
		gene.gff, 
		exon.gff, 
		y.pos, 
		include.gene.names = T,
		text.cex = 0.5,
		gene.thickness.fraction = 15
	)
	# Set the minimum SNP point size
	hpc <- par('cex')/2
	# Pull out the SNPs that have at least some alt alleles in the cluster
	snp.sub.table <- snp.table[n_Alt > 0]
	# Set up the colours and pch
	snp.sub.table[, 
		`:=`(bg = 'deepskyblue', col = 'deepskyblue3', pch = 21)]
	snp.sub.table[(!non.synonymous & grepl('synonymous_variant', INFO)), 
		`:=`(bg = 'yellow', col = 'goldenrod3', pch = 21)]
	snp.sub.table[(non.synonymous), 
		`:=`(bg = 'violetred1', col = 'violetred4', pch = 23)]
	if (remove.background.snps)
		snp.sub.table <- snp.sub.table[isdiff == T, ]
	else 
		# Lighten the colour of SNPs that are not different from the genetic background of the population
		snp.sub.table[isdiff == F, `:=`(bg = lighten.col(bg, 0.2), col = lighten.col(col, 0.2))]
	# Now add the points to the plot
	snp.sub.table[, points(
		Pos, 
		rep(y.pos, .N), 
		bg = bg,
		col = col,
		cex = hpc + hpc * n_Alt / ..cluster.size, 
		pch = pch
	)]
	# Now indicate cluster size and P-value
	label.x <- start.pos - (end.pos - start.pos)/50
	arrow.code <- ifelse(cluster.dir == 'alive', '\u2191', '\u2193')
	text(label.x, y.pos + 0.04, paste('n =', cluster.size), adj = 1, cex = 0.6)
	p.col <- ifelse(cluster.p < 0.05, 'red', 'black')
	text(label.x, y.pos - 0.04, paste('p =', format(cluster.p, digits = 2, scientific = -2)), adj = 1, col = p.col, cex = 0.6)
	text(label.x, y.pos - 0.04, arrow.code, adj = 0, cex = 0.8)
	
	plot(c(start.pos, end.pos), c(0,0), type = 'n', bty = 'n', xlab = paste('Position on', chrom), yaxt = 'n', ylab = '', cex.lab = 0.8, cex.axis = 0.6)
	if (!missing(filename))
		dev.off()
}

plot.hap.cluster(snp.table, sig.test, 'Madina_Obuasi_shared_Cyp6_haplotype.svg', point.cex = 0.8)



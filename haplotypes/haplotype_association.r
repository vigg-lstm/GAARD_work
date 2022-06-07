library(data.table)
library(lme4)
library(stringr)
library(ape)

study.pops <- c('Avrankou_coluzzii_Delta',
                'Baguida_gambiae_Delta',
                'Baguida_gambiae_PM',
                'Korle-Bu_coluzzii_Delta',
                'Korle-Bu_coluzzii_PM',
                'Madina_gambiae_Delta',
                'Madina_gambiae_PM',
                'Obuasi_gambiae_Delta',
                'Obuasi_gambiae_PM')

python.env <- '~/miniconda3/envs/gaard/bin/python3'

# Function to perform logistic regression for genotype-phenotype association.
logreg.test <- function(genotype, phenotype){
	genotype <- as.numeric(genotype)
	null.model <- glm(phenotype ~ 1, family = 'binomial')
	model1 <- glm(phenotype ~ genotype, family = 'binomial')
	pval <- anova(model1, null.model, test = 'Chisq')[['Pr(>Chi)']][2]
	direction <- levels(phenotype)[sign(model1$coefficients['genotype'])/2+1.5]
	data.table(logregP = pval, direction = direction)
}

# Get the window ranges 
window.ranges <- lapply(
	setNames(nm = study.pops), 
	function(pop) 
		fread(paste('../focal_gwas', pop, 'focal_window_ranges.csv', sep = '/'))
)

sig.test.list <- list()
for (pop in study.pops){
	sig.tests.this.pop <- list()
	# Loop through the windows. This code means nrow = 0 results in not entering the loop
	for (i in seq(length.out = nrow(window.ranges[[pop]]))){
		window.range <- paste(
			window.ranges[[pop]][i, chrom],
			paste(window.ranges[[pop]][i, .(start, end)], collapse = '-'),
			sep = ':'
		)
		# Assign haplotypes using python script
		assign.hap.command <- paste(python.env, 'assign_haplotypes.py', pop, window.range, '..')
		system(assign.hap.command)
		# Load the results
		haplotypes <- fread(paste(pop, '_', window.range, '.csv', sep = ''), sep = '\t', key = 'sample_name')
		# Perform association analysis
		sig.test <- grep('cluster_', names(haplotypes), value = T) %>%
		            setNames(nm = .) %>%
		            lapply(function(cluster) haplotypes[, cbind(cluster.size = sum(get(cluster)), logreg.test(get(cluster), as.factor(phenotype)))]) %>%
		            rbindlist(idcol = 'cluster.name')
		sig.tests.this.pop[[window.range]] <- sig.test
	}
	sig.test.list[[pop]] <- rbindlist(sig.tests.this.pop, idcol = 'window')
}
sig.tests <- rbindlist(sig.test.list, idcol = 'population')

fwrite(sig.tests, file = 'haplotype_significance_tests.csv', sep = '\t')

# A function that will perform a fisher test on cluster allele frequencies and report whether the mutant is 
# more frequency in the haplotype, with a p.value that is smaller than a threshold
test.fisher.p <- function(cluster_ref, cluster_alt, non_cluster_ref, non_cluster_alt, threshold = 0.01){
	freq.test <- (cluster_alt / cluster_ref) > (non_cluster_alt / non_cluster_ref)
	d <- c(cluster_ref, cluster_alt, non_cluster_ref, non_cluster_alt)
	p.test <- fisher.test(matrix(d, 2, 2))$p.value <= threshold
	freq.test & p.test
}

# Write a function to load all the haplotype data for each window. This will then create a list which we can 
# save as an R object that we can then put in the .gitignore
load.snps <- function(pop.window){
	csv.name <- paste(pop.window, '_hap_sequence.csv', sep = '')
	vcf.name <- paste(pop.window, '.vcf', sep = '')
	snp.eff <- merge(
		fread(csv.name), 
		fread(vcf.name)[, -c('ID', 'QUAL', 'FILTER')], 
		by.x = c('Chrom', 'Pos', 'Ref', 'Alt'), 
		by.y = c('#CHROM', 'POS', 'REF', 'ALT')
	)
	# Note which snps are non-synonymous
	snp.eff[, non.synonymous := grepl('missense', INFO)]
	# Note which snps are different between clusters and non-clusters. 
	clusters <- str_extract(names(snp.eff), 'cluster_\\d+$') %>%
	            unlist %>%
	            unique %>%
	            {.[!is.na(.)]}
	for (cluster in clusters){
		snp.eff[, paste(cluster, '_isdiff', sep = '') := mapply(test.fisher.p, 
		                                                        get(paste('n_Ref_', cluster, sep = '')),
		                                                        get(paste('n_Alt_', cluster, sep = '')),
		                                                        n_Ref_non_cluster,
		                                                        n_Alt_non_cluster)
		       ]
	}
	snp.eff
}

# Now load all the haplotype data
sig.tests$pop.window <- paste(sig.tests$population, sig.tests$window, sep = '_')
setkey(sig.tests, pop.window, cluster.name)
snp.tables <- lapply(setNames(nm = sig.tests$pop.window), load.snps)

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
plot.hap.cluster <- function(pop.win, filename, point.cex, remove.background.snps = T, add = F){
	chrom <- str_extract(pop.win, '(?<=_)[23]?[LRX](?=:)')
	start.pos <- as.integer(str_extract(pop.win, '(?<=:)\\d+(?=-)'))
	end.pos <- as.integer(str_extract(pop.win, '(?<=-)\\d+$'))
	# Get the SNPs
	snp.table <- snp.tables[[pop.win]]
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
	clusters <- str_extract(names(snp.table), 'cluster_\\d+$') %>%
	            unlist %>%
	            unique %>%
	            {.[!is.na(.)]}
	
	plot.height = 0.2 * (length(clusters) - 1) + 1.5
	plot.width = 6.5
	if (!missing(filename)){
		if (grepl('\\.eps', filename))
			postscript(filename, width = plot.width, height = plot.height, horizontal = F, onefile = FALSE, paper = "special")
		else if (grepl('\\.png', filename))
			png(filename, width = plot.width, height = plot.height, units = 'in', res = 200)
		else if (grepl('\\.tif', filename))
			tiff(filename, width = plot.width, height = plot.height, units = 'in', res = 200)
	}
	else if (!add){
		x11(width = plot.width, height = plot.height)
	}
	layout(
		matrix(1:2, 2, 1), 
		heights = c(0.4 + 0.2 * (length(clusters) - 1), 0.1)
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
	plot(c(start.pos, end.pos), c(0, 0.3 + length(clusters) * 0.2), type = 'n', bty = 'n', xaxt = 'n', xlab = '', yaxt = 'n', ylab = '')
	for (i in 1:length(clusters)){
		cluster <- clusters[i]
		alt.column <- paste('n_Alt', cluster, sep = '_')
		ref.column <- paste('n_Ref', cluster, sep = '_')
		cluster.size <- sig.tests[.(pop.win, cluster), cluster.size]
		cluster.p <- sig.tests[.(pop.win, cluster), logregP]
		cluster.dir <- sig.tests[.(pop.win, cluster), direction]
		y.pos <- 0.3 + (length(clusters) - i) * 0.2
		draw.gene.model(
			start.pos, 
			end.pos, 
			gene.gff, 
			exon.gff, 
			y.pos, 
			include.gene.names = (i == length(clusters)),
			text.cex = 0.5,
			gene.thickness.fraction = 15
		)
		# Set the minimum SNP point size
		hpc <- par('cex')/2
		# Pull out the SNPs that have at least some alt alleles in the cluster
		snp.sub.table <- snp.table[get(alt.column) > 0]
		# Set up the colours and pch
		snp.sub.table[, 
			`:=`(bg = 'deepskyblue', col = 'deepskyblue3', pch = 21)]
		snp.sub.table[(!non.synonymous & grepl('synonymous_variant', INFO)), 
			`:=`(bg = 'yellow', col = 'goldenrod3', pch = 21)]
		snp.sub.table[(non.synonymous), 
			`:=`(bg = 'violetred1', col = 'violetred4', pch = 23)]
		if (remove.background.snps)
			snp.sub.table <- snp.sub.table[get(paste(cluster, '_isdiff', sep = '')), ]
		else 
			# Lighten the colour of SNPs that are not different from the genetic background of the population
			snp.sub.table[get(paste(cluster, '_isdiff', sep = '')) == F, 
				`:=`(bg = lighten.col(bg, 0.2), col = lighten.col(col, 0.2))]
		# Now add the points to the plot
		snp.sub.table[, points(
			Pos, 
		    rep(y.pos, .N), 
		    bg = bg,
		    col = col,
		    cex = hpc + hpc * get(alt.column) / ..cluster.size, 
		    pch = pch
		)]
		# Now indicate cluster size and P-value
		label.x <- start.pos - (end.pos - start.pos)/50
		arrow.code <- ifelse(cluster.dir == 'alive', '\u2191', '\u2193')
		text(label.x, y.pos + 0.04, paste('n =', cluster.size), adj = 1, cex = 0.6)
		p.col <- ifelse(cluster.p < 0.05, 'red', 'black')
		text(label.x, y.pos - 0.04, paste('p =', format(cluster.p, digits = 2, scientific = -2)), adj = 1, col = p.col, cex = 0.6)
		text(label.x, y.pos - 0.04, arrow.code, adj = 0, cex = 0.8)
	}
	plot(c(start.pos, end.pos), c(0,0), type = 'n', bty = 'n', xlab = paste('Position on', chrom), yaxt = 'n', ylab = '', cex.lab = 0.8, cex.axis = 0.6)
	mtext(pop.win, outer = T, cex = 0.8)
	if (!missing(filename))
		dev.off()
}

for(pop.win in names(snp.tables))
	plot.hap.cluster(pop.win, paste(pop.win, '.png', sep = ''), point.cex = 0.8)



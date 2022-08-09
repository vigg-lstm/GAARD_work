library(data.table)
#library(plotrix)

# The aim of this script is to correlate the presence of the bellwether SNPs and the significant Asaia reads
# to see whether they are both driven by the same signal. 

# Function to lighten a colour
lighten.col <- function(colour, lightness){
	col.rgb <- col2rgb(colour)/255
	red <- 1-(1-col.rgb[1, ])*lightness
	blue <- 1-(1-col.rgb[2, ])*lightness
	green <- 1-(1-col.rgb[3, ])*lightness
	rgb(cbind(red, blue, green))
}

# Function to add transparency to a vector of colours
transp <- function(col, alpha = 0.5) {
	rgb.value <- col2rgb(col)/255
	t.col <- rgb(t(rgb.value), alpha = alpha)
	invisible(t.col)
}

# First, load the top bellwether SNPs for Korle-Bu Delta
kb.delta.bw.filename <- '../bellwhether_SNPs/Korle-Bu_coluzzii_Delta/classical_analysis/filtered_top_snp_table.csv'
kb.delta.bw <- fread(kb.delta.bw.filename)
sample.names <- colnames(kb.delta.bw)[grepl('WA', colnames(kb.delta.bw))]
num.snps <- nrow(kb.delta.bw)
par(mfrow = rep(num.snps, 2), cex = 0.5, mar = c(1,1,1,1), mgp = c(2, 0.3, 0))
for (i in 1:num.snps){
	for (j in 1:num.snps){
		if (i == j)
			plot(0,0, type = 'n', bty = 'n', xlab = '', ylab = '', xaxt = 'n', yaxt = 'n')
		else{
			genotype.counts <- table(as.numeric(kb.delta.bw[i, ..sample.names]), as.numeric(kb.delta.bw[j, ..sample.names]))
			x.coords <- as.numeric(colnames(genotype.counts)) + 1
			y.coords <- as.numeric(rownames(genotype.counts)) + 1
			plot(rep(x.coords, each = length(y.coords)), rep(y.coords, length(x.coords)), cex = log10(genotype.counts), pch = 19, ylim = c(0.8, max(y.coords) + 0.2), xlim = c(0.8, max(x.coords) + 0.2), xaxt = 'n', yaxt = 'n')
			axis(1, x.coords, labels = c('wt', 'het', 'mut')[x.coords]) 
			axis(2, y.coords, labels = c('wt', 'het', 'mut')[y.coords]) 
		}
	}
}

# Next, load the Asaia reads data. 
kb.delta.asaia.filename <- '~/Liverpool/ML/GAARD/Korle-Bu_coluzzii_Delta/acetobacter_checking/acetobacter_matching_kmers.csv'
kb.delta.asaia <- fread(kb.delta.asaia.filename)
shared.sample.names <- intersect(sample.names, colnames(kb.delta.asaia))
mean.asaia.reads <- as.numeric(kb.delta.asaia[nrow(kb.delta.asaia), ..shared.sample.names])

x11()
par(mfrow = rep(ceiling(sqrt(num.snps)), 2), cex = 0.6, mar = c(3,3,1,0), mgp = c(1.5, 0.5, 0))
for (i in 1:num.snps){
	jitter.y = 0
	jitter.x = runif(length(shared.sample.names), -0.2, 0.2)
	plot(as.numeric(kb.delta.bw[1, ..shared.sample.names]) + jitter.x, mean.asaia.reads + jitter.y, cex = 0.5)
}

# Need to control for the fact that both are correlated with phenotype. So need to calculate the residuals of 
# each of the genotypes within each phenotype, and check that these are still correlated. 

# Load the phenotypes
phenotypes.filename <- '/home/eric/Desktop/Liverpool/ML/GAARD_SNP/pre_QC/metadata/sample_phenotypes.csv'
phenotypes.table <- fread(phenotypes.filename, key = 'specimen')
phenotypes <- setNames(phenotypes.table[shared.sample.names, phenotype], shared.sample.names)

# Get the residuals for the measures of interest after accounting for phenotype
combined.table <- data.table(sample.name = shared.sample.names, asaia = mean.asaia.reads, t(kb.delta.bw[, ..shared.sample.names]), phenotype = phenotypes)
var.names <- colnames(combined.table)[c(-1, -ncol(combined.table))]
residuals.table <- combined.table[, c(.(sample.name = sample.name), sweep(.SD, 2, colMeans(.SD))), by = phenotype, .SDcols = var.names]

# Now plot those residuals against each-other
x11()
par(mfrow = rep(ceiling(sqrt(num.snps)), 2), cex = 0.6, mar = c(3,3,1,0), mgp = c(1.5, 0.5, 0))
for (i in 1:num.snps){
#	jitter.y = runif(length(shared.sample.names), -0.2, 0.2)
#	jitter.x = runif(length(shared.sample.names), -0.2, 0.2)
	jitter.x = 0
	jitter.y = 0
	which.var <- var.names[i+1]
	plot(residuals.table[[which.var]] + jitter.x, residuals.table[, asaia] + jitter.y, cex = 0.5, col = c('red', 'blue')[(phenotypes == 'alive') + 1])
}

# Those plots aren't so clear. Let's instead just use the non-residual values but make the plots separately for dead 
# and alive. 

splitplot <- function(snp.name){
	asaia.phen.split <- split(combined.table$asaia, interaction(combined.table[[snp.name]], phenotypes))
	group.order <- c("0.dead", "1.dead", "0.alive", "1.alive")
	group.pos <- setNames(c(1,2,4,5), group.order)
	boxcols <- setNames(c('orange3', 'red3','blue', 'purple4'), group.order)
	boxplot(asaia.phen.split[group.order], border = boxcols, lwd = 2, boxwex = 0.5, at = group.pos, names = rep('', length(group.order)), cex = 0, ylab = 'Mean Asaia reads', cex.lab = 1.2)
	for (group in group.order){
		this.col <- transp(boxcols[group], 0.3)
		points(group.pos[group] + runif(length(asaia.phen.split[[group]]), -0.2, 0.2), asaia.phen.split[[group]], pch = 21, col = this.col, bg = this.col, cex = 1)
	}
	axis(1, group.pos, labels = rep(c('wt', 'het'), 2), cex.axis = 1.2)
	axis(1, group.pos[c(1,3)]+0.5, labels = c('Dead', 'Alive'), tick = F, padj = 2, cex.axis = 1.2)
}

x11()
par(mfrow = rep(ceiling(sqrt(num.snps)), 2), cex = 0.6, mar = c(3,3,1,0), mgp = c(1.5, 0.5, 0))
for (snp in var.names[-1]){
#	jitter.y = runif(length(shared.sample.names), -0.2, 0.2)
#	jitter.x = runif(length(shared.sample.names), -0.2, 0.2)
	splitplot(snp)
}

# For saving, let's just pick one SNP and plot that one
png('Asaia_vs_bellwether.png', width = 4, height = 4, units = 'in', res = 300)
par(cex = 0.7)
splitplot('V3')
dev.off()
# So that's pretty conclusive. 

# Now we want to look at Asaia reads as a function of plate position. For this, we need to know what well each 
# sample occupied on its original plate and on the sanger plate, although since these bacterial reads are present
# in our GAARD samples, I don't think it's Sanger contamination. 

# Load the metadata
meta.filename <- '/home/eric/Liverpool/ML/GAARD_SNP/pre_QC/metadata/sample_metadata.csv'
meta <- fread(meta.filename, key = 'MalGEN_ID')

# There are 3 well layouts we need to deal with. First is the layout on our GAARD plates. Then there is the layout
# on the plates that we sent to Sanger. Finally, there is the layout that Sanger used when they combined plates. 
# The first two layouts can be found in the following files:
kmer.sample.names <- names(kb.delta.asaia)[grepl('WA', colnames(kb.delta.asaia))]
lstm.layouts.filename <- '/home/eric/Liverpool/ML/GAARD/data/GAARD_WestAfrica_Sanger_submission.csv'
lstm.layouts <- fread(lstm.layouts.filename, key = 'SampleID')
layout1 <- lstm.layouts[kmer.sample.names, .(SampleID, GAARDPlate, GAARDWell)]
layout2 <- lstm.layouts[kmer.sample.names, .(SampleID, WGSPlate, WGSWell)]
# The third layout can be found in the manifest that Chris Clarkson sent me
sanger.layout.filename <- '5563stdy_manifest_12930_210119.csv'
sanger.layout <- fread(sanger.layout.filename)
colnames(sanger.layout) <- gsub(' ', '_', colnames(sanger.layout))
sanger.layout$SampleID <- meta[sanger.layout[['SUPPLIER_SAMPLE_NAME']], External_ID]
setkey(sanger.layout, SampleID)
layout3 <- sanger.layout[kmer.sample.names, .(SampleID, SANGER_PLATE_ID, WELL)]
colnames(layout3) <- c('SampleID', 'SangerPlate', 'SangerWell')
layouts <- merge(merge(layout1, layout2), layout3)
layouts$asaia <- as.numeric(kb.delta.asaia[nrow(kb.delta.asaia), ..kmer.sample.names])
# Normalise the asaia values so the max is 1
layouts$asaia.norm <- layouts$asaia / max(layouts$asaia)
# Add the phenotypes
layouts$phenotype <- phenotypes.table[layouts$SampleID, phenotype]
# Outside the scope of this script, I should try the generic Asaia primers on the GAARD samples with qPCR, to see
# if Asaia in general follow the same pattern as that specific sequence. 

# Write a function to plot values of a sample based on well position
well.plot <- function(wells, values, phenotypes, transform = 5, ...){
	well.x <- setNames(rep(1:12, each = 8), outer(LETTERS[1:8], 1:12, paste, sep = ''))
	well.y <- setNames(rep(8:1, 12), outer(LETTERS[1:8], 1:12, paste, sep = ''))
	phen.col <- c('red3', 'deepskyblue4')[(phenotypes == 'alive') + 1]
	if(transform)
		values <- values^(1/transform)
	plot(well.x, well.y, type = 'n', xlab = '', ylab = '', xaxt = 'n', yaxt = 'n')
	title(..., line = 3)
	points(well.x[wells], well.y[wells], cex = 2*(0.5 + values), pch = 21, bg = lighten.col(phen.col, (values/2)+0.5), col = lighten.col(phen.col, values))
	missing.wells <- setdiff(names(well.x), wells)
	points(well.x[missing.wells], well.y[missing.wells], pch = 4, cex = 2, col = 'grey50', lwd = 2)
	axis(2, at = 8:1, labels = LETTERS[1:8], tick = F, las = 1)
	axis(3, at = 1:12, labels = 1:12, tick = F)
}

# Now a function to plot this for all the plates in a given layout
plate.well.plot <- function(plates, wells, values, phenotypes){
	w <- split(wells, plates)
	v <- split(values, plates)
	p <- split(phenotypes, plates)
	num.plates <- length(w)
	num.plot.cols <- ceiling(sqrt(num.plates)) + 1
	num.plot.rows <- (num.plates %% num.plot.cols) + 1
	par(mfrow = c(num.plot.rows, num.plot.cols))
	for (pl in names(w))
		well.plot(w[[pl]], v[[pl]], p[[pl]], main = pl)
}
png('plate_layouts.png', width = 12, height = 6, units = 'in', res = 300)
with(layouts, plate.well.plot(GAARDPlate, GAARDWell, asaia.norm, phenotype))
dev.off()

x11()
with(layouts, plate.well.plot(WGSPlate, WGSWell, asaia.norm, phenotype))
x11()
with(layouts, plate.well.plot(SangerPlate, SangerWell, asaia.norm, phenotype))


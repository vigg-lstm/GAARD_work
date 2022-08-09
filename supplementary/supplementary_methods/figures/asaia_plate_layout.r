library(data.table)
library(magrittr)
#library(plotrix)

# The aim of this script is to correlate the presence of the bellwether SNPs and the significant Asaia reads
# to see whether they are both driven by the same signal. 

# Function to lighten a colour and add transparency
lighten.col <- function(colour, lightness, alpha = 1){
	col.rgb <- col2rgb(colour)/255
	new.rgb <- 1-(1-col.rgb)*lightness
	t.col <- rgb(t(new.rgb), alpha = alpha)
	t.col
}

# A function to load the bracken data for asaia
load.asaia <- function(sample.name, bracken.folder){
	fn <- paste(bracken.folder, '/', sample.name, '/', sample.name, '_genus.bracken', sep = '')
	fread(fn)[name == 'Asaia', fraction_total_reads]
}

bracken.folder = '~/data/ML/GAARD_bracken/Korle-Bu_coluzzii_Delta'
kb.delta.samples <- list.dirs(bracken.folder, full.names = F) %>% .[. != '']

asaia <- setNames(nm = kb.delta.samples) %>%
         sapply(load.asaia, bracken.folder = bracken.folder)

# Load the phenotype data
phen.filename <- '../../../data/combined/sample_phenotypes.csv'
phenotypes <- fread(phen.filename, key = 'specimen')

# Load the Sanger submission data, which contains information on well position 
layouts.filename <- 'GAARD_WestAfrica_Sanger_submission.csv'
layouts <- fread(layouts.filename, key = 'SampleID')[kb.delta.samples, 
                                                     .(SampleID, GAARDPlate, GAARDWell, asaia)]

# Normalise the asaia values so the max is 1
layouts$asaia.norm <- layouts$asaia / max(layouts$asaia)
# Add the phenotypes
layouts$phenotype <- phenotypes[layouts$SampleID, phenotype]

# A function to plot values of a sample based on well position
well.plot <- function(wells, values, phenotypes, transform = 5, ...){
	well.x <- setNames(rep(1:12, each = 8), outer(LETTERS[1:8], 1:12, paste, sep = ''))
	well.y <- setNames(rep(8:1, 12), outer(LETTERS[1:8], 1:12, paste, sep = ''))
	phen.col <- c('red3', 'deepskyblue4')[(phenotypes == 'alive') + 1]
	if(transform)
		values <- values^(1/transform)
	plot(well.x, well.y, type = 'n', xlab = '', ylab = '', xaxt = 'n', yaxt = 'n')
	title(..., line = 2)
	points(well.x[wells], well.y[wells], cex = 3*(values), pch = 21, bg = lighten.col(phen.col, (values/2)+0.5), col = lighten.col(phen.col, values))
	missing.wells <- setdiff(names(well.x), wells)
	points(well.x[missing.wells], well.y[missing.wells], pch = 4, cex = 2, col = 'grey50', lwd = 2)
	axis(2, at = 8:1, labels = LETTERS[1:8], tick = F, las = 1)
	axis(3, at = 1:12, labels = 1:12, tick = F)
}

# Now a function to plot this for all the plates in a given layout
plate.well.plot <- function(plates, wells, values, phenotypes, transform = 5){
	w <- split(wells, plates)
	v <- split(values, plates)
	p <- split(phenotypes, plates)
	num.plates <- length(w)
	# Choose how many rows and columns for the display
	num.row.col <- matrix(c(1,1,1,2,2,2,3,3,3,3,3,3,1,2,3,2,3,3,3,3,3,4,4,4), 12, 2)
	if (num.plates <= 12)
		par(mfrow = num.row.col[num.plates, ])
	else { 
		num.plot.cols <- ceiling(sqrt(num.plates)) + 1
		num.plot.rows <- (num.plates %% num.plot.cols) + 1
		par(mfrow = c(num.plot.rows, num.plot.cols))
	}
	par(mar = c(1,2,4,1), mgp = c(2,0.5,0))
	for (pl in names(w))
		well.plot(w[[pl]], v[[pl]], p[[pl]], transform = transform, main = pl)
}

png('plate_layouts2.png', width = 8, height = 6, units = 'in', res = 300)
with(layouts, plate.well.plot(GAARDPlate, GAARDWell, asaia.norm, phenotype, transform = 5))
dev.off()


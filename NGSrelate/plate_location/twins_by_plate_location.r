library(data.table)
library(plotrix)
library(magrittr)

# The aim of this script is to correlate the presence of the bellwether SNPs and the significant Asaia reads
# to see whether they are both driven by the same signal. 

# We want to look for the position of twins as a function of plate position. For this, we need to know what well each 
# sample occupied on its original plate and on the sanger plate.

# Load the twins data
ngsrelate.filename <- '../ag3_gaard1244.3L.ngsRelate.tsv'
useful.columns <- c('specimen.x', 'location.x', 'species.x', 'specimen.y', 'location.y', 'species.y', 'KING')
ngsrelate <- fread(ngsrelate.filename, sep = '\t')[, ..useful.columns]
all.samples <- c(ngsrelate$specimen.x, ngsrelate$specimen.y) %>%
               unique() %>%
               sort()
twin.threshold <- 0.4
twins <- ngsrelate[KING >= twin.threshold, ]

# Load the metadata
meta.filename <- '/home/eric/Liverpool/ML/GAARD_SNP/pre_QC/metadata/sample_metadata.csv'
meta <- fread(meta.filename, key = 'External_ID')[all.samples, ]

# Add the phenotypes
phenotypes.filename <- '/home/eric/Desktop/Liverpool/ML/GAARD_SNP/pre_QC/metadata/sample_phenotypes.csv'
phenotypes.table <- fread(phenotypes.filename, key = 'specimen')

# There are 3 well layouts we need to deal with. First is the layout on our GAARD plates. Then there is the layout
# on the plates that we sent to Sanger. Finally, there is the layout that Sanger used when they combined plates. 
# The first two layouts can be found in the following files:
lstm.layouts.filename <- '/home/eric/Liverpool/ML/GAARD/data/GAARD_WestAfrica_Sanger_submission.csv'
lstm.layouts <- fread(lstm.layouts.filename, key = 'SampleID')
layout1 <- lstm.layouts[all.samples, .(SampleID, GAARDPlate, GAARDWell)]
layout2 <- lstm.layouts[all.samples, .(SampleID, WGSPlate, WGSWell)]
# The third layout can be found in the manifest that Chris Clarkson sent me
sanger.layout.filename <- '5563stdy_manifest_12930_210119.csv'
sanger.layout <- fread(sanger.layout.filename, key = 'SUPPLIER SAMPLE NAME')[meta$MalGEN_ID, ]
colnames(sanger.layout) <- gsub(' ', '_', colnames(sanger.layout))
setkey(meta, 'MalGEN_ID')
sanger.layout$SampleID <- meta[sanger.layout[['SUPPLIER_SAMPLE_NAME']], External_ID]
setkey(sanger.layout, SampleID)
setkey(meta, 'External_ID')
layout3 <- sanger.layout[, .(SampleID, SANGER_PLATE_ID, WELL)]
colnames(layout3) <- c('SampleID', 'SangerPlate', 'SangerWell')

layouts <- merge(merge(layout1, layout2), layout3)
layouts <- cbind(layouts, phenotypes.table[layouts$SampleID, .(Insecticide = insecticide, Phenotype = phenotype, Location = location)])

# A function to draw an arc
draw.arc <- function(x, y, ra, theta1, theta2){
	# Get the points required. First, get a points for a circle of radius 1 centered on 0
	Npoints <- 500
	basic.circle <- exp(pi * 1i * seq(0, 2, length.out = Npoints+1)[-1])
	# Then with the right radius
	rad.circle <- basic.circle * ra
	# And adjust the centre
	centered.circle <- rad.circle + x + y*1i
	# Now with the angle we want
	prop.t1 <- min(theta1, theta2) / (2*pi)
	prop.t2 <- max(theta1, theta2) / (2*pi)
	temp <- ((0:Npoints) / Npoints) 
	which.points <- which((temp > prop.t1) & (temp < prop.t2))
	arc.range <- centered.circle[which.points]
	lines(arc.range, col = 'magenta', lwd = 1.5, xpd = NA)
}

# A function to join to points by an arc
join.by.arc <- function(x1, y1, x2, y2, bow = 1){
	# The centre of the arc should sit on a line of slope
	sl <- (x1 - x2) / (y2 - y1)
	# that passes between the two points
	midpoint <- c(mean(c(x1, x2)), mean(c(y1, y2)))
	# The distance from that midpoint should be proportional to the distance between the original points. 
	di <- sqrt(sum(c(x2 - x1, y2 - y1)^2))
	# Where do you sit on that line to get a distance of di from midpoint?
	if (sl^2 == Inf){
		x.offset <- 0
		y.offset <- di * bow
	}
	else {
		x.offset <- bow * sqrt((di^2) / (1 + sl^2))
		y.offset <- x.offset * sl
	}
	# So where is the centre
	ce <- c(midpoint[1] + x.offset, midpoint[2] + y.offset)
	# And calculate the radius required to touch both points. 
	ra <- sqrt(sum(c(bow * di, di/2)^2))
	# Angle to the points
	theta <- function(p1, p2){
		x <- p2[1] - p1[1]
		y <- p2[2] - p1[2]
		theta <- atan2(y, x)
		if (theta < 0)
			theta <- 2*pi + theta
		theta
	}
	theta1 <- theta(ce, c(x1, y1))
	theta2 <- theta(ce, c(x2, y2))
	# Draw the arc
	draw.arc(ce[1], ce[2], ra, theta1, theta2)
}

# Write a function to plot twin connections between plates.
well.plot <- function(samples, wells, phenotypes, insecticide, twins, ...){
	well.x <- setNames(rep(1:12, each = 8), outer(LETTERS[1:8], 1:12, paste, sep = ''))
	well.y <- setNames(rep(8:1, 12), outer(LETTERS[1:8], 1:12, paste, sep = ''))
	phen.col <- c('red3', 'deepskyblue4')[(phenotypes == 'alive') + 1]
	insecticide.pch = c(21, 24)[(insecticide == 'Delta') + 1]
	plot(well.x, well.y, type = 'n', xlab = '', ylab = '', xaxt = 'n', yaxt = 'n')
	title(..., line = 3)
	points(well.x[wells], well.y[wells], cex = 2, pch = insecticide.pch, bg = phen.col, col = 'grey50')
	missing.wells <- setdiff(names(well.x), wells)
	points(well.x[missing.wells], well.y[missing.wells], pch = 4, cex = 2, col = 'grey50', lwd = 2)
	twins.present <- twins[specimen.x %in% samples | specimen.y %in% samples, ]
	# Keep track of number of off plates per sample
	off.plates <- setNames(numeric(length(samples)), samples)
	# The next line makes 1:nrow(twins.present) while dealing with the possibility of nrow = 0.
	for (i in seq(0, nrow(twins.present))[-1]){
		if (twins.present$specimen.x[i] %in% samples){
			if (twins.present$specimen.y[i] %in% samples){
				well1 <- wells[which(samples == twins.present$specimen.x[i])]
				well2 <- wells[which(samples == twins.present$specimen.y[i])]
				join.by.arc(well.x[well1], well.y[well1], well.x[well2], well.y[well2])
				#segments(well.x[well1], well.y[well1], well.x[well2], well.y[well2])
			}
			# Otherwise, the other well is on a different plate.
			else{
				# Increment by 1
				well1 <- wells[which(samples == twins.present$specimen.x[i])]
				x.offset <- off.plates[twins.present$specimen.x[i]] %>%
				            {((.+1) %/% 2) * ((. %% 2) - 0.5) * 0.5}
				points(well.x[well1] + x.offset, well.y[well1] - 0.2, pch = 19, cex = 0.8, col = 'orange')
				off.plates[twins.present$specimen.x[i]] <- off.plates[twins.present$specimen.x[i]] + 1
			}
		}
		else{
			well2 <- wells[which(samples == twins.present$specimen.y[i])]
			x.offset <- off.plates[twins.present$specimen.y[i]] %>%
						{((.+1) %/% 2) * ((. %% 2) - 0.5)}
			points(well.x[well2] + x.offset, well.y[well2] - 0.2, pch = 19, cex = 0.8, col = 'orange')
			off.plates[twins.present$specimen.x[i]] <- off.plates[twins.present$specimen.x[i]] + 1
		}
	}
	axis(2, at = 8:1, labels = LETTERS[1:8], tick = F, las = 1)
	axis(3, at = 1:12, labels = 1:12, tick = F)
}

# Now a function to plot this for all the plates in a given layout
plate.well.plot <- function(plates, samples, wells, phenotypes, insecticide, fn = ''){
	s <- split(samples, plates)
	w <- split(wells, plates)
	p <- split(phenotypes, plates)
	i <- split(insecticide, plates)
	num.plates <- length(w)
	num.plot.cols <- ceiling(sqrt(num.plates))
	num.plot.rows <- ((num.plates + num.plot.cols - 1) %/% num.plot.cols)
	if (fn == '')
		x11()
	else 
		png(fn, width = 4 * num.plot.cols, height = 3 * num.plot.rows, units = 'in', res = 300)
	par(mfrow = c(num.plot.rows, num.plot.cols), mar = c(1,3, 5, 1))
	for (pl in names(w))
		well.plot(s[[pl]], w[[pl]], p[[pl]], i[[pl]], main = pl, twins = twins)
	if (fn != '')
		dev.off()
}

layout.by.sampleset <- split(layouts, layouts$Location)
for (lay in names(layout.by.sampleset)){
	plot.fn1 <- paste(lay, '_GAARDPlate.png', sep = '')
	with(layout.by.sampleset[[lay]], plate.well.plot(GAARDPlate, SampleID, GAARDWell, Phenotype, Insecticide, plot.fn1))
	#
	plot.fn2 <- paste(lay, '_WGSPlate.png', sep = '')
	with(layout.by.sampleset[[lay]], plate.well.plot(WGSPlate, SampleID, WGSWell, Phenotype, Insecticide, plot.fn2))
	#
	plot.fn3 <- paste(lay, '_SangerPlate.png', sep = '')
	with(layout.by.sampleset[[lay]], plate.well.plot(SangerPlate, SampleID, SangerWell, Phenotype, Insecticide, plot.fn3))
}


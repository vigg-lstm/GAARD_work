library(data.table)
library(stringr)
library(fdrtool)
library(magrittr)

# Load the data
load('../../randomisations/Fst/fst_randomisations_Avrankou.Rdata')

# Use Obuasi gambiae Delta as an example
fst <- windowed.data.1000$Obuasi_gambiae_Delta$moving.Fst

# Use the negative part of the distribution to tell us how far to look in the positive part. 
find.outlier.threshold <- function(fst.vector, modal.increments = 3){
	fst.modal.bin <- which.max(table(cut(fst.vector, breaks = seq(min(fst.vector, na.rm = T), max(fst.vector, na.rm = T), 0.0005))))
	fst.mode <- mean(as.numeric(strsplit(sub(']$', '', sub('^\\(', '', names(fst.modal.bin))), ',')[[1]]))
	min.point <- min(fst.vector, na.rm = T)
	neg.diff <- fst.mode - min.point
	thresh <- fst.mode + neg.diff * modal.increments
	c(mode = fst.mode, min.point = min.point, thresh = thresh)
	
}

ot <- find.outlier.threshold(fst)

svg('Fst_example_histogram.svg', width = 5, height = 5)
layout(matrix(1:2, 2, 1), heights = c(3,1))
par(mar = c(0, 2.7, 0, 1))
hist.col = 'goldenrod'
hist(fst, breaks = 100, col = hist.col, xlab = '', xaxt = 'n', mgp = c(1.8, 0.7, 0.1), cex.axis = 0.6, main = '')
annot.offset = -25
segments(ot['mode'], annot.offset, ot['min.point'], annot.offset, col = 'firebrick2', lwd = 2)
segments(ot['thresh'], annot.offset, ot['mode'], annot.offset, col = 'firebrick2', lwd = 2, lty = 2)
points(ot['mode'], annot.offset, pch = 19, col = 'red', cex = 0.7)
par(mar = c(2.5, 2.7, 1, 1))
hist(fst, breaks = 100, col = hist.col, xlab = 'Fst', ylim = c(0,10), main = '', yaxt = 'n', ylab = '', mgp = c(1.4,0.5,0.1), cex.axis = 0.6)
hist(fst[fst > ot['thresh']], breaks = 50, col = 'dodgerblue', add = T)
axis(2, at = c(0, 5, 10), mgp = c(1.6,0.7,0.1), cex.axis = 0.6)
dev.off()




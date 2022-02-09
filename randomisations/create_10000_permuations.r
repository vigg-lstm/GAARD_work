# Add 9000 randomisations to the previously created 1000, for a total of 10000
library(data.table)
library(stringr)

# Load the phenotypes (original GAARD metadata) and the VObs metadata (provided by Sanger, contains only 
# sequenced samples. 
phenotypes <- fread('phenotype_randomisations.csv', key = 'specimen')

# A function that will return a list of length k, where each k is a random shuffling of the input vector
shuffle <- function(x, k)
	data.table(replicate(k, sample(x, length(x), replace = F)))

# Set a different seed to the original script, in order to avoid getting overlapping randomisations
set.seed(1)
num.randomisations <- 10000
replicate.names <- paste('r', str_pad(1:num.randomisations, nchar(as.character(num.randomisations)), pad = 0), sep = '')
# Add a digit to the old randomisation names to fit with the new nomenclature
colnames(phenotypes) <- sub('^r', 'r0', colnames(phenotypes))
# Add the new replicates
phenotypes[, (replicate.names[1001:num.randomisations]) := shuffle(phenotype, num.randomisations - 1000), by = population]

# Write the table to file
fwrite(phenotypes, 'phenotype_10000_randomisations.csv', sep = '\t', quote = F)

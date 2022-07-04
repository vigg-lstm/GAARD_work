library(data.table)
library(stringr)

# Load the phenotypes (original GAARD metadata) and the VObs metadata (provided by Sanger, contains only 
# sequenced samples. 
phenotypes <- fread('../data/combined/sample_phenotypes.csv', key = 'specimen')
vobs.meta <- fread('../data/combined/all_samples.samples.meta.csv', key = 'partner_sample_id')
# Remove males
vobs.meta <- vobs.meta[sex_call == 'F', ]

# Reduced the phenotype table to only include sequenced samples
phenotypes <- phenotypes[vobs.meta$partner_sample_id, ]

# Have a column to indicate population (location x species x insecticide)
phenotypes$population <- with(phenotypes, paste(location, species, insecticide, sep = '.'))

# Reorder and remove columns
column.order <- c('specimen', 'country', 'location', 'species', 'insecticide', 'population', 'phenotype')
phenotypes <- phenotypes[, ..column.order]

# A function that will return a list of length k, where each k is a random shuffling of the input vector
shuffle <- function(x, k)
	data.table(replicate(k, sample(x, length(x), replace = F)))

# Create 1000 randomisations of the phenotype labels, stratified by population, and add them to the phenotypes
# table
set.seed(42)
num.randomisations <- 1000
replicate.names <- paste('r', str_pad(1:num.randomisations, nchar(as.character(num.randomisations)), pad = 0), sep = '')
phenotypes[, (replicate.names) := shuffle(phenotype, num.randomisations), by = population]

# Write the table to file
fwrite(phenotypes, 'phenotype_randomisations.csv', sep = '\t', quote = F)

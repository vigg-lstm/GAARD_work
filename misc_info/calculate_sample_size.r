library(data.table)

meta <- fread('../data/combined/all_samples.samples.meta.csv', key = 'sample_id')
meta.sample.names <- meta$partner_sample_id
phen <- fread('../data/combined/sample_phenotypes.csv', key = 'specimen')[sort(meta.sample.names), ]

# Load the sib groups information
sib.groups <- fread('../NGSrelate/full_relatedness/sib_group_table.csv', sep = '\t')
sibs.to.remove <- sib.groups[keep == F, sample.name]

# Identify the males and add them to the list of samples to remove
males <- meta[sex_call == 'M', partner_sample_id]
samples.to.remove <- c(sibs.to.remove, males)

sample.sizes <- phen[, .(dead = sum(phenotype == 'dead'), alive = sum(phenotype == 'alive')), 
				       keyby = c('location', 'species', 'insecticide')]

nomales.phen <- phen[!(males)]
nomales.sample.sizes <- nomales.phen[, .(dead = sum(phenotype == 'dead'), alive = sum(phenotype == 'alive')), 
				                       keyby = c('location', 'species', 'insecticide')]

final.phen <- phen[!(samples.to.remove)]
final.sample.sizes <- final.phen[, .(dead = sum(phenotype == 'dead'), alive = sum(phenotype == 'alive')), 
				                   keyby = c('location', 'species', 'insecticide')]

cat('\nSample sizes before removing males and sibs:\n\n')
print(sample.sizes)

cat('\nSample sizes after removing males:\n\n')
print(nomales.sample.sizes)

cat('\n\nSample sizes after removing males and sibs:\n\n')
print(final.sample.sizes)

fwrite(nomales.sample.sizes, file = 'sample_sizes.csv', sep = '\t')
fwrite(final.sample.sizes, file = 'final_sample_sizes.csv', sep = '\t')


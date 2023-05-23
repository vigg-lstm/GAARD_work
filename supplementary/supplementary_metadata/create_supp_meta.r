library(data.table)
library(magrittr)

sequenced.meta <- fread('../../data/combined/all_samples.samples.meta.csv')[sex_call == 'F'] %>%
                  setnames('partner_sample_id', 'sample.id')
submission.meta <- fread('../../data/combined/sample_metadata.csv') %>%
                   setnames(c('External_ID', 'Date_of_Collection'), c('sample.id', 'collection.date'))
phen <- fread('../../data/combined/sample_phenotypes.csv') %>%
                  setnames(c('specimen', 'exposure_time'), c('sample.id', 'exposure.time'))
sibs <- fread('../../NGSrelate/full_relatedness/sib_group_table.csv') %>%
                  setnames('sample.name', 'sample.id')

output.table <- merge(sequenced.meta[, .(sample.id, country, location, latitude, longitude)],
                      submission.meta[, .(sample.id, collection.date)],
					  all = F) %>%
                merge(., phen[, !c('location', 'country', 'plate')],
					  all = F) %>%
				merge(., sibs[, .(sample.id, sib.group.id, exclude.as.sib = !keep)],
					  all = T)

output.table[is.na(exclude.as.sib), exclude.as.sib := F]
			
# Add accession numbers
accessions <- fread('../../data/sample_sequencing_runs.csv')

combined.accessions <- accessions[, .(sample.accession = unique(V13), 
                                      run.accession = paste(V14, collapse = ',')), 
                                    by = 'V1'] %>%
                       setnames('V1', 'sample.id')

output.table <- merge(output.table, combined.accessions, all.x = T, all.y = F, by = 'sample.id')

fwrite(output.table, 'metadata.csv', sep = '\t')


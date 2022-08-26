library(data.table)
library(magrittr)

sequenced.meta <- fread('../../data/combined/all_samples.samples.meta.csv')[sex_call == 'F'] %>%
                  setnames('partner_sample_id', 'sample.id')
submission.meta <- fread('../../data/combined/sample_metadata.csv') %>%
                   setnames(c('External_ID', 'Date_of_Collection'), c('sample.id', 'collection.date'))
phen <- fread('../../data/combined/sample_phenotypes.csv') %>%
                  setnames('specimen', 'sample.id')
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
			
fwrite(output.table, 'metadata.csv', sep = '\t')

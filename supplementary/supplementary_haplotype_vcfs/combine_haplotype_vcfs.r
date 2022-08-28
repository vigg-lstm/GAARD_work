library(data.table)
library(magrittr)

combine.haplotypes <- function(sample.set, hap.files){
	haplotype.vcfs <- grep(paste(sample.set, '.*\\.vcf', sep = ''), hap.files, value = T) 
	haplotype.csvs <- grep(paste(sample.set, '.*hap_sequence\\.csv', sep = ''), hap.files, value = T)
	allsites.vcf <- lapply(haplotype.vcfs, fread, select = c('#CHROM', 'POS', 'INFO')) %>%
	                rbindlist() %>%
	                setnames(c('#CHROM', 'POS', 'INFO'), c('Chrom', 'Pos', 'Info'))
	allsites.csv <- lapply(haplotype.csvs, fread) %>%
	                rbindlist(fill = T)
	full.table <- merge(allsites.csv, allsites.vcf)
	fwrite(full.table, file = paste(sample.set, '_haplotype_SNPs.csv', sep = ''), sep = '\t')
}

sample.sets <- c('Avrankou_coluzzii_Delta', 
                 'Baguida_gambiae_Delta', 
                 'Baguida_gambiae_PM', 
                 'Korle-Bu_coluzzii_Delta', 
                 'Korle-Bu_coluzzii_PM', 
                 'Madina_gambiae_Delta', 
                 'Madina_gambiae_PM', 
                 'Obuasi_gambiae_Delta', 
                 'Obuasi_gambiae_PM')

sapply(sample.sets, 
       combine.haplotypes, 
       hap.files = list.files('../../haplotypes/', full.names = T)
)



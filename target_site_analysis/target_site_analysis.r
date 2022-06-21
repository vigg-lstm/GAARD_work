library(reticulate)

use_condaenv('gaard')
zarr <- import('zarr')
expanduser <- import('os')$path$expanduser


# Second, pull out the allelic read depths from the zarr.
zarr_folder_avrankou <- '~/scratch/VObs_GAARD/1253-VO-TG-DJOGBENOU-VMF00052'
genotypes <- zarr$open(expanduser(paste(zarr_folder_avrankou, '/2R/calldata/GT', sep = '')), mode = 'r')
# Get the sampple names. First we get the VObs sample names from the zarr, then we convert to GAARD


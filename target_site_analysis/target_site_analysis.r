library(data.table)
library(magrittr)
library(reticulate)
use_condaenv('gaard')
zarr <- import('zarr')
np <- import('numpy', convert = F)
expanduser <- import('os')$path$expanduser

meta <- fread('../data/combined/sample_metadata.csv', key = 'MalGEN_ID')
phen <- fread('../data/combined/sample_phenotypes.csv', key = 'specimen')

pos <- list('2L' = py_to_r(np$array(zarr$open(expanduser('~/scratch/VObs_GAARD/sites/2L/variants/POS'), mode = 'r'))),
            '2R' = py_to_r(np$array(zarr$open(expanduser('~/scratch/VObs_GAARD/sites/2R/variants/POS'), mode = 'r'))),
            '3R' = py_to_r(np$array(zarr$open(expanduser('~/scratch/VObs_GAARD/sites/3R/variants/POS'), mode = 'r'))))

# Identify the index of the kdr SNPs
kdr.pos <- c('Vgsc.254K' = 2390177,
             'Vgsc.402L' = 2391228,
             'Vgsc.466M' = 2399997,
             'Vgsc.490I' = 2400071,
             'Vgsc.531V' = 2402466,
             'Vgsc.697P' = 2407967,
             'Vgsc.791M' = 2416980,
             'Vgsc.995S' = 2422651,
             'Vgsc.995F' = 2422652,
             'Vgsc.1507I' = 2429556,
             'Vgsc.1527T' = 2429617,
             'Vgsc.1570Y' = 2429745,
             'Vgsc.1597G' = 2429897,
             'Vgsc.1603T' = 2429915,
             'Vgsc.1746S' = 2430424,
             'Vgsc.1853I' = 2430817,
             'Vgsc.1868T' = 2430863,
             'Vgsc.1874S' = 2430880,
             'Vgsc.1874L' = 2430881,
             'Vgsc.1934V' = 2431061,
             'Vgsc.1940T' = 2431079
)
# In theory, we don't need "sort" here because we have already put the vector in index order, but since 
# the names will end up incorrectly associated with position if the markers are not in position order, we
# included it just in case.
kdr.pos <- sort(kdr.pos)
kdr.indices <- setNames(which(pos[['2L']] %in% kdr.pos), names(kdr.pos))

rdl.pos <- c('Rdl.296S' = 25429235,
             'Rdl.296G' = 25429236
)
rdl.pos <- sort(rdl.pos)
rdl.indices <- setNames(which(pos[['2L']] %in% rdl.pos), names(rdl.pos))


ace1.pos <- c('Ace1.119S' = 3492074)
ace1.pos <- sort(ace1.pos)
ace1.indices <- setNames(which(pos[['2R']] %in% ace1.pos), names(ace1.pos))

gste2.pos <- c('Gste2.114T' = 28598166,
              'Gste2.119V' = 28598062
)
gste2.pos <- sort(gste2.pos)
gste2.indices <- setNames(which(pos[['3R']] %in% gste2.pos), names(gste2.pos))


# Write a function to pull out genotypes at a list of positions
get.zarr.genotypes <- function(zarr_folder, chrom, indices){
	genotypes <- paste(zarr_folder, chrom, 'calldata/GT', sep = '/') %>%
	             expanduser() %>%
	             zarr$open(mode = 'r')
	vobs.sample.names <- zarr$open(expanduser(paste(zarr_folder, 'samples', sep = '/')), mode = 'r')$get_basic_selection()
	gaard.sample.names <- meta[sub('-.*', '', vobs.sample.names), External_ID]
	
	genotypes.array <- genotypes$get_orthogonal_selection(tuple(np_array(as.integer(indices - 1)), 
	                                                            np_array(as.integer((1:genotypes$shape[[2]]) - 1)), 
	                                                            np_array(as.integer(0:1)))
	                   )
	
	# For some reason, the genotypes array loads fine in reticulate, but the ALT alleles get loaded in the 
	# wrong row/column order (it outputs a N x 3 table, as it should, but the first row contains what should
	# be the first three values of the first column, etc...). I know that R and python use different row/column
	# prorities for indexing, but I don't know why the genotypes.array doesn't have this problem. Anyway, we
	# need to sort it out. I think this is the last time I'll use reticulate for zarr, too many issues that 
	# could end up flying under the radar. The last line of this block fixes the problem. 
	alt.alleles <- paste('~/scratch/VObs_GAARD/sites', chrom, 'variants/ALT', sep = '/') %>%
	               expanduser() %>%
	               zarr$open(mode = 'r') %>%
				   {.$get_orthogonal_selection(np_array(as.integer(indices - 1)))} %>%
				   matrix(., nrow(.), ncol(.), byrow = T)
	
	# The following block is a bit complicated as it involved nested magrittr pipes, but basically it goes 
	# through the loci, ignores ones where there are no mutants, splits ones where there are more than one
	# mutants, and spits out the resulting matrix
	extracted.genotypes <- 1:dim(genotypes.array)[1] %>%
	lapply(function(i){
		locus <- names(indices)[i]
		mutant.alleles <- table(genotypes.array[i,,]) %>%
		.[names(.) != '0']
		as.integer(names(mutant.alleles)) %>%
		lapply(function(ma){
			this.genotype <- rowSums(genotypes.array[i,,] == ma)
			this.name <- paste(locus, alt.alleles[i, ma], sep = '_')
			matrix(this.genotype, 1, length(this.genotype), dimnames = list(this.name, gaard.sample.names))
		}) %>%
		do.call(rbind, .)
	}) %>%
	do.call(rbind, .)
	
	extracted.genotypes
}

# Pull out the genotypes from the zarr.
zarr.folders <- c(Avrankou = '~/scratch/VObs_GAARD/1237-VO-BJ-DJOGBENOU-VMF00050',
                  Baguida = '~/scratch/VObs_GAARD/1253-VO-TG-DJOGBENOU-VMF00052',
                  Ghana = '~/scratch/VObs_GAARD/1244-VO-GH-YAWSON-VMF00051'
)

kdr.genotypes <- lapply(zarr.folders, get.zarr.genotypes, '2L', kdr.indices) %>%
                 Reduce(function(x, y) {xy <- merge(x, y, by = 'row.names', all = T); rownames(xy) <- xy[,1];  xy[,-1]}, .) %>%
                 .[sort(rownames(.)), sort(colnames(.))]
kdr.genotypes[is.na(kdr.genotypes)] <- 0

rdl.genotypes <- lapply(zarr.folders, get.zarr.genotypes, '2L', rdl.indices) %>%
                 Reduce(function(x, y) {xy <- merge(x, y, by = 'row.names', all = T); rownames(xy) <- xy[,1];  xy[,-1]}, .) %>%
                 .[sort(rownames(.)), sort(colnames(.))]
rdl.genotypes[is.na(rdl.genotypes)] <- 0

ace1.genotypes <- lapply(zarr.folders, get.zarr.genotypes, '2R', ace1.indices) %>%
                 Reduce(function(x, y) {xy <- merge(x, y, by = 'row.names', all = T); rownames(xy) <- xy[,1];  xy[,-1]}, .) %>%
                 .[sort(rownames(.)), sort(colnames(.))]
ace1.genotypes[is.na(ace1.genotypes)] <- 0

gste2.genotypes <- lapply(zarr.folders, get.zarr.genotypes, '3R', gste2.indices) %>%
                 Reduce(function(x, y) {xy <- merge(x, y, by = 'row.names', all = T); rownames(xy) <- xy[,1];  xy[,-1]}, .) %>%
                 .[sort(rownames(.)), sort(colnames(.))]
gste2.genotypes[is.na(gste2.genotypes)] <- 0


# Now calculate the allele frequency in each of the populations
wgs.phen <- phen[colnames(kdr.genotypes)]
wgs.phen$phenotype <- as.factor(wgs.phen$phenotype)
wgs.phen[, (rownames(kdr.genotypes)) := data.table(t(kdr.genotypes))]
wgs.phen[, (rownames(rdl.genotypes)) := data.table(t(rdl.genotypes))]
wgs.phen[, (rownames(ace1.genotypes)) := data.table(t(ace1.genotypes))]
wgs.phen[, (rownames(gste2.genotypes)) := data.table(t(gste2.genotypes))]
setkey(wgs.phen, location, insecticide)
all.markers <- c(rownames(kdr.genotypes),
                 rownames(rdl.genotypes),
                 rownames(ace1.genotypes),
                 rownames(gste2.genotypes))
pop.freqs <- wgs.phen[, lapply(.SD, function(x) mean(x)/2), by = .(location, species), .SDcols = all.markers]
t.pop.freqs <- cbind(data.table(all.markers),
                     data.table(pop.freqs[, t(.SD), .SDcols =  all.markers])) %>%
               setattr('names', c('Locus', paste(pop.freqs$location, pop.freqs$species, sep = '.')))

fwrite(t.pop.freqs, 'target_site_genotype_frequencies.csv', sep = '\t')

study.freqs <- wgs.phen[, lapply(.SD, function(x) mean(x)/2), by = .(location, species, insecticide), .SDcols = all.markers] %>%
               setkey(location, insecticide)

# GLMs:
glmodelling <- function(input.table, list.of.markers = markers, rescolumn = 'AliveDead', control.for = character(), glm.function = NULL, verbose = T){
	# Check whether the markers and random effects are present in the data.frame
	if (sum(list.of.markers %in% colnames(input.table)) != length(list.of.markers))
		stop('Some of the requested markers were not found in genotypes table.')
	if (sum(control.for %in% colnames(input.table)) != length(control.for))
		stop('At least one random effect was not found in genotypes table.')
	if (!(rescolumn %in% colnames(input.table)))
		stop('Resistance column not found in genotypes table.')
	# Remove any requested random effects that have only 1 level
	random.effects <- character()
	for (this.control in control.for){
		level.counts <- tapply(input.table[,this.control], input.table[,this.control], length)
		if (max(level.counts, na.rm = T) < sum(level.counts, na.rm = T)) 
			random.effects <- c(random.effects, this.control)
		else
			cat('Removing random effect ', this.control, ' was removed because it is invariable.')
	}
	# If you did not set a glm.function, decide which glm function you are going to use, based on whether mixed 
	# modelling will be necessary
	if (is.null(glm.function))
		glm.function <- ifelse(length(random.effects) > 0, 'glmer', 'glm')
	if (glm.function == 'glm'){
		# Create the random effect string, which will be empty if we have no random effects
		random.effects.string <- ifelse(length(random.effects) > 0, paste(' +', paste(random.effects, collapse = ' + ', sep = '')), '')
		# Set the name of the column containing the P.value in the anova function (which is different depending
		# on the glm function you use
		P.val.column <- 'Pr(>Chi)'
	}
	else if (glm.function == 'glmer'){
		random.effects.string <- ifelse(length(random.effects) > 0, paste(' +', paste('(1|', random.effects, ')', collapse = ' + ', sep = '')), '')
		P.val.column <- 'Pr(>Chisq)'
	}
	if (verbose)
		cat('Using the following string to control for confounding factors: ', random.effects.string, '\n', sep = '')
	# We remove markers for which there is no variation in the dataset or for which some alleles are too rare. 
	if (verbose)
		cat('\nDetecting invariable markers.\n')
	kept.markers <- character()
	invariable.markers <- character()
	for (this.marker in list.of.markers){
		allele.counts <- tapply(input.table[,this.marker], input.table[,this.marker], length)
		if (max(allele.counts, na.rm = T) <= (sum(allele.counts, na.rm = T) - 2)) 
			kept.markers <- c(kept.markers, this.marker)
		else
			invariable.markers <- c(invariable.markers, this.marker)
	}
	if (length(kept.markers) == 0)
		stop('Fail. None of the markers provided were variable.')
	# We check whether there are any ordered factors and recode them as numeric
	converted.table <- input.table
	has.ordered.factors <- F
	for (this.marker in kept.markers){
		if ('ordered' %in% class(converted.table[, this.marker])){
			if (verbose)
				cat('Converting ordered factor ', this.marker, ' to numeric.\n', sep = '')
			converted.table[, this.marker] <- as.numeric(converted.table[, this.marker])
			has.ordered.factors <- T
		}
	}
	# We do the glm analysis directly on the table from the global environment rather than the argument, this 
	# way the table that was used is recorded in the output. If we had to convert the ordered factors, then
	# we are forced to use a new table
	if (has.ordered.factors){
		working.table.name <- make.names(paste(deparse(substitute(input.table)), '_numeric_conversion', sep = ''))
		eval(parse(text = paste(working.table.name, '<- converted.table')))
	}
	else
		working.table.name <- make.names(deparse(substitute(input.table)))
	
	# For each marker, we calculate its pseudo-R2 and P-value compared to the null model.
	if (verbose)
		cat('\nAnalysing markers independently.\n')
	individual.markers <- data.frame(P = numeric(), pseudo.R2 = numeric())
	for (this.marker in kept.markers){
		# Remove the Na values for this marker
		this.table <- converted.table[!is.na(converted.table[,this.marker]),]
		# Build the model 
		this.model <- eval(parse(text = paste(glm.function, '(', rescolumn, ' ~ ', this.marker, random.effects.string, ', data = this.table, family = binomial)', sep = '')))
		# Build the null model
		this.null.model <- eval(parse(text = paste(glm.function, '(', rescolumn, ' ~ 1 ', random.effects.string, ', data = this.table, family = binomial)', sep = '')))
		# Get the stats
		this.p <- anova(this.model, this.null.model, test = 'Chisq')[[P.val.column]][2]
		# Report pseudo Rsquared if we used GLM and if we have the modEvA package
		if (('modEvA' %in% (.packages())) & (glm.function == 'glm'))
			this.pseudo.r2 <- mean(unlist(RsqGLM(this.model)))
		else
			this.pseudo.r2 <- NA
		individual.markers[this.marker, ] <- c(this.p, this.pseudo.r2)
	}
	
	# We now build the full model and remove markers one by one until all markers are significant
	markers.remaining <- kept.markers
	full.model.text <- paste(glm.function, '(', rescolumn, ' ~ ', paste(markers.remaining, collapse = ' + '), random.effects.string, ', data = ', working.table.name, ', family = binomial)', sep = '')
	# We'll keep track of the markers with perfect correlation
	correlated.markers <- character()
	if (verbose)
		cat('\nRunning commentary on model optimisation:\n\n')
	while(length(markers.remaining)){
		# It is possible that when building the model, some of the markers become monomorphic (because the other 
		# alleles are removed when other markers are NA). The model building will raise an error if this is the 
		# case, so let's check for this and temporarily remove such markers before we go ahead.
		temporarily.removed <- character()
		markers.remaining.genotypes <- input.table[, markers.remaining, drop = F]
		markers.remaining.genotypes <- markers.remaining.genotypes[complete.cases(markers.remaining.genotypes), , drop = F]
		for (this.marker in markers.remaining){
			allele.counts <- tapply(markers.remaining.genotypes[ ,this.marker], markers.remaining.genotypes[,this.marker], length)
			if (max(allele.counts, na.rm = T) == sum(allele.counts, na.rm = T)) {
				if (verbose){
					cat('Temporarily removing marker ', this.marker, ' as it has no variation left when NAs at other ',
					    'loci are removed.\n', sep = '')
				}
				temporarily.removed <- c(temporarily.removed, this.marker)
			}
		}
		# If any markers need to be temporarily removed, do so here.
		working.markers <- markers.remaining[!(markers.remaining %in% temporarily.removed)]
		# Build the model using the working.markers
		old.model.text <- paste(glm.function, '(', rescolumn, ' ~ ', paste(working.markers, collapse = ' + '), random.effects.string, ', data = ', working.table.name, ', family = binomial)', sep = '')
		if (verbose)
			cat('Building model:', old.model.text, '\n')
		old.model <- eval(parse(text = old.model.text))
		p.values <- numeric()
		# This will tell us whether we have broken out of the marker loop and therefore should go to the next model loop
		next.loop <- F
		for (this.marker in working.markers){
			new.model <- eval(parse(text = paste('update(old.model, .~.-', this.marker, ')', sep = '')))
			# Check that the new model doesn't have more rows than the old one (if the removed marker had unique NAs)
			if (length(fitted(new.model)) == length(fitted(old.model))){
				this.p.value <- anova(old.model, new.model, test = 'Chisq')[[P.val.column]][2]
			}
			# Otherwise, we need to rebuild the models with those samples removed.
			else{
				if (has.ordered.factors)
					reduced.input.table <- converted.table[!is.na(converted.table[,this.marker]),]
				else 
					reduced.input.table <- input.table[!is.na(input.table[,this.marker]),]
				temp.old.model <- eval(parse(text = paste(glm.function, '(', rescolumn, ' ~ ', paste(working.markers, collapse = ' + '), random.effects.string, ', data = reduced.input.table, family = binomial)', sep = '')))
				new.model <- eval(parse(text = paste('update(temp.old.model, .~.-', this.marker, ')', sep = '')))
				this.p.value <- anova(temp.old.model, new.model, test = 'Chisq')[[P.val.column]][2]
			}
			# If the p.value was NA, then there is perfect correlation or something else was wrong. Remove the marker with a comment
			if (is.na(this.p.value)){
				if (verbose)
					cat('\tRemoving marker ', this.marker, ' due to perfect correlation.\n\n', sep = '')
				markers.remaining <- markers.remaining[markers.remaining != this.marker]
				correlated.markers <- c(correlated.markers, this.marker)
				if (length(markers.remaining) == 0){
					final.model <- paste(glm.function, '(', rescolumn, ' ~ ', '1 ', random.effects.string, ', data = ', working.table.name, ', family = binomial)', sep = '')
					if (verbose)
						cat('All variables removed.\n')
				}
				next.loop <- T
				break
			}
			else{
				p.values[this.marker] <- this.p.value
			}
		}
		if (next.loop)
			next
		else if (max(p.values) > 0.05){
			# Remove the highest non-significant p-value
			if (verbose)
				cat('\tRemoving marker ', names(p.values)[which.max(p.values)], ' as the highest, non-significant marker (P = ', max(p.values), ').\n\n', sep = '')
			marker.to.remove <- names(p.values)[which.max(p.values)]
			markers.remaining <- markers.remaining[markers.remaining != marker.to.remove]
			if (length(markers.remaining) == 0){
				# If there are no markers remaining, then the final model is the null model
				final.model <- eval(parse(text = paste(glm.function, '(', rescolumn, ' ~ ', '1 ', random.effects.string, ', data = ', working.table.name, ', family = binomial)', sep = '')))
				if (verbose)
					cat('\tAll variables removed.\n')
			}
		}
		else {
			final.model <- old.model
			if (verbose){
				cat('\tAll remaining variables significant, keeping final model:\n')
				print(final.model)
			}
			break
		}
	}
	if (verbose){
		if (length(markers.remaining) == 0)
			cat('Final model was the null model.\n\n')
		else {
			cat('Final model contained ', length(markers.remaining), ' parameters: ', paste(markers.remaining, collapse = ','), '.\n\n', sep = '')
			if (length(temporarily.removed) > 0){
				cat('\tMarkers ', paste(temporarily.removed, collapse = ','), ' were temporarily removed but could ',
					'never be reintegrated into the model.\n', sep = '')
			}
		}
	}
	# Now get the p-values and pseudo R-squared value for all the variables in the final model, when added as the 
	# last variable
	if (length(markers.remaining) > 0){
		deviance.effect <- numeric()
		for (this.marker in names(p.values)){
			reduced.model <- eval(parse(text = paste('update(final.model, .~.-', this.marker, ')', sep = '')))
			deviance.effect[this.marker] <- deviance(reduced.model) - deviance(final.model)
		}
	}
	else {
		p.values <- NA
		deviance.effect <- NA
		temporarily.removed <- NA
	}
	cat('\n')
	#
	final.model.sig <- data.frame(P = p.values, deviance = deviance.effect)
	print(final.model.sig)
	list('invariable.markers' = invariable.markers, 'correlated.markers' = correlated.markers, 'sig.alone' = individual.markers, 'final.model' = final.model, 'final.sig' = final.model.sig, 'temporarily.removed' = temporarily.removed)
}

# Let's test things in Avrankou Delta
cat('\n\nAvrankou Delta:\n')
avrankou.delta.seg.markers <- study.freqs[.('Avrankou', 'Delta'), ..all.markers] %>%
                              {(. > 0.05) & (. < 0.95)} %>%
                              colnames(.)[.]
avrankou.delta.table <- as.data.frame(wgs.phen[.('Avrankou', 'Delta')])
glmodelling(avrankou.delta.table, avrankou.delta.seg.markers, 'phenotype')
cat('\nNote that the direction of effect of 995F, though non-significant, is that 995F makes you less',
    'resistant (ie: 402L affords more resistance:\n')
print(glm(phenotype ~ Vgsc.995F_T, data = avrankou.delta.table, family = 'binomial')$coefficients)

# Baguida Delta
cat('\n\nBaguida Delta:\n')
baguida.delta.seg.markers <- study.freqs[.('Baguida', 'Delta'), ..all.markers] %>%
                        {(. > 0.05) & (. < 0.95)} %>%
                        colnames(.)[.]
baguida.delta.table <- as.data.frame(wgs.phen[.('Baguida', 'Delta')])
glmodelling(baguida.delta.table, baguida.delta.seg.markers, 'phenotype')
cat('\nNothing is significant in Baguida Delta\n')

# Korle-Bu Delta
cat('\n\nKorle-Bu Delta:\n')
kb.delta.seg.markers <- study.freqs[.('Korle-Bu', 'Delta'), ..all.markers] %>%
                        {(. > 0.05) & (. < 0.95)} %>%
                        colnames(.)[.]
# We kick out 402L_C 
kb.delta.seg.markers <- setdiff(kb.delta.seg.markers, c('Vgsc.402L_C'))
kb.delta.table <- as.data.frame(wgs.phen[.('Korle-Bu', 'Delta')])
glmodelling(kb.delta.table, kb.delta.seg.markers, 'phenotype')
cat('\nNothing is significant in Korle-Bu Delta\n')

# Madina Delta
cat('\n\nMadina Delta:\n')
madina.delta.seg.markers <- study.freqs[.('Madina', 'Delta'), ..all.markers] %>%
                        {(. > 0.05) & (. < 0.95)} %>%
                        colnames(.)[.]
madina.delta.table <- as.data.frame(wgs.phen[.('Madina', 'Delta')])
glmodelling(madina.delta.table, madina.delta.seg.markers, 'phenotype')
cat('\nGste2-119V negatively associated with Deltamethrin resistance in Madina, but only just.\n')

# Obuasi Delta
cat('\n\nObuasi Delta:\n')
obuasi.delta.seg.markers <- study.freqs[.('Obuasi', 'Delta'), ..all.markers] %>%
                        {(. > 0.05) & (. < 0.95)} %>%
                        colnames(.)[.]
obuasi.delta.table <- as.data.frame(wgs.phen[.('Obuasi', 'Delta')])
glmodelling(obuasi.delta.table, obuasi.delta.seg.markers, 'phenotype')
cat('\nVgsc.791M and Gste2-119V positively associated with Deltamethrin resistance in Obuasi. Vgsc.1746S negatively associated. None are particularly strongly significant.\n')

# Baguida PM
cat('\n\nBaguida PM:\n')
baguida.pm.seg.markers <- study.freqs[.('Baguida', 'PM'), ..all.markers] %>%
                        {(. > 0.05) & (. < 0.95)} %>%
                        colnames(.)[.]
baguida.pm.table <- as.data.frame(wgs.phen[.('Baguida', 'PM')])
glmodelling(baguida.pm.table, baguida.pm.seg.markers, 'phenotype')
cat('\nVgsc.1570Y positively associated with PM resistance in Baguida.\n')

# Korle-Bu PM
cat('\n\nKorle-Bu PM:\n')
kb.pm.seg.markers <- study.freqs[.('Korle-Bu', 'PM'), ..all.markers] %>%
                        {(. > 0.05) & (. < 0.95)} %>%
                        colnames(.)[.]
kb.pm.table <- as.data.frame(wgs.phen[.('Korle-Bu', 'PM')])
glmodelling(kb.pm.table, kb.pm.seg.markers, 'phenotype')
cat('\nAce1 very strongly associated with PM resistance in Korle-Bu. Vgsc-995F (ie: absence of 402L) and Gste2-114T both negatively associated.\n')

# Madina PM
cat('\n\nMadina PM:\n')
madina.pm.seg.markers <- study.freqs[.('Madina', 'PM'), ..all.markers] %>%
                        {(. > 0.05) & (. < 0.95)} %>%
                        colnames(.)[.]
madina.pm.table <- as.data.frame(wgs.phen[.('Madina', 'PM')])
glmodelling(madina.pm.table, madina.pm.seg.markers, 'phenotype')
cat('\nAce1 very strongly associated with PM resistance in Madina. Vgsc-1853Y negatively associated.\n')

# Obuasi PM
cat('\n\nObuasi PM:\n')
obuasi.pm.seg.markers <- study.freqs[.('Obuasi', 'PM'), ..all.markers] %>%
                        {(. > 0.05) & (. < 0.95)} %>%
                        colnames(.)[.]
obuasi.pm.table <- as.data.frame(wgs.phen[.('Obuasi', 'PM')])
glmodelling(obuasi.pm.table, obuasi.pm.seg.markers, 'phenotype')
cat('\nAce1 very strongly associated with PM resistance in Obuasi.\n')

# Now let's try comining locations:
pm.table <- as.data.frame(wgs.phen[.('Obuasi', 'PM')])
# We kick out 402L since it's perfectly associated with 995 and the presence of two alleles is confusing
non.402.markers <- setdiff(all.markers, c('Vgsc.402L_C', 'Vgsc.402L_T')) 
glmodelling(pm.table, non.402.markers, 'phenotype', control.for = 'location')
cat('\nAce1 is the only significant SNP. It is positively associated with resistance.\n')

delta.table <- as.data.frame(wgs.phen[.('Obuasi', 'Delta')])
glmodelling(delta.table, non.402.markers, 'phenotype', control.for = 'location')
cat('\nVgsc-791M and Gste2-119V are both positively associated with deltamethrin resistance. Vgsc-1746S is negatively associated.')
cat('\nNone of the P-value are particularly small.\n')



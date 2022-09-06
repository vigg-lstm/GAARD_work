#!/usr/bin/env python
# coding: utf-8

"""
"""

import sys
sys.stderr = open(snakemake.log[0], "w")

import probetools as probe
import numpy as np
import pandas as pd
import allel
import dask.array as da
from pathlib import Path

randomisations = pd.read_csv("phenotype_randomisations.csv", sep="\t")

# Garuds Selection Scans # 
cloud = False #snakemake.params['cloud']
#ag3_sample_sets = "" #snakemake.params['ag3_sample_sets']
stat = 'H12' #snakemake.params['GarudsStat']
windowSize = 1200 #snakemake.params['windowSize']
windowStep = 600 #snakemake.params['windowStep']
contig = snakemake.params['contig']
cohort = snakemake.params['cohort']

#cohort = '.'.join(np.array(cohort_out.split("."))[[2,0,1]])

haplotypePath = f"/home/snagi/lstm_projects/VObs_GAARD/ag3_gaard_haplotypes/{contig}/GT/"
positionsPath=  f"/home/snagi/lstm_projects/VObs_GAARD/ag3_gaard_haplotypes/sites/{contig}/POS/" #snakemake.input['haplotypes'] if stat in ['H1', 'H12', 'H2/1'] else []


metadata = pd.read_csv("~/lstm_projects/VObs_GAARD/ag3_gaard_haplotypes/gaard_metadata_haplotypes.tsv", sep="\t")
female_bool = metadata.eval("sex_call == 'F'")
meta = metadata.query("sex_call == 'F'")
sibs = pd.read_csv("resources/sib_group_table.csv", sep="\t")
exclude = sibs.query("keep == False")['sample.name']
sibs_bool = meta.eval("partner_sample_id not in @exclude")
metadata = meta.query("partner_sample_id not in @exclude").reset_index(drop=True)


haps, pos = probe.loadZarrArrays(haplotypePath, positionsPath, siteFilterPath=None, haplotypes=True, cloud=cloud, contig=contig)
haps = allel.HaplotypeArray(haps).to_genotypes(ploidy=2).compress(female_bool, axis=1).compress(sibs_bool, axis=1).to_haplotypes()


def garudsStat(stat, geno, pos, cut_height=None, metric='euclidean', window_size=1200, step_size=600):
    
    """
    Calculates G12/G123/H12
    """
        
    # Do we want to cluster the Multi-locus genotypes (MLGs), or just group MLGs if they are identical
    if stat == "G12":
        garudsStat = allel.moving_statistic(geno, clusterMultiLocusGenotypes, size=window_size, step=step_size, metric=metric, cut_height=cut_height, g=2)
    elif stat == "G123":
        garudsStat = allel.moving_statistic(geno, clusterMultiLocusGenotypes, size=window_size, step=step_size, metric=metric, cut_height=cut_height, g=3)
    elif stat == "H12":
        garudsStat,_,_,_ = allel.moving_garud_h(geno, size=window_size, step=step_size)
    else:
        raise ValueError("Statistic is not G12/G123/H12")

    midpoint = allel.moving_statistic(pos, np.median, size=window_size, step=step_size)
    
    return(garudsStat, midpoint)



#### Load cohort data and their indices in genotype data
### run garudStat for that query. already loaded contigs 

#cohorts = probe.getCohorts(metadata=metadata, 
#                    columns=['species', 'phenotype', 'location', 'insecticide'], 
#                    minPopSize=15, excludepath="resources/sib_group_table.csv")

random = randomisations.query("population == @cohort")

alive_dead_dict = {}
for i in np.arange(1, 200):
    i = '%0.4d' % i

    alive_dead_dict['alive'] = np.where(pd.factorize(random[f'r{i}'])[0])[0]
    alive_dead_dict['dead']  = np.where(pd.factorize(random[f'r{i}'])[0] == False)[0]

    for pheno, idxs in alive_dead_dict.items():
# Loop through each cohort, manipulate genotype arrays and calculate chosen Garuds Statistic

        if Path("randomisations/H12/H12_{cohort}_{pheno}{i}.{contig}.tsv").exists():
            probe.log(f"skipping {cohort} {pheno} {i} {contig}")
            continue

        if stat in ['H1', 'H12', 'H123']:
            # get indices for haplotype Array and filter
            hapInds = np.sort(np.concatenate([idxs*2, (idxs*2)+1]))
            gt_cohort = haps.take(hapInds, axis=1)
        else:
            raise ValueError("Statistic is not G12/G123/H1/H12")

        probe.log(f"--------- Running {stat} on {cohort} | Chromosome {contig} ----------")
        probe.log("filter to biallelic segregating sites")

        ac_cohort = gt_cohort.count_alleles(max_allele=3)
        # N.B., if going to use to_n_alt later, need to make sure sites are 
        # biallelic and one of the alleles is the reference allele
        ref_ac = ac_cohort[:, 0]
        loc_sites = ac_cohort.is_biallelic() & (ref_ac > 0)
        gt_seg = da.compress(loc_sites, gt_cohort, axis=0)
        pos_seg = da.compress(loc_sites, pos, axis=0)

        probe.log(f"compute input data for {stat}")
        pos_seg = pos_seg.compute(numworkers=1)

        # calculate G12/G123/H12 and plot figs 
        gStat, midpoint = garudsStat(stat=stat,
                                    geno=gt_seg, 
                                    pos=pos_seg, 
                                    cut_height=" ",
                                    metric='euclidean',
                                    window_size=windowSize,
                                    step_size=windowStep)

        probe.windowedPlot(statName=stat, 
                    cohortText = cohort + "_" + pheno + str(i),
                    cohortNoSpaceText= cohort + "_" + pheno +str(i),
                    values=gStat, 
                    midpoints=midpoint,
                    prefix=f"randomisations/{stat}", 
                    contig=contig,
                    colour="dodgerblue",
                    ymin=0,
                    ymax=1)

#!/usr/bin/env python
# coding: utf-8


from probetools import *
import pandas as pd 

import sys
sys.stderr = open(snakemake.log[0], "w")
from pathlib import Path


randomisations = pd.read_csv("phenotype_randomisations.csv", sep="\t")

meta = pd.read_csv("/home/snagi/lstm_data/gaard-probe/config/gaard_metadata.tsv", sep="\t")
female_bool = meta.eval("sex_call == 'F'")
meta = meta.query("sex_call == 'F'")
sibs = pd.read_csv("resources/sib_group_table.csv", sep="\t")
exclude = sibs.query("keep == False")['sample.name']
sibs_bool = meta.eval("partner_sample_id not in @exclude")
meta = meta.query("partner_sample_id not in @exclude")
contig = snakemake.params['contig']
cohort = snakemake.params['cohort']

#cohort = '.'.join(np.array(cohort_out.split("."))[[2,0,1]])

# Outgroup data
outgroupMetaPath = f"resources/AG1000G-ML-B/samples.meta.csv"
Mali2004Meta = pd.read_csv(outgroupMetaPath)
species = pd.read_csv("resources/AG1000G-ML-B/samples.species_aim.csv")
Mali2004Meta = Mali2004Meta.merge(species)


stat = 'PBS'
random = randomisations.query("population == @cohort")
spp = random['species'].unique()[0]

print(contig, f"PBS Running {cohort}")

# Load arrays
snps, pos = loadZarrArrays(genotypePath=f"/home/snagi/lstm_projects/VObs_GAARD/ag3_gaard/{contig}/GT/",
                            positionsPath=f"/home/snagi/lstm_projects/VObs_GAARD/sites/{contig}/variants/POS", 
                            siteFilterPath=f"/home/snagi/lstm_projects/VObs_GAARD/site_filters/dt_20200416/gamb_colu/{contig}/variants/filter_pass/")

# filter out sibs
snps = snps.compress(female_bool, axis=1)
snps = snps.compress(sibs_bool, axis=1)

### Load outgroup Arrays and subset to each species, storing
snpsOutgroup, pos = loadZarrArrays(genotypePath=f"resources/AG1000G-ML-B/{contig}/calldata/GT/", 
                                positionsPath=f"/home/snagi/lstm_projects/VObs_GAARD/sites/{contig}/variants/POS/", 
                                siteFilterPath=f"/home/snagi/lstm_projects/VObs_GAARD/site_filters/dt_20200416/gamb_colu/{contig}/variants/filter_pass/")
snpsOutgroupDict = {}

for sp in ['gambiae', 'coluzzii']:
    sp_bool = Mali2004Meta['species_gambiae_coluzzii'] == sp
    snpsOutgroupDict[sp] =  snpsOutgroup.compress(sp_bool, axis=1)

for i in np.arange(1,200):
    
    i = '%0.4d' % i
    
    if Path("randomisations/PBS/PBS_PBS.{cohort}.{i}.{contig}.tsv").exists():
        print(f"skipping {cohort} {i} {contig}")
        continue

    alive = np.where(pd.factorize(random[f'r{i}'])[0])[0]
    dead = np.where(pd.factorize(random[f'r{i}'])[0] == False)[0]

    ac_cohort = snps.count_alleles(max_allele=3).compute(numworkers=4)
    # N.B., if going to use to_n_alt later, need to make sure sites are 
    # biallelic and one of the alleles is the reference allele
    ref_ac = ac_cohort[:, 0]
    loc_sites = ac_cohort.is_biallelic() & (ref_ac > 0)
    gt_seg = da.compress(loc_sites, snps, axis=0)
    pos_seg = da.compress(loc_sites, pos, axis=0)

    log(f"compute input data for {stat}")
    pos_seg = pos_seg.compute(numworkers=4)

    ac_out = allel.GenotypeArray(snpsOutgroupDict[spp].compress(loc_sites, axis=0)).count_alleles()
    ac_pheno1 = allel.GenotypeArray(gt_seg).take(alive, axis=1).count_alleles()
    ac_pheno2 = allel.GenotypeArray(gt_seg).take(dead, axis=1).count_alleles()

    assert ac_out.shape[0] == pos_seg.shape[0], "Array Outgroup/POS are the wrong length"
    assert ac_pheno1.shape[0] == pos_seg.shape[0], "Array phenotype1/POS are the wrong length"
    assert ac_pheno2.shape[0] == pos_seg.shape[0], "Arrays phenotype2/POS the wrong length"

    log("calculate PBS and plot figs")
    # calculate PBS and plot figs 
    pbsArray = allel.pbs(ac_pheno1, ac_pheno2, ac_out, 
                window_size=1000, window_step=500, normed=True)
    midpoint = allel.moving_statistic(pos_seg, np.mean, size=1000, step=500)

    

    windowedPlot(statName=stat, 
                cohortText = f"PBS {cohort} {i}",
                cohortNoSpaceText= f"PBS.{cohort}.{i}",
                values=pbsArray, 
                midpoints=midpoint,
                prefix=f"randomisations/{stat}", 
                contig=contig,
                colour="Grey",
                ymin=-0.3,
                ymax=0.3,
                save=True)
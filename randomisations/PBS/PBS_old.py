from sys import argv
import numpy as np
import pandas as pd 
from pathlib import Path
import malariagen_data
import allel

if (len(argv) == 4):
    cohort = argv[1]
    contig = argv[2]
    num_randomisations = int(argv[3])
else:
    raise Exception("Fail. There should be three command line arguments (cohort, contig, num_randomisations).")

ag3 = malariagen_data.Ag3(pre=True)
metadata = ag3.sample_metadata("3.2")
metadata = metadata.query("sex_call == 'F'")
sibs = pd.read_csv("../../NGSrelate/full_relatedness/sib_group_table.csv", sep="\t")
exclude = sibs.query("keep == False")['sample.name']
metadata = metadata.query("partner_sample_id not in @exclude").reset_index(drop=True)

randomisations = pd.read_csv("../phenotype_randomisations.csv", sep="\t")
#### make sure all metadata samples are in the randomisations file (so GAARD data)
metadata = metadata.query(f"partner_sample_id in {randomisations['specimen'].to_list()}")
ids = metadata['partner_sample_id']
# make sure randomisations have only samples in the metadata, so excluding sibs and females
random = randomisations.query("specimen in @ids")
spp = random.query("population == @cohort")['species'].unique()[0]
print(f"Species of this cohort is {spp}.")
cohort_ids = random.query(f"population == @cohort")['specimen'].to_list()

print(f"Downloading genotypes for {cohort} {contig}")
ds_snps = ag3.snp_calls(
    region=contig,
    site_mask="gamb_colu",
    sample_query=f"partner_sample_id in {cohort_ids}",
    sample_sets=None,
)

gt = allel.GenotypeArray(ds_snps["call_genotype"].data)

print(f"Getting biallelic mask for {cohort} {contig}")
ac = gt.count_alleles()
biallelic_mask = ac.is_biallelic()

print(f"removing biallelics for {cohort} {contig}, {biallelic_mask.sum()} snps, array is shape: {ds_snps['call_genotype'].shape}")
ds_snps = ds_snps.isel(variants=biallelic_mask)
print(f" for {cohort} {contig} after, array is shape: {ds_snps['call_genotype'].shape}")
gt = gt[biallelic_mask, :]
metadata.set_index('sample_id', inplace = True)
snp_ids = metadata.loc[ds_snps.sample_id,'partner_sample_id'].values

pos = ds_snps["variant_position"].values
midpoints = allel.moving_statistic(pos, statistic=np.mean, size=1000)

outgroup_query=f"country == 'Mali' & year == 2004 & taxon == '{spp}'"
ds_snps_outgroup = ag3.snp_calls(
    region=contig,
    site_mask='gamb_colu',
    sample_query=outgroup_query,
    sample_sets=None,
).isel(variants=biallelic_mask)

gt_outgroup = allel.GenotypeArray(ds_snps_outgroup["call_genotype"].data)
ac_outgroup = gt_outgroup.count_alleles()

alive_ids = random.query(f"phenotype == 'alive' & population == @cohort")['specimen'].to_list()
dead_ids = random.query(f"phenotype == 'dead' & population == @cohort")['specimen'].to_list()

cohort1_gt = gt[:, np.isin(snp_ids, alive_ids)]
cohort2_gt = gt[:, np.isin(snp_ids, dead_ids)]
ac1 = cohort1_gt.count_alleles()
ac2 = cohort2_gt.count_alleles()

print(f"--------- Running PBS on {cohort} | Chromosome {contig} ----------")
pbs = allel.pbs(ac1=ac1, ac2=ac2, ac3=ac_outgroup, window_size=1000, normed=True)

pbs_df = pd.DataFrame({'midpoint':midpoints, 'pbs':pbs})
pbs_df.to_csv(f"PBS_outputs/PBS_{cohort}.{contig}.tsv", sep="\t")

for i in np.arange(1,num_randomisations+1):
    i = '%0.4d' % i
    
    if Path(f"PBS_outputs/PBS_{cohort}_{i}.{contig}.tsv").exists():
        print(f"skipping {cohort} {i} {contig}")
        continue

    alive_ids = random.query(f"r{i} == 'alive' & population == @cohort")['specimen'].to_list()
    dead_ids = random.query(f"r{i} == 'dead' & population == @cohort")['specimen'].to_list()

    cohort1_gt = gt[:, np.isin(snp_ids, alive_ids)]
    cohort2_gt = gt[:, np.isin(snp_ids, dead_ids)]
    ac1 = cohort1_gt.count_alleles()
    ac2 = cohort2_gt.count_alleles()

    print(f"--------- Running PBS on {cohort} randomisation {i} | Chromosome {contig} ----------")
    pbs = allel.pbs(ac1=ac1, ac2=ac2, ac3=ac_outgroup, window_size=1000, normed=True)

    pbs_df = pd.DataFrame({'midpoint':midpoints, 'pbs':pbs})
    pbs_df.to_csv(f"PBS_outputs/PBS_{cohort}_{i}.{contig}.tsv", sep="\t")


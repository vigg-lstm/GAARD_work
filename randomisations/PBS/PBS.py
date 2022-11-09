from sys import argv
import numpy as np
import pandas as pd
import malariagen_data
import allel
from sys import stdout 
from pathlib import Path

if (len(argv) == 4):
    cohort = argv[1]
    contig = argv[2]
    num_randomisations = int(argv[3])
else:
    raise Exception("Fail. There should be three command line arguments (cohort, contig, num_randomisations).")

window_size = 1000

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
stdout.flush()
cohort_ids = random.query(f"population == @cohort")['specimen'].to_list()

print(f"Downloading SNP data for {cohort}")
stdout.flush()
ds_snps = ag3.snp_calls(
    region=contig,
    site_mask="gamb_colu",
    sample_query=f"partner_sample_id in {cohort_ids}",
    sample_sets=None,
)
gt = allel.GenotypeArray(ds_snps["call_genotype"].data)

print(f"Identifying non-segregating alleles in {cohort}")
stdout.flush()
ac = gt.count_alleles()
segregating_mask = ac.is_segregating()

print("Downloading SNP data for Mali as outgroup")
stdout.flush()
outgroup_query=f"country == 'Mali' & year == 2004 & taxon == '{spp}'"
ds_snps_outgroup = ag3.snp_calls(
    region=contig,
    site_mask='gamb_colu',
    sample_query=outgroup_query,
    sample_sets=None,
).isel(variants=segregating_mask)

print("Identifying non-segregating alleles in Mali")
stdout.flush()
gt_outgroup = allel.GenotypeArray(ds_snps_outgroup["call_genotype"].data)
outgroup_segregating_mask = gt_outgroup.count_alleles().is_segregating()
gt_outgroup = gt_outgroup[outgroup_segregating_mask]
ac_outgroup = gt_outgroup.count_alleles()

ds_snps = ds_snps.isel(variants=segregating_mask).isel(variants=outgroup_segregating_mask)
gt = gt[segregating_mask][outgroup_segregating_mask]
metadata.set_index('sample_id', inplace = True)
gaard_ids = metadata.loc[ds_snps.sample_id,'partner_sample_id'].values

pos = ds_snps["variant_position"].values
midpoints = allel.moving_statistic(pos, statistic=np.mean, size=window_size)
startpoints = allel.moving_statistic(pos, statistic=np.min, size=window_size)
endpoints = allel.moving_statistic(pos, statistic=np.max, size=window_size)

alive_ids = random.query(f"phenotype == 'alive' & population == @cohort")['specimen'].to_list()
dead_ids = random.query(f"phenotype == 'dead' & population == @cohort")['specimen'].to_list()

cohort1_gt = gt[:, np.isin(gaard_ids, alive_ids)]
cohort2_gt = gt[:, np.isin(gaard_ids, dead_ids)]
ac1 = cohort1_gt.count_alleles()
ac2 = cohort2_gt.count_alleles()

print(f"--------- Running PBS on {cohort} | Chromosome {contig} ----------")
stdout.flush()
pbs = allel.pbs(ac1=ac1, ac2=ac2, ac3=ac_outgroup, window_size=window_size, normed=True)
true_pbs_df = pd.DataFrame({'startpoint': startpoints, 
                            'endpoint': endpoints,
                            'midpoint':midpoints,
                            'pbs':np.round(pbs, 3)})
true_pbs_df.to_csv(f"PBS_outputs/PBS_{cohort}_true.{contig}.tsv", sep="\t", index = False)

# Now the randomisations
randomised_pbs = {}
for i in np.arange(1,num_randomisations+1):
    i = '%0.4d' % i
        
    if Path(f"PBS_outputs/PBS_{cohort}_{i}.{contig}.tsv").exists():
        print(f"skipping {cohort} {i} {contig}")
        continue

    alive_ids = random.query(f"r{i} == 'alive' & population == @cohort")['specimen'].to_list()
    dead_ids = random.query(f"r{i} == 'dead' & population == @cohort")['specimen'].to_list()
        
    cohort1_gt = gt[:, np.isin(gaard_ids, alive_ids)]
    cohort2_gt = gt[:, np.isin(gaard_ids, dead_ids)]
    ac1 = cohort1_gt.count_alleles()
    ac2 = cohort2_gt.count_alleles()
        
    print(f"--------- Running PBS on {cohort} randomisation {i} | Chromosome {contig} ----------")
    stdout.flush()
    pbs = allel.pbs(ac1=ac1, ac2=ac2, ac3=ac_outgroup, window_size=window_size, normed=True)
    pbs_df = pd.DataFrame({'midpoints':midpoints, 'pbs':np.round(pbs, 3)})
    pbs_df.to_csv(f"PBS_outputs/PBS_{cohort}_{i}.{contig}.tsv", sep="\t")
    randomised_pbs[f"r{i}"] = np.round(pbs, 3)
    
randomised_pbs_df = pd.DataFrame(randomised_pbs)
pbs_df = pd.concat([true_pbs_df, randomised_pbs_df], axis = 1)

pbs_df.to_csv(f"PBS_outputs/PBS_{cohort}.{contig}.tsv", sep="\t", index = False)
    
    
    

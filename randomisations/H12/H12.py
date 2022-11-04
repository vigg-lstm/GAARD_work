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

def h12_gwss(
        ht,
        pos,
        window_size,
    ):

        h1, h12, h123, h2_h1 = allel.moving_garud_h(ht, size=window_size)

        x = allel.moving_statistic(pos, statistic=np.mean, size=window_size)

        return x, h12


ag3 = malariagen_data.Ag3(pre=True)
metadata = ag3.sample_metadata("3.2")
meta = metadata.query("sex_call == 'F'")
sibs = pd.read_csv("../../NGSrelate/full_relatedness/sib_group_table.csv", sep="\t")
exclude = sibs.query("keep == False")['sample.name']
metadata = meta.query("partner_sample_id not in @exclude").reset_index(drop=True)

randomisations = pd.read_csv("../phenotype_randomisations.csv", sep="\t")
#### make sure all metadata samples are in the randomisations file (so GAARD data)
metadata = metadata.query(f"partner_sample_id in {randomisations['specimen'].to_list()}")
ids = metadata['partner_sample_id']
# make sure randomisations have only samples in the metadata, so excluding sibs and females
random = randomisations.query("specimen in @ids")
cohort_ids = random.query(f"population == @cohort")['specimen'].to_list()

print(f"Getting biallelic mask for {cohort} {contig}")
ds_haps = ag3.haplotypes(
    region=contig,
    analysis="gamb_colu",
    sample_query=f"partner_sample_id in {cohort_ids}",
    sample_sets=None,
)

gt = allel.GenotypeArray(ds_haps["call_genotype"].data)

ac = gt.count_alleles()
biallelic_mask = ac.is_biallelic()

print(f"removing biallelics for {cohort} {contig}, {biallelic_mask.sum()} snps, array is shape: {ds_haps['call_genotype'].shape}")
ds_haps = ds_haps.isel(variants=biallelic_mask)
print(f" for {cohort} {contig} after, array is shape: {ds_haps['call_genotype'].shape}")
pos = ds_haps["variant_position"].values
ht = gt[biallelic_mask, :].to_haplotypes()
metadata.set_index('sample_id', inplace = True)
hap_ids = np.repeat(metadata.loc[ds_haps.sample_id,'partner_sample_id'].values, 2)

alive_dead_dict = {}
alive_dead_dict['alive'] = random.query(f"phenotype == 'alive' & population == @cohort")['specimen']
alive_dead_dict['dead']  = random.query(f"phenotype == 'dead' & population == @cohort")['specimen']
for pheno in ['alive', 'dead']:

    if Path(f"H12_outputs/H12_{cohort}_{pheno}.{contig}.tsv").exists():
        print(f"skipping {cohort} {pheno} {contig}")
        continue

    subsample_ht = ht[:, np.isin(hap_ids, alive_dead_dict[pheno])]
    print(f"--------- Running H12 on {cohort} {pheno} | Chromosome {contig} ----------")
    x, h12 = h12_gwss(
        subsample_ht,
        pos,
        window_size=1000,
    )

    h12_df = pd.DataFrame({'midpoint':x, 'h12':h12})
    h12_df.to_csv(f"H12_outputs/H12_{cohort}.{pheno}.{contig}.tsv", sep="\t")

# And now the randomisations
for i in np.arange(1, num_randomisations+1):
    i = '%0.4d' % i

    alive_dead_dict['alive'] = random.query(f"r{i} == 'alive' & population == @cohort")['specimen']
    alive_dead_dict['dead']  = random.query(f"r{i} == 'dead' & population == @cohort")['specimen']

    for pheno in ['alive', 'dead']:

        if Path(f"H12_outputs/H12_{cohort}_{pheno}{i}.{contig}.tsv").exists():
            print(f"skipping {cohort} {pheno} {i} {contig}")
            continue

        print(f"--------- Running H12 on {cohort} {pheno} randomisation {i} | Chromosome {contig} ----------")
        subsample_ht = ht[:, np.isin(hap_ids, alive_dead_dict[pheno])]
        x, h12 = h12_gwss(
            subsample_ht,
            pos,
            window_size=1000,
        )

        h12_df = pd.DataFrame({'midpoint':x, 'h12':h12})
        h12_df.to_csv(f"H12_outputs/H12_{cohort}_{pheno}{str(i)}.{contig}.tsv", sep="\t")

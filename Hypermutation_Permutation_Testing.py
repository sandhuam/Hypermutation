import pandas as pd
import numpy as np
impact_data = pd.read_csv(
    '/Users/amarsandhu/Documents/hypermutation_pr/data_mutations_extended_20190403_prepped.txt', sep='\t'
)

# Separate silent from non-silent
svars = 'Silent'
nsvars = ['Frame_Shift_Del', 'Frame_Shift_Ins', 'In_Frame_Del', 'In_Frame_Ins', 'Missense_Mutation',
          'Nonsense_Mutation', 'Nonstop_Mutation', 'Splice_Site', 'Translation_Start_Site']
silent_data = impact_data[impact_data["Variant_Classification"] == svars]
nonsilent_data = impact_data[impact_data["Variant_Classification"].isin(nsvars)]
nonsilent_data.reset_index(drop=True, inplace=True)

# Find hypermutated samples
hyperm = nonsilent_data.copy()
hyperm.loc[:, "tmz_signature":"msi_high"] = hyperm.loc[:, "tmz_signature":"msi_high"].fillna(value=False)
# hyperm = hyperm[hyperm.loc[:, "tmz_signature"]]
# hyperm = hyperm[hyperm.loc[:, "pol_signature"]]
hyperm = hyperm[hyperm.loc[:, "mmr_signature"]]
hyperm = hyperm[hyperm.loc[:, "msi_high"]]
hyperm.to_csv('/Users/amarsandhu/Documents/hypermutation_pr/hypermutation_samples.csv')

# Load gene lengths
gene_lengths = pd.read_csv('/Users/amarsandhu/Documents/hypermutation_pr/impact_gene_cds_lengths.tsv', sep='\t')

# Find mutation rate
hyper_gene_lengths = gene_lengths[gene_lengths["Hugo_Symbol"].isin(hyperm["Hugo_Symbol"])]
total_length = hyper_gene_lengths["cds_length"].sum()
samples = hyperm["Tumor_Sample_Barcode"].value_counts()
# ultra_mutation_rates = sample_mutation_rates[sample_mutation_rates > 100]
# sample_mutation_rates = samples/total_length * 1000000
# rate_mb = sample_mutation_rates.mean()
rate = samples.mean()

# Find mutational signatures and proportions of trinucleotide substitutions
# Done in R
hypermutation_tnm = pd.read_csv('/Users/amarsandhu/Documents/hypermutation_pr/hypermutation_trinucleotides.csv')
hypermutation_tnm.set_index("Unnamed: 0", inplace=True)
trinucleotide_sum = hypermutation_tnm.sum(axis=0)
trinucleotide_rate = trinucleotide_sum/sum(trinucleotide_sum)
substitutions = trinucleotide_rate.index

# Drop the mutations - choosing genes selected against
sample_gene_frequency = pd.crosstab(index=hyperm["Tumor_Sample_Barcode"], columns=hyperm["Hugo_Symbol"])
mean_gene = sample_gene_frequency.mean()
test_table = pd.DataFrame(0, index=gene_lengths["Hugo_Symbol"], columns=["mean_gene", "gene_mutation_freq", "counter"])
test_table.loc[mean_gene.index, "mean_gene"] = mean_gene
test_table["counter"] = 0
for i in range(100000):
    # mutation_order = np.random.choice(substitutions, size=int(round(rate)), p=trinucleotide_rate)
    cum_lengths = np.cumsum(gene_lengths["cds_length"])
    rnd = np.random.randint(1, cum_lengths.iloc[-1], int(round(rate)))
    selected_genes = np.searchsorted(cum_lengths, rnd, side='right')
    mutated_genes = gene_lengths.iloc[selected_genes]
    gene_mutation_freq = mutated_genes["Hugo_Symbol"].value_counts()
    test_table["gene_mutation_freq"] = 0
    test_table.loc[gene_mutation_freq.index, "gene_mutation_freq"] = gene_mutation_freq
    temp_count = np.where(test_table["mean_gene"] <= test_table["gene_mutation_freq"], 1, 0)
    test_table["counter"] += temp_count

print(test_table[test_table.counter >= 95000])

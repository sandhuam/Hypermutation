import pandas as pd
import numpy as np
import pyreadr

impact_data = pd.read_csv(
    '/Users/amarsandhu/Documents/hypermutation_pr/data_mutations_extended_20190403_prepped.txt', sep='\t')

# Separate silent from non-silent
svars = 'Silent'
nsvars = ['Frame_Shift_Del', 'Frame_Shift_Ins', 'In_Frame_Del', 'In_Frame_Ins', 'Missense_Mutation',
          'Nonsense_Mutation', 'Nonstop_Mutation', 'Splice_Site', 'Translation_Start_Site']
silent_data = impact_data[impact_data["Variant_Classification"] == svars]
nonsilent_data = impact_data[impact_data["Variant_Classification"].isin(nsvars)]

# Find hypermutated samples
hyperm = nonsilent_data.copy()
hyperm.loc[:, "tmz_signature":"msi_high"] = hyperm.loc[:, "tmz_signature":"msi_high"].fillna(value=False)
hyperm = hyperm[hyperm.loc[:, "mmr_signature"]]
hyperm.drop(134752, inplace=True)
hyperm.to_csv('/Users/amarsandhu/Documents/hypermutation_pr/hypermutation_samples.csv')

# Load gene lengths
gene_lengths = pd.read_csv('/Users/amarsandhu/Documents/hypermutation_pr/impact_gene_cds_lengths.tsv', sep='\t')

# Find mutation rate
samples = hyperm["Tumor_Sample_Barcode"].value_counts()
rate = samples.mean()

# Find average number of mutations per gene per sample and initialize final table
sample_gene_frequency = pd.crosstab(index=hyperm["Tumor_Sample_Barcode"], columns=hyperm["Hugo_Symbol"])
mean_gene = sample_gene_frequency.mean()
test_table = pd.DataFrame(0, index=gene_lengths["Hugo_Symbol"], columns=["mean_gene", "gene_mutation_freq", "counter"])
test_table.loc[mean_gene.index, "mean_gene"] = mean_gene

# Find proportions of trinucleotide substitutions - substitution matrix made in R
hypermutation_tnm = pd.read_csv('/Users/amarsandhu/Documents/hypermutation_pr/hypermutation_trinucleotides.csv')
hypermutation_tnm.set_index("Unnamed: 0", inplace=True)
trinucleotide_sum = hypermutation_tnm.sum(axis=0)
trinucleotide_rate = trinucleotide_sum/sum(trinucleotide_sum)
substitutions = trinucleotide_rate.index

# Load trinucleotide sequences
read_file = pyreadr.read_r('/Users/amarsandhu/Documents/hypermutation_pr/impact_genes_coding_sequences.rds')
cseq = read_file[None]
cs = cseq.set_index("Hugo_Symbol")
cs = cs.reindex(gene_lengths["Hugo_Symbol"])
sequences = cs["coding"]


# Drop mutations on the IMPACT genes based on available trinucleotides
def gene_freqs(substitution):
    good_chars = substitution.translate({ord('['): None, ord(']'): None, ord('>'): None})
    original = good_chars[:2] + good_chars[-1]
    count = sequences.apply(lambda x: x.count(original))
    count /= sum(count)
    return count


def mutated_genes(subst_list):
    freq_list = subst_list.apply(gene_freqs)
    mut_genes = pd.Series(0)
    for i in range(int(round(rate))):
        mut_genes[i] = np.random.choice(gene_lengths.Hugo_Symbol, p=freq_list.iloc[i])
    return mut_genes


for i in range(100):
    mutation_order = pd.Series(np.random.choice(substitutions, size=int(round(rate)), p=trinucleotide_rate))
    mgenes = mutated_genes(mutation_order)
    gene_mutation_freq = mgenes.value_counts()
    test_table["gene_mutation_freq"] = 0
    test_table.loc[gene_mutation_freq.index, "gene_mutation_freq"] = gene_mutation_freq
    test_table["counter"] += np.where(test_table["mean_gene"] <= test_table["gene_mutation_freq"], 1, 0)

print(test_table[test_table.counter >= 95])

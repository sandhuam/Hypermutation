import pandas as pd
import numpy as np
import pyreadr
import regex as re
import matplotlib.pyplot as plt
import seaborn as sns
import math
import statistics
import rpy2.robjects as ro
import subprocess

impact_data = pd.read_csv('/Users/amarsandhu/Documents/hypermutation_pr/all_impact_mutations_annotated_cohort.maf',
                          sep='\t')  # Load impact maf

genes = pd.Series(pd.read_csv('/Users/amarsandhu/Documents/hypermutation_pr/impact341_gene_panel2.txt',
                              sep='\t').columns)  # Load list of genes in 341 gene panel
filtered_data = impact_data[impact_data["Hugo_Symbol"].isin(genes)]  # Only look at 341 genes
filtered_data = filtered_data[filtered_data["is-a-hotspot"].isna()]  # Ignore Hotspots

# Load Noah's cluster results
cluster_results = pd.read_csv("/Users/amarsandhu/Documents/hypermutation_pr/hypermutators_clustering_results.txt",
                              sep='\t')

# Find hypermutated samples
hyperm = filtered_data.copy()  # Insert how you would separate the hypermutated samples
# hyperm = hyperm[hyperm.PATIENT_ID.isin(cluster_results.patient_id)][hyperm[hyperm.PATIENT_ID.isin(
#     cluster_results.patient_id)].tmb >= 10]  # Trying to filter by cluster results (filters by tmb >= 10)

# The hypermutated samples should be split into cohorts
hyperm_tmz_glioma = hyperm
hyperm_pol_endometrial = hyperm
hyperm_msi_colorectal = hyperm

# Load trinucleotide sequences - 341 genes
read_file = pyreadr.read_r('/Users/amarsandhu/Documents/hypermutation_pr/impact_genes_coding_sequences.rds')
cseq = read_file[None]
cs = cseq.set_index("Hugo_Symbol")
sequences = cs["coding"]

# Generate Series of trinucleotide substitution occurences (executes slowly)
trinucleotide_matrix = hyperm[hyperm.Variant_Type == "SNP"].apply(lambda x: (x.Ref_Tri[0] + "[" + x.Reference_Allele + ">" +
                                                      x.Tumor_Seq_Allele2 + "]" + x.Ref_Tri[2]), axis=1).value_counts()

# Load gene lengths
gene_lengths = pd.read_csv('/Users/amarsandhu/Documents/hypermutation_pr/impact_gene_cds_lengths.tsv', sep='\t')

# Load Tumor Suppressors
tumor_suppressors = list(pd.read_csv('/Users/amarsandhu/Documents/hypermutation_pr/tumor_suppressors.txt'))

# Load Gene Covariates
gene_covariates = pd.read_csv("/Users/amarsandhu/Documents/hypermutation_pr/gene_covariates.txt", sep='\t')


def perm_test(hyperm, genes, gene_covariates, trinucleotide_matrix=None, sequences=None,
              gene_lengths=None, tumor_suppressors=None, multiple_mutations=True, mutational_signature=False,
              gc_content_norm=True, reptime_norm=True):
    genes = genes[genes.isin(hyperm.Hugo_Symbol)]

    # Find number of mutations per gene and initialize final table
    sum_gene = hyperm[hyperm.Variant_Type == "SNP"].Hugo_Symbol.value_counts()
    test_table = pd.DataFrame(0, index=genes, columns=["sum_gene", "gene_mutation_freq", "counter1",
                                                       "counter2", "tumor_suppressor", "sum_oncogene"]
                              )
    test_table.loc[sum_gene.index, "sum_gene"] = sum_gene

    if gc_content_norm or reptime_norm:  # Load genomic correlates
        gene_covariates = gene_covariates.set_index("Hugo_Symbol")
        gene_covariates = gene_covariates.reindex(genes)
        test_table["reptime"] = gene_covariates.reptime
        test_table["percentage_gene_gc_content"] = gene_covariates.percentage_gene_gc_content

        k = round(1 + 3.322 * math.log(genes.__len__()))

        if gc_content_norm:
            def gc_norm(table):
                gc_table = table[["gene_mutation_freq", "percentage_gene_gc_content"]]
                gc_table = gc_table.sort_values(by="percentage_gene_gc_content")
                gc_table.gene_mutation_freq[
                    gc_table.gene_mutation_freq == 0] = 1  # Add pseudocount of 1 to get rid of 0s
                gc_cont = np.array_split(gc_table, k)
                gc_c = [gc_cont[i].gene_mutation_freq * statistics.median(gc_table.gene_mutation_freq) /
                        statistics.median(gc_cont[i].gene_mutation_freq) for i in range(gc_cont.__len__())]
                gc_corr_freq = pd.concat(gc_c)
                gc_corr_freq = gc_corr_freq.reindex(table.index)
                return gc_corr_freq

        if reptime_norm:
            def rep_norm(table):
                rep_table = table[["gene_mutation_freq", "reptime"]]
                rep_table = rep_table.sort_values(by="reptime")
                rep_table.gene_mutation_freq[
                    rep_table.gene_mutation_freq == 0] = 1  # Add pseudocount of 1 to get rid of 0s
                rep_cont = np.array_split(rep_table, k)
                rep_c = [rep_cont[i].gene_mutation_freq * statistics.median(rep_table.gene_mutation_freq) /
                         statistics.median(rep_cont[i].gene_mutation_freq) for i in range(rep_cont.__len__())]
                rep_corr_freq = pd.concat(rep_c)
                rep_corr_freq = rep_corr_freq.reindex(table.index)
                return rep_corr_freq

    if mutational_signature:  # Trinucleotides
        if trinucleotide_matrix is None:
            # Generate matrix of trinucleotide substitutions
            trinucleotide_sum = hyperm[hyperm.Variant_Type == "SNP"].apply(
                lambda x: (x.Ref_Tri[0] + "[" + x.Reference_Allele + ">" +
                           x.Tumor_Seq_Allele2 + "]" + x.Ref_Tri[2]), axis=1).value_counts()
        else:
            trinucleotide_sum = trinucleotide_matrix
        substitutions = trinucleotide_sum.index

        assert sequences is not None, "No Sequence Data Passed"
        # Organize Sequences
        sequences = sequences.reindex(genes)

        # Calculate probability of hitting a gene for each trinucleotide substitution
        def gene_freqs(substitution):
            good_chars = substitution.translate({ord('['): None, ord(']'): None, ord('>'): None})
            original = good_chars[:2] + good_chars[-1]
            count = sequences.apply(lambda x: len(re.findall(original, x, overlapped=True)))
            count /= sum(count)
            return count

        # Initialize Gene Probability Table
        trin_genes = pd.Series(substitutions).apply(gene_freqs).set_index(substitutions)

        # Run simulation
        for i in range(1000):
            test_table["gene_mutation_freq"] = 0
            mgenes = trin_genes.apply(
                lambda x: pd.Series(np.random.choice(sequences.index, size=trinucleotide_sum[x.name],
                                                     p=x)).value_counts(), axis=1)
            gene_mutation_freq = mgenes.sum()
            test_table.loc[gene_mutation_freq.index, "gene_mutation_freq"] = gene_mutation_freq
            if gc_content_norm: test_table.gene_mutation_freq = gc_norm(test_table)  # GC-Content Normalization
            if reptime_norm: test_table.gene_mutation_freq = rep_norm(test_table)  # Replication Timing Normalization
            test_table["counter1"] += np.where(test_table["sum_gene"] <= test_table["gene_mutation_freq"], 1, 0)
            test_table["counter2"] += np.where(test_table["sum_gene"] >= test_table["gene_mutation_freq"], 1, 0)

    else:  # Gene Lengths
        assert gene_lengths is not None, "Gene Lengths not Passed as Argument"

        # Load gene lengths
        gene_lengths = gene_lengths.set_index("Hugo_Symbol")
        gene_lengths = gene_lengths.reindex(genes)

        # Run simulation
        for i in range(10000):
            test_table["gene_mutation_freq"] = 0
            mgenes = pd.Series(np.random.choice(gene_lengths.index, size=sum_gene.sum(),
                                                p=gene_lengths.cds_length/sum(gene_lengths.cds_length)))
            gene_mutation_freq = mgenes.value_counts()
            if gc_content_norm: test_table.gene_mutation_freq = gc_norm(test_table)  # GC-Content Normalization
            if reptime_norm: test_table.gene_mutation_freq = rep_norm(test_table)  # Replication Timing Normalization
            test_table.loc[gene_mutation_freq.index, "gene_mutation_freq"] = gene_mutation_freq
            test_table["counter1"] += np.where(test_table["sum_gene"] <= test_table["gene_mutation_freq"], 1, 0)
            test_table["counter2"] += np.where(test_table["sum_gene"] >= test_table["gene_mutation_freq"], 1, 0)

    # Annotate Tumor Suppressors and find oncogenic mutation frequency/gene
    if tumor_suppressors is not None:
        test_table.tumor_suppressor = test_table.index.isin(tumor_suppressors)
    sum_oncogene = hyperm[hyperm.Variant_Type == "SNP"][
        ~hyperm[hyperm.Variant_Type == "SNP"].oncogenic.isna()].Hugo_Symbol.value_counts()
    test_table.loc[sum_oncogene.index, "sum_oncogene"] = sum_oncogene

    if multiple_mutations:
        # Check for double mutations in the data
        sample_gene = pd.crosstab(index=hyperm[hyperm.Variant_Type == "SNP"].Tumor_Sample_Barcode,
                                  columns=hyperm[hyperm.Variant_Type == "SNP"].Hugo_Symbol)
        double_muts = (sample_gene >= 2).sum()

        # Initialize table
        test_table2 = pd.DataFrame(0, index=genes,
                                   columns=["double_muts", "double_sim", "counter1", "counter2",
                                            "tumor_suppressor"])
        test_table2.loc[double_muts.index, "double_muts"] = double_muts

        if mutational_signature:  # Trinucleotides
            # Find trinucleotide substitutions per sample
            trinucleotide_rate = trinucleotide_sum / trinucleotide_sum.sum()
            sig_genes = trin_genes.multiply(trinucleotide_rate, axis="index").sum()
            test_table2 = test_table2.reindex(sig_genes.index)
            test_table2 = test_table2.fillna(0)

            # Run Simulation
            for i in range(1000):
                sim_sample_gene = sample_gene.sum(axis=1).apply(
                    lambda x: pd.Series(np.random.choice(sig_genes.index, size=x, p=sig_genes)).value_counts())
                double_sim = (sim_sample_gene >= 2).sum()
                test_table2.loc[double_sim.index, "double_sim"] = double_sim
                test_table2["counter1"] += np.where(test_table2["double_muts"] <= test_table2["double_sim"], 1, 0)
                test_table2["counter2"] += np.where(test_table2["double_muts"] >= test_table2["double_sim"], 1, 0)

        else:  # Gene Length
            for i in range(1000):
                sim_sample_gene = sample_gene.sum(axis=1).apply(
                    lambda x: pd.Series(np.random.choice(gene_lengths.index, size=x, p=gene_lengths.cds_length/sum(
                        gene_lengths.cds_length))).value_counts())
                double_sim = (sim_sample_gene >= 2).sum()
                test_table2.loc[double_sim.index, "double_sim"] = double_sim
                test_table2["counter1"] += np.where(test_table2["double_muts"] <= test_table2["double_sim"], 1, 0)
                test_table2["counter2"] += np.where(test_table2["double_muts"] >= test_table2["double_sim"], 1, 0)

        # Annotate Tumor Suppressors
        if tumor_suppressors is not None:
            test_table2.tumor_suppressor = test_table2.index.isin(tumor_suppressors)
    else:
        test_table2 = []
    return test_table, test_table2


tuple_tmz_glioma = perm_test(hyperm_tmz_glioma, genes, gene_covariates, trinucleotide_matrix, sequences, gene_lengths,
                           tumor_suppressors, mutational_signature=True, gc_content_norm=False, reptime_norm=False,
                           multiple_mutations=False)

tuple_pol_endometrial = perm_test(hyperm_pol_endometrial, genes, gene_covariates, trinucleotide_matrix, sequences,
                                  gene_lengths, tumor_suppressors, mutational_signature=True, gc_content_norm=False,
                                  reptime_norm=False, multiple_mutations=False)

tuple_msi_colorectal = perm_test(hyperm_msi_colorectal, genes, gene_covariates, trinucleotide_matrix, sequences,
                           gene_lengths, tumor_suppressors, mutational_signature=True, gc_content_norm=False,
                           reptime_norm=False, multiple_mutations=False)

df_tmz_glioma = tuple_tmz_glioma[0].reset_index()
df_pol_endometrial = tuple_pol_endometrial[0].reset_index()
df_msi_colorectal = tuple_msi_colorectal[0].reset_index()

df_tmz_glioma["cohort"] = "TMZ Glioma"
df_pol_endometrial["cohort"] = "POLE Endometrial Cancer"
df_msi_colorectal["cohort"] = "MSI-High Colorectal Cancer"

returned_dataframe = df_tmz_glioma.append([df_pol_endometrial, df_msi_colorectal], ignore_index=True)

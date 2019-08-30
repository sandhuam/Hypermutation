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

impact_data = pd.read_csv(
    '/Users/amarsandhu/Documents/hypermutation_pr/data_mutations_extended_20190403_prepped.txt', sep='\t')

tcga_data = pd.read_csv('/Users/amarsandhu/Documents/hypermutation_pr/mc3.v0.2.8.PUBLIC.maf', sep='\t')

# Filter data
genes = pd.Series(
    pd.read_csv('/Users/amarsandhu/Documents/hypermutation_pr/impact341_gene_panel2.txt', sep='\t').columns)
filtered_data = impact_data[impact_data["Hugo_Symbol"].isin(genes)]
filtered_data = filtered_data[filtered_data["is-a-hotspot"].isna()]

filtered_tcga = tcga_data.copy()
# Initialize Silent Mutations
# unfiltered_data = pd.read_csv("/Users/amarsandhu/Documents/hypermutation_pr/unfilteredMafWithTriunc.txt", sep='\t')
# silent_data = unfiltered_data[unfiltered_data.Variant_Classification == "Silent"]
# silent_data = silent_data[silent_data.Hugo_Symbol.isin(genes)]
# silent_data = silent_data[silent_data.Tumor_Sample_Barcode.isin(
#     silent_data.Tumor_Sample_Barcode.value_counts()[silent_data.Tumor_Sample_Barcode.value_counts() >= 20].index)]

# Initialize Clustering Results
cluster_results = pd.read_csv("/Users/amarsandhu/Documents/hypermutation_pr/hypermutators_clustering_results.txt",
                              sep='\t')
# Find hypermutated samples
hyperm = filtered_tcga.copy()
# hyperm = hyperm[hyperm.PATIENT_ID.isin(cluster_results.patient_id)][hyperm[hyperm.PATIENT_ID.isin(
#     cluster_results.patient_id)].tmb >= 10]
hyperm.loc[:, "tmz_signature":"msi_high"] = hyperm.loc[:, "tmz_signature":"msi_high"].fillna(value=False)
# hyperm = silent_data.copy()
hyperm.to_csv('/Users/amarsandhu/Documents/hypermutation_pr/hypermutation_samples.csv')
filepath = '/Users/amarsandhu/Documents/hypermutation_pr/hypermutation_samples.csv'

# Load trinucleotide sequences
read_file = pyreadr.read_r('/Users/amarsandhu/Documents/hypermutation_pr/impact_genes_coding_sequences.rds')
cseq = read_file[None]
cs = cseq.set_index("Hugo_Symbol")
sequences = cs["coding"]

# Load Matrix of Subsitutions
subprocess.call(['Rscript', '/Users/amarsandhu/Documents/R_Documents/trinucleotide_matrix.R'], shell=False)
trinucleotide_matrix = pd.read_csv('/Users/amarsandhu/Documents/hypermutation_pr/hypermutation_trinucleotides.csv')

# Load Gene Lengths
gene_lengths = pd.read_csv('/Users/amarsandhu/Documents/hypermutation_pr/gencode.v19.basic.exons.bed', sep='\t', header=None)
gene_lengths[3] = gene_lengths[3].str.split(pat=':', n=1, expand=True)[0]
gene_lengths.rename(columns={3: 'Hugo_Symbol'}, inplace=True)
gene_lengths = gene_lengths.drop(columns=[4, 5])
gene_lengths['cds_length'] = gene_lengths[2] - gene_lengths[1] + 1
gene_lengths = gene_lengths.groupby("Hugo_Symbol").sum()
gene_lengths.reset_index(inplace=True)

# Load Tumor Suppressors
tumor_suppressors = list(pd.read_csv('/Users/amarsandhu/Documents/hypermutation_pr/tumor_suppressors.txt'))

# Load Covariates
gene_covariates = pd.read_csv("/Users/amarsandhu/Documents/hypermutation_pr/gene_covariates.txt", sep='\t')


def perm_test(hyperm, genes, gene_covariates, trinucleotide_matrix=None, sequences=None, filepath=None,
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
            assert filepath is not None, "No Trinucleotide Matrix Passed and No Filepath Passed"
            # Generate matrix of trinucleotide substitutions
            rstring = """
                function(filepath){
                    overall <- read.csv(filepath, sep = ",", stringsAsFactors = F)
                    library("maftools")
                    overallMaf = read.maf(overall, vc_nonSyn = 'Silent')
                    overall_tnm = trinucleotideMatrix(overallMaf, ref_genome = "/Users/amarsandhu/Documents/R_Documents/b37.fasta")
                    write.csv(overall_tnm$nmf_matrix, '/Users/amarsandhu/Documents/hypermutation_pr/hypermutation_trinucleotides.csv')
                }
            """
            rfunc = ro.r(rstring)
            rfunc(filepath)
            tnm = pd.read_csv('/Users/amarsandhu/Documents/hypermutation_pr/hypermutation_trinucleotides.csv')
            trinucleotide_sum = tnm.set_index("Unnamed: 0").sum(axis=0)
        else:
            trinucleotide_sum = trinucleotide_matrix.set_index("Unnamed: 0").sum(axis=0)
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


result = perm_test(hyperm, hyper_filepath, gene_list, gene_covars, mutational_signature=True)

# Plot results
# Clonality
sns.set_style("whitegrid")
ax = sns.barplot(["Overall", "TMZ Glioma", "POLE Endometrial", "MSI Colorectal"],
                 [hyperm.clonal.sum() / len(hyperm.clonal),
                  hyperm_tmz_glioma.clonal.sum() / len(hyperm_tmz_glioma.clonal),
                  hyperm_pol_endometrial.clonal.sum() / len(hyperm_pol_endometrial.clonal),
                  hyperm_msi_colorectal.clonal.sum() / len(hyperm_msi_colorectal.clonal)])
ax.text(0, hyperm.clonal.sum() / len(hyperm.clonal), round(hyperm.clonal.sum() / len(hyperm.clonal), 2), ha="center")
ax.text(1, hyperm_tmz_glioma.clonal.sum() / len(hyperm_tmz_glioma.clonal),
        round(hyperm_tmz_glioma.clonal.sum() / len(hyperm_tmz_glioma.clonal), 2), ha="center")
ax.text(2, hyperm_pol_endometrial.clonal.sum() / len(hyperm_pol_endometrial.clonal),
        round(hyperm_pol_endometrial.clonal.sum() / len(hyperm_pol_endometrial.clonal), 2), ha="center")
ax.text(3, hyperm_msi_colorectal.clonal.sum() / len(hyperm_msi_colorectal.clonal),
        round(hyperm_msi_colorectal.clonal.sum() / len(hyperm_msi_colorectal.clonal), 2), ha="center")
ax.set(xlabel="Fraction of Mutations that are Clonal", ylabel="Subset", title="Clonality in Hypermutation")

# Hotspots
hots = pd.crosstab(impact_data["is-a-hotspot"].fillna(value='N'), impact_data.tmb)
plt.scatter(hots.columns, hots.iloc[1] / (hots.iloc[0] + hots.iloc[1]))
plt.xlabel("Tumor Mutational Burden (nmut/Mb)")
plt.ylabel("Fraction that are Hotspot Mutations")
plt.title("Proportion of Mutations that are Hotspot Mutations by Tumor Mutational Burden")

# Genomic covariates
plt.figure(1)
plt.subplots_adjust(hspace=0.4)
plt.subplot(222)
t1 = test_table.copy()
t1 = t1.sort_values(by="counter1")
polyt1 = np.poly1d(np.polyfit(t1.counter1 / 1000, t1.reptime, deg=2))
t2 = test_table.copy()
t2 = t2.sort_values(by="counter2", ascending=False)
polyt2 = np.poly1d(np.polyfit(-t2.counter2 / 1000, t2.reptime, deg=2))
plt.plot(t1.counter1 / 1000, polyt1(t1.counter1 / 1000))
plt.scatter(t1.counter1 / 1000, t1.reptime, s=5)
plt.xlabel("Deviation from Null - Enrichment")
plt.ylabel("Replication Timing")
plt.title("Deviation from Null vs. Replication Timing")

plt.subplot(221)
plt.plot(-t2.counter2 / 1000, polyt2(-t2.counter2 / 1000))
plt.scatter(-t2.counter2 / 1000, t2.reptime, s=5)
plt.xlabel("Deviation from Null - Depletion (Absolute Value = p-value)")
plt.ylabel("Replication Timing")
plt.title("Deviation from Null vs. Replication Timing")

plt.subplot(224)
t1 = test_table.copy()
t1 = t1.sort_values(by="counter1")
polyt1 = np.poly1d(np.polyfit(t1.counter1 / 1000, t1.percentage_gene_gc_content, deg=2))
t2 = test_table.copy()
t2 = t2.sort_values(by="counter2", ascending=False)
polyt2 = np.poly1d(np.polyfit(-t2.counter2 / 1000, t2.percentage_gene_gc_content, deg=2))
plt.plot(t1.counter1 / 1000, polyt1(t1.counter1 / 1000))
plt.scatter(t1.counter1 / 1000, t1.percentage_gene_gc_content, s=5)
plt.xlabel("Deviation from Null - Enrichment p-value")
plt.ylabel("GC-Content")
plt.title("Deviation from Null vs. GC-Content")

plt.subplot(223)
plt.plot(-t2.counter2 / 1000, polyt2(-t2.counter2 / 1000))
plt.scatter(-t2.counter2 / 1000, t2.percentage_gene_gc_content, s=5)
plt.xlabel("Deviation from Null - Depletion (Absolute Value = p-value)")
plt.ylabel("GC-Content")
plt.title("Deviation from Null vs. GC-Content")

# Selection by oncogenic mutation frequency
plt.figure(2)
plt.subplots_adjust(hspace=0.4)
y = test_table_msi_colorectal[test_table_msi_colorectal.tumor_suppressor]
y = y.sort_values(by="sum_oncogene").dropna()
# plt.subplot(211)
plt.bar(y.index, y.counter1 / 1000, color='b')
plt.axhline(y=0.05, color='r')
plt.text(0, 0.05, '0.05', va='bottom', color='r')
plt.xticks(range(len(y)), y.index, rotation=60, fontsize=6)
# plt.axis(ymin=0, ymax=0.15)
plt.xlabel("Gene (In Order of Increasing Mutation Frequency)")
plt.ylabel("p-value")
plt.title("Selection for Tumor Suppressors by Oncogenic Mutation Frequency in MSI-High Colorectal Cancer "
          "n={}".format(len(y)))
plt.show()

polytmzglioma = np.poly1d(np.polyfit(test_table_tmz_glioma[test_table_tmz_glioma.tumor_suppressor].sum_oncogene,
                                     test_table_tmz_glioma[test_table_tmz_glioma.tumor_suppressor].counter1 / 1000,
                                     deg=1))
polypolendometrial = np.poly1d(
    np.polyfit(test_table_pol_endometrial[test_table_pol_endometrial.tumor_suppressor].sum_oncogene,
               test_table_pol_endometrial[test_table_pol_endometrial.tumor_suppressor].counter1 / 1000, deg=1))
polymsicolorectal = np.poly1d(
    np.polyfit(test_table_msi_colorectal[test_table_msi_colorectal.tumor_suppressor].sum_oncogene,
               test_table_msi_colorectal[test_table_msi_colorectal.tumor_suppressor].counter1 / 1000, deg=1))

plt.plot(range(16), polytmzglioma(range(16)), label="TMZ Glioma")
plt.plot(range(53), polypolendometrial(range(53)), label="POLE Endometrial Cancer")  # 56, 56, 57
plt.plot(range(104), polymsicolorectal(range(104)), label="MSI-High Colorectal Cancer")  # 120
plt.axhline(y=0.05, color='m')
plt.xlabel("Number of Oncogenic Mutations")
plt.ylabel("p-value")
plt.title("Selection for Tumor Suppressors vs Oncogenic Mutation Frequency")
plt.axis(ymin=0, xmax=30)
plt.legend()

h = test_table[test_table.tumor_suppressor]
h = h.sort_values(by="sum_oncogene").dropna()
plt.subplot(212)
plt.bar(h.index, h.counter1 / 1000, color='b')
plt.axhline(y=0.05, color='r')
plt.text(0, 0.05, '0.05', va='bottom', color='r')
plt.xticks(range(len(h)), h.index, rotation=60, fontsize=6)
# plt.axis(ymin=0, ymax=0.15)
plt.xlabel("Gene (In Order of Increasing Mutation Frequency)")
plt.ylabel("p-value")
plt.title("Selection Against Tumor Suppressors by Oncogenic Mutation Frequency in TMZ Glioma  n={}".format(len(h)))
plt.show()

# Normal Selection
plt.figure(3)
plt.subplots_adjust(hspace=0.4)
plt.subplot(211)
y = result[0].counter1.copy()
y = y.sort_values()
mask1 = y / 1000 <= 0.05
mask2 = y / 1000 > 0.05
plt.xticks(range(len(y)), y.index, rotation=60, fontsize=8)
plt.axhline(y=0.05, color='r')
plt.xlabel("Gene")
plt.ylabel("p-value")
plt.title("Selection For Genes with Silent Mutations n={}".format(len(y[mask1])))
plt.bar(y[mask1].index, y[mask1] / 1000, color='m')
plt.bar(y[mask2].index, y[mask2] / 1000, color='b')

plt.subplot(212)
h = result[0].counter2.copy()
h = h.sort_values()
mask1 = h / 1000 <= 0.05
mask2 = h / 1000 > 0.05
plt.xticks(range(len(h)), h.index, rotation=60, fontsize=6)
plt.axhline(y=0.05, color='r')
plt.xlabel("Gene")
plt.ylabel("p-value")
plt.title("Selection Against Genes with Silent Mutations n={}".format(len(h[mask1])))
plt.bar(h[mask1].index, h[mask1] / 1000, color='m')
plt.bar(h[mask2].index, h[mask2] / 1000, color='b')

# Gene Commonality
plt.figure(4)
gene_commonality = pd.read_csv("/Users/amarsandhu/Documents/hypermutation_pr/gene_commonality.txt", sep='\t')
y = test_table_msi_colorectal[test_table_msi_colorectal.tumor_suppressor]
y = y.reindex(gene_commonality.colorectal).dropna()
cts = pd.Series(range(11)).apply(lambda x: round(x * len(y) / 10)).diff().dropna()
cts = [int(i) for i in cts]
bands = sum([[s] * n for s, n in zip(range(1, 11), cts)], [])
y["bands"] = bands
y = y.sort_values(by="counter1")
y.counter1 = y.counter1 / 1000
color_dict = {}
for i in range(1, 11):
    color_dict[i] = sns.color_palette("RdBu_r", 10).as_hex()[i - 1]
y["colors"] = y.bands.apply(lambda x: color_dict[x])
plt.bar(y.index, y.counter1, color=y.colors)
plt.xticks(range(len(y)), y.index, rotation=60, fontsize=6)
plt.xlabel("Gene")
plt.ylabel("p-value")
plt.legend(["Top 10% Most Related"])
plt.title("Tumor Suppressors Selected for in TMZ Glioma")

y = y.sort_values(by="bands")
y.counter1 = y.counter1 <= 0.05
props = [y.counter1[round(len(y) / 10 * i):round(len(y) / 10 * (i + 1))].sum() /
         y.counter1[round(len(y) / 10 * i):round(len(y) / 10 * (i + 1))].__len__() for i in range(10)]
plt.bar(["{}th-{}th".format(i * 10, (i + 1) * 10) for i in range(10)], props)
plt.xlabel("Percentile of Most Common Tumor Suppressors")
plt.ylabel("Proportion of Genes with p-value <= 0.05 for Positive Selection")
sns.set_style("whitegrid")

# Multiple Mutations
plt.figure(5)
plt.subplots_adjust(hspace=0.4)
y = test_table2[test_table2.tumor_suppressor].counter1.sort_values()
y = y[y / 1000 <= 0.15]
mask1 = y / 1000 > 0.05
mask2 = y / 1000 <= 0.05
plt.subplot(211)
plt.bar(y[mask2].index, y[mask2] / 1000, color='m')
plt.bar(y[mask1].index, y[mask1] / 1000, color='b')
plt.axhline(y=0.05, color='r')
plt.text(0, 0.05, '0.05', va='bottom', color='r')
plt.xticks(range(len(y)), y.index, rotation=60, fontsize=8)
plt.axis(ymin=0, ymax=0.15)
plt.xlabel("Gene")
plt.ylabel("p-value")
plt.title(
    "Tumor Suppressors With More Multiple Oncogenic Mutations in MSI-High Colorectal Cancer n={}".format(len(y[mask2])))
plt.show()

h = test_table2.counter2.sort_values()
h = h[h / 1000 <= 0.15]
mask1 = h / 1000 > 0.05
mask2 = h / 1000 <= 0.05
plt.subplot(212)
plt.bar(h[mask2].index, h[mask2] / 1000, color='m')
plt.bar(h[mask1].index, h[mask1] / 1000, color='b')
plt.axhline(y=0.05, color='r')
plt.text(0, 0.05, '0.05', va='bottom', color='r')
plt.xticks(range(len(h)), h.index, rotation=60, fontsize=8)
plt.axis(ymin=0, ymax=0.15)
plt.xlabel("Gene")
plt.ylabel("p-value")
plt.title("Genes With Fewer Multiple Oncogenic Mutations in MSI-High Colorectal Cancer n={}".format(len(h[mask2])))
plt.show()

# Extra Plot
plt.figure(6)
y = hyperm.copy()
ycounts = y.Hugo_Symbol.value_counts()
ycounts = ycounts.reindex(gene_list.isin(hyperm.Hugo_Symbol))
yvalue = sigg / sum(sigg) * ycounts.sum()
x = np.linspace(0, 55, 10)
plt.scatter(rr, yvalue)
plt.plot(x, x, '-r', label='y=x')
plt.xlabel("Corrected Number of Mutations per Gene")
plt.ylabel(
    "Probability of Hitting a Gene (Mutational Signature & Available Trinucleotides)*Number of Mutations in the Data")
plt.legend()

y = hyperm.copy()
ycounts = y.Hugo_Symbol.value_counts()
ycounts = ycounts.reindex(gene_lengths.drop(index=472).Hugo_Symbol)
yvalue = gene_lengths.drop(index=472).cds_length / sum(gene_lengths.drop(index=472).cds_length) * ycounts.sum()

import statsmodels.formula.api as smf

gene_covariates = gene_covariates.reindex(ycounts.index)
yy = pd.DataFrame(0, index=ycounts.index, columns=["counts", "percentage_gene_gc_content", "reptime"])
yy["counts"] = ycounts.apply(lambda x: math.log10(x))
yy[["percentage_gene_gc_content", "reptime"]] = gene_covariates[["percentage_gene_gc_content", "reptime"]]
reg = smf.ols("counts ~ percentage_gene_gc_content + reptime", data=yy).fit()
rr = yy.counts - reg.params[1] * yy.percentage_gene_gc_content - reg.params[2] * yy.reptime + reg.params[0]

plt.figure(7)
x = np.linspace(0, 2.5, 10)
plt.scatter(rr, yvalue.apply(math.log10))
plt.plot(x, x, '-r', label='y=x')
plt.xlabel("Corrected Number of Mutations per Gene")
plt.ylabel("Gene Length/Sum of Gene Lengths*Number of Mutations in the Data")
plt.legend()

plt.figure(8)
plt.subplot(211)
plt.plot(gene_covariates.percentage_gene_gc_content.sort_values())
plt.axhline(y=gene_covariates.percentage_gene_gc_content["CARD11"], color="r")
plt.title("GC-Content")
plt.subplot(212)
plt.plot(gene_covariates.reptime.sort_values())
plt.axhline(y=gene_covariates.reptime["CARD11"], color='r')
plt.title("Replication Timing")

import statsmodels.formula.api as smf

gene_covariates = gene_covariates.reindex(ycounts.index)
yy = pd.DataFrame(0, index=ycounts.index, columns=["counts", "percentage_gene_gc_content", "reptime"])
yy["counts"] = ycounts.apply(lambda x: math.log10(x))
yy[["percentage_gene_gc_content", "reptime"]] = gene_covariates[["percentage_gene_gc_content", "reptime"]]
reg = smf.ols("counts ~ percentage_gene_gc_content + reptime", data=yy).fit()
rr = 10 ** (yy.counts - reg.params[1] * yy.percentage_gene_gc_content - reg.params[2] * yy.reptime)
yy.reptime.hist()

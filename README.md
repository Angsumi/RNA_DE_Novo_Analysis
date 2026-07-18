# De Novo RNA-Seq Transcriptome Assembly and Downstream Analysis

This repository contains the complete processed results, scientific figures, parsed CSV tables, and reproducible scripts from the de novo transcriptome assembly and differential expression workflow executed on the GCE virtual machine `rna-denovo-vm`.

---

## 📂 Repository Structure

The data in this repository is structured as follows:

```
RNA_DE_Novo_Analysis/
├── README.md                              # This detailed workflow documentation
├── TrinityStats.txt                      # Contig length distribution and N50 statistics
├── de_samples.txt                        # edgeR sample relationship matrix
├── trinotate_DE_annotation_report.xls     # Full Excel sheet of DE genes functional annotation
│
├── figures/                              # Publication-ready plots (.png)
│   ├── busco_completeness.png            # Assembly ortholog completeness stacked bar chart
│   ├── de_overlap_venn.png               # 3-way Venn diagram of DE gene overlaps
│   ├── de_heatmap.png                    # Hierarchical clustering heatmap of 432 DE genes
│   ├── pca_plot.png                      # Principal Component Analysis (PCA) of samples
│   ├── candidate_genes_expression.png    # Expression profiles of top 6 candidate genes
│   ├── top_go_categories.png             # Top 15 most frequent Gene Ontology (GO) terms
│   └── top_pfam_domains.png              # Top 15 most frequent Pfam domain designations
│
├── csv_data/                             # Comma-Separated Values (CSV) of all figures & results
│   ├── AZ1_vs_AZ2_gene_DE_results.csv     # Pairwise gene differential expression values (logFC, FDR)
│   ├── AZ1_vs_AZ3_gene_DE_results.csv     #
│   ├── AZ2_vs_AZ3_gene_DE_results.csv     #
│   ├── AZ1_vs_AZ2_isoform_DE_results.csv  # Pairwise transcript-level DE values (logFC, FDR)
│   ├── AZ1_vs_AZ3_isoform_DE_results.csv  #
│   ├── AZ2_vs_AZ3_isoform_DE_results.csv  #
│   ├── busco_completeness_percentages.csv # Raw percentages/counts for BUSCO chart
│   ├── de_overlap_counts.csv              # Exact overlaps counts for the Venn diagram
│   ├── de_genes_expression_matrix.csv     # Clean matrix of TMM expression values for the 432 DE genes
│   ├── isoform_switching_candidates.csv   # List of 228 genes showing isoform switching
│   ├── pca_coordinates.csv                # PC1 and PC2 coordinates of samples
│   ├── top_candidate_genes_expression.csv # Expressions of top 6 candidate genes
│   ├── top_go_categories.csv              # Counts/descriptions for the GO term bar chart
│   └── top_pfam_domains.csv               # Counts/descriptions for the Pfam domains chart
│
└── Scripts_and_Logs/                     # Automation scripts and run console logs
    ├── run_trinity_norm_12cpu_fix.sh     # Trinity de novo assembler run script
    ├── run_salmon_manual.sh              # Salmon abundance quantification automation script
    ├── run_de_analysis_no_replicates.sh  # edgeR differential expression automation script
    ├── run_trinotate_subset.sh           # Patched Trinotate subset annotation pipeline script
    ├── extract_de_sequences.py           # Custom Python script to identify and extract DE sequences
    ├── generate_figures.py               # Python visualization pipeline script
    ├── generate_csv_data.py              # Python script compiling TSV/CSV raw tables
    ├── additional_analyses.py            # Python script executing PCA, TFs, candidate expression
    ├── downstream_run.log                # Master downstream analysis execution log
    ├── de_analysis.log                   # edgeR execution console log
    ├── trinotate_subset.log              # Database building & annotation compilation log
    └── infernal.log                      # cmscan Rfam non-coding RNA search log
```

---

## 🧬 Summary of Completed Analyses

### 1. De Novo Assembly (Trinity)
* **Total Genes**: 2,682,805
* **Total Transcripts**: 2,944,151
* **Contig N50**: 363 bp
* **Total Assembled Bases**: 1.07 Gbp
* **Assembly Completeness (BUSCO)**: 68.0% Complete (85/125 orthologs), 28.8% Fragmented (36/125), 3.2% Missing (4/125).

### 2. Abundance Estimation (Salmon)
* Abundance quantification mapped raw reads back to the assembled contigs:
  * **AZ1**: 80.37% alignment rate (10.75M mapped reads)
  * **AZ2**: 79.47% alignment rate (10.78M mapped reads)
  * **AZ3**: 80.28% alignment rate (6.73M mapped reads)

### 3. Differential Expression (edgeR)
* Pairwise DE genes identified at FDR < 0.05 (dispersion fixed at 0.1):
  * **AZ1 vs AZ2**: 993 genes (441 with |logFC| > 2)
  * **AZ1 vs AZ3**: 1,103 genes (417 with |logFC| > 2)
  * **AZ2 vs AZ3**: 839 genes (453 with |logFC| > 2)
* A union of **432 unique differentially expressed genes** (mapping to 2,840 transcripts) was extracted for downstream annotation.

### 4. Functional Annotation (Trinotate)
* Functional profiling of the 2,840 DE transcripts:
  * Open Reading Frames (ORFs) predicted by **TransDecoder**.
  * Homologous proteins annotated by **DIAMOND** blastp (against SwissProt).
  * Functional domains annotated by **HMMER** (against Pfam).
  * Non-coding RNAs scanned by **Infernal (cmscan)** (against Rfam).
* Compiled into the comprehensive report sheet `trinotate_DE_annotation_report.xls`.

---

## 🛠 How to Use the Raw CSV Data
All files in `csv_data/` are standard Comma-Separated Values (CSV). You can open them directly in:
* **Microsoft Excel** or **Google Sheets** for sorting and custom filtering.
* **GraphPad Prism** or **OriginPro** for plotting.
* **R** (`read.csv()`) or **Python** (`pandas.read_csv()`) for advanced statistics and custom clustering.
